"""
Test AntiFold_runner with mocked AntiFold model.

TDD approach: test the interface and workflow before real model integration.
"""

import importlib.util
import os
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest


class MockArgs:
    """Minimal args object for testing."""
    def __init__(self, pdbdir="", outquiver="", outpdbdir="outputs"):
        self.pdbdir = pdbdir
        self.quiver = ""
        self.outquiver = outquiver
        self.outpdbdir = outpdbdir
        self.runlist = ""
        self.checkpoint_name = "check.point"
        self.debug = False
        self.deterministic = False
        self.loop_string = "H1,H2,H3,L1,L2,L3"
        self.seqs_per_struct = 2
        self.temperature = 0.2
        self.augment_eps = 0
        self.protein_features = "full"
        self.omit_AAs = "CX"
        self.num_connections = 48
        self.allow_x = False


def _load_sequence_design_script_module():
    script_path = Path(__file__).resolve().parents[2] / "scripts" / "proteinmpnn_interface_design.py"
    spec = importlib.util.spec_from_file_location("proteinmpnn_interface_design", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestAntiFold_runnerWorkflow:
    """Test AntiFold_runner workflow with mocked model."""

    def test_runner_loads_converter(self):
        """Runner should instantiate HLTtoAntiFoldConverter."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager

        args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test")
        struct_manager = StructManager(args)
        runner = AntiFold_runner(args, struct_manager)
        assert runner.converter is not None

    def test_runner_sets_device(self):
        """Runner should detect GPU/CPU and set device."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager

        args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test")
        struct_manager = StructManager(args)
        runner = AntiFold_runner(args, struct_manager)
        assert runner.device in ("cpu", "cuda:0")

    def test_sequence_optimize_returns_list_of_tuples(self):
        """sequence_optimize must return list of (seq, score) tuples."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager
        from rfantibody.proteinmpnn.sample_features import SampleFeatures
        from rfantibody.util.pose import Pose

        args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test")
        struct_manager = StructManager(args)
        runner = AntiFold_runner(args, struct_manager)

        # Mock the actual AntiFold call
        runner._run_antifold = MagicMock(return_value=[
            ("QVQLQESGPGLVKPSETLSLTCAVSGYSISSGYYWGWIRQPPGKGLEWIGSIYHSGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAGLTQSSHNDANWGQGTLVTVSS/" +
             "VLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKRLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLNPVFGGGTKLEIKR", 0.5),
        ])

        # Create minimal pose
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        pose = Pose.from_pdb(input_pdb)
        sample_feats = SampleFeatures(pose, "test_tag")
        sample_feats.loop_string2fixed_res("H1,H2,H3,L1,L2,L3")

        results = runner.sequence_optimize(sample_feats)

        assert isinstance(results, list)
        assert len(results) > 0
        for seq, score in results:
            assert isinstance(seq, str)
            assert isinstance(score, float)
            assert len(seq) > 0

    def test_run_model_creates_output(self):
        """run_model should process one structure and write output."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager

        with tempfile.TemporaryDirectory() as tmpdir:
            args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test", outpdbdir=tmpdir)
            struct_manager = StructManager(args)
            runner = AntiFold_runner(args, struct_manager)

            # Mock sequence generation
            runner._run_antifold = MagicMock(return_value=[
                ("QVQLQESGPGLVKPSETLSLTCAVSGYSISSGYYWGWIRQPPGKGLEWIGSIYHSGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAGLTQSSHNDANWGQGTLVTVSS/" +
                 "VLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKRLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLNPVFGGGTKLEIKR", 0.5),
            ])

            # load_pose expects full path, not just tag
            pdb_path = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
            runner.run_model(pdb_path, args)

            # Should have created output file(s)
            output_files = list(Path(tmpdir).glob("*.pdb"))
            assert len(output_files) > 0, f"No output PDB files in {tmpdir}"


class TestSequenceDesignLoop:
    def test_zero_persisted_outputs_exit_nonzero(self, tmp_path: Path):
        module = _load_sequence_design_script_module()

        class FakeStructManager:
            output_pdb = True
            output_quiver = False

            def __init__(self):
                self.outpdbdir = str(tmp_path / "outputs")
                self.outquiver = SimpleNamespace(get_tags=lambda: [])
                self.checkpoints: list[str] = []

            def iterate(self):
                yield "design_a"
                yield "design_b"

            def record_checkpoint(self, pdb):
                self.checkpoints.append(pdb)

        class FakeRunner:
            def run_model(self, _pdb, _args):
                return None

        struct_manager = FakeStructManager()
        args = SimpleNamespace(debug=False, backend="antifold")

        with pytest.raises(SystemExit, match="produced zero persisted outputs"):
            module._run_design_loop(args, struct_manager, FakeRunner())

        assert struct_manager.checkpoints == []

    def test_only_successful_outputs_are_checkpointed(self, tmp_path: Path):
        module = _load_sequence_design_script_module()

        class FakeStructManager:
            output_pdb = True
            output_quiver = False

            def __init__(self):
                self.outpdbdir = str(tmp_path / "outputs")
                self.outquiver = SimpleNamespace(get_tags=lambda: [])
                self.checkpoints: list[str] = []

            def iterate(self):
                yield "design_success"
                yield "design_failure"

            def record_checkpoint(self, pdb):
                self.checkpoints.append(pdb)

        class FakeRunner:
            def run_model(self, pdb, _args):
                output_dir = tmp_path / "outputs"
                output_dir.mkdir(parents=True, exist_ok=True)
                if pdb == "design_success":
                    (output_dir / "design_success_dldesign_0.pdb").write_text("ATOM\n", encoding="utf-8")
                    return None
                raise RuntimeError("boom")

        struct_manager = FakeStructManager()
        args = SimpleNamespace(debug=False, backend="antifold")

        module._run_design_loop(args, struct_manager, FakeRunner())

        assert struct_manager.checkpoints == ["design_success"]

    def test_overwriting_existing_output_still_counts_as_success(self, tmp_path: Path):
        module = _load_sequence_design_script_module()

        class FakeStructManager:
            output_pdb = True
            output_quiver = False

            def __init__(self):
                self.outpdbdir = str(tmp_path / "outputs")
                self.outquiver = SimpleNamespace(get_tags=lambda: [])
                self.checkpoints: list[str] = []

            def iterate(self):
                yield "design_success"

            def record_checkpoint(self, pdb):
                self.checkpoints.append(pdb)

        class FakeRunner:
            def run_model(self, _pdb, _args):
                output_dir = tmp_path / "outputs"
                output_dir.mkdir(parents=True, exist_ok=True)
                (output_dir / "design_success_dldesign_0.pdb").write_text(
                    "UPDATED_OUTPUT_WITH_NEW_SIZE\n",
                    encoding="utf-8",
                )

        output_dir = tmp_path / "outputs"
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "design_success_dldesign_0.pdb").write_text("OLD\n", encoding="utf-8")

        struct_manager = FakeStructManager()
        args = SimpleNamespace(debug=False, backend="antifold")

        module._run_design_loop(args, struct_manager, FakeRunner())

        assert struct_manager.checkpoints == ["design_success"]
