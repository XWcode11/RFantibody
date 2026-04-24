"""
Integration test for real AntiFold inference.

Runs actual model inference (requires GPU + weights).
Skipped if weights not found or GPU unavailable.
"""

import os
import tempfile

import pytest
import torch


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


@pytest.mark.skipif(
    not torch.cuda.is_available(),
    reason="Requires GPU for AntiFold inference"
)
@pytest.mark.skipif(
    not os.path.exists("weights/antifold_model.pt"),
    reason="Requires AntiFold model weights"
)
class TestAntiFoldRealInference:
    """Test real AntiFold model inference end-to-end."""

    def test_model_loads_successfully(self):
        """AntiFold model should load without errors."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager

        args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test")
        struct_manager = StructManager(args)
        runner = AntiFold_runner(args, struct_manager)

        # Trigger lazy loading
        model = runner.model
        assert model is not None

    def test_run_antifold_returns_sequences(self):
        """_run_antifold should return real sequences with scores."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager
        from rfantibody.antifold.imgt_converter import HLTtoAntiFoldConverter

        args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test")
        struct_manager = StructManager(args)
        runner = AntiFold_runner(args, struct_manager)

        # Convert test PDB to AntiFold format
        converter = HLTtoAntiFoldConverter()
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.writelines(converter.convert(input_pdb, "H", "L"))
            tmp_pdb = f.name

        try:
            regions = converter.infer_antifold_regions(input_pdb)
            results = runner._run_antifold(tmp_pdb, regions)

            assert isinstance(results, list)
            assert len(results) == args.seqs_per_struct
            for seq, score in results:
                assert isinstance(seq, str)
                assert "/" in seq  # AntiFold FASTA format
                assert isinstance(score, float)
                assert score != 0.0
        finally:
            os.unlink(tmp_pdb)

    def test_sequence_optimize_threads_to_pose(self):
        """sequence_optimize should produce threadable sequences."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        from rfantibody.proteinmpnn.struct_manager import StructManager
        from rfantibody.proteinmpnn.sample_features import SampleFeatures
        from rfantibody.util.pose import Pose

        args = MockArgs(pdbdir="test/proteinmpnn/inputs_for_test")
        struct_manager = StructManager(args)
        runner = AntiFold_runner(args, struct_manager)

        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        pose = Pose.from_pdb(input_pdb)
        sample_feats = SampleFeatures(pose, "test_tag")
        sample_feats.loop_string2fixed_res("H1,H2,H3,L1,L2,L3")

        results = runner.sequence_optimize(sample_feats)

        assert isinstance(results, list)
        assert len(results) > 0
        for seq, score in results:
            assert isinstance(seq, str)
            assert "/" not in seq  # Should be cleaned (no chain separator)
            assert len(seq) > 0
            assert isinstance(score, float)
