"""
Test suite for AntiFold integration.

TDD approach: write tests first, then implement.
"""

import os
import tempfile
from pathlib import Path

import pytest

from rfantibody.antifold.imgt_converter import HLTtoAntiFoldConverter


class TestHLTtoAntiFoldConversion:
    """Test converting HLT-format PDB to AntiFold-compatible input."""

    def test_extract_hl_chains_removes_target(self):
        """HLT PDB should have T chain removed, keeping only H and L."""
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        converter = HLTtoAntiFoldConverter()
        output_lines = converter.convert(input_pdb, heavy_chain="H", light_chain="L")

        # Verify no T chain atoms in output
        for line in output_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21:22].strip()
                assert chain != "T", f"Target chain T found in output: {line[:30]}"

    def test_preserves_h_and_l_chains(self):
        """Both H and L chains must be present in output."""
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        converter = HLTtoAntiFoldConverter()
        output_lines = converter.convert(input_pdb, heavy_chain="H", light_chain="L")

        chains = set()
        for line in output_lines:
            if line.startswith("ATOM"):
                chains.add(line[21:22].strip())

        assert "H" in chains, "Heavy chain H missing from output"
        assert "L" in chains, "Light chain L missing from output"

    def test_cdrs_aligned_to_imgt_anchors(self):
        """CDR residues must land in AntiFold's IMGT ranges (27-38/56-65/105-117)."""
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        converter = HLTtoAntiFoldConverter()
        output_lines = converter.convert(input_pdb, heavy_chain="H", light_chain="L")

        cdr_old_labels: dict[tuple[str, int], str] = {}
        with open(input_pdb) as f:
            for line in f:
                if line.startswith("REMARK PDBinfo-LABEL"):
                    parts = line.strip().split()
                    if len(parts) >= 4 and parts[-1][:1] in ("H", "L") and parts[-1][1:] in ("1", "2", "3"):
                        cdr_old_labels[(parts[-1][0], int(parts[2]))] = parts[-1]

        original_ca, converted_ca = [], []
        with open(input_pdb) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    original_ca.append((line[21:22].strip(), int(line[22:26].strip())))
        for line in output_lines:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                converted_ca.append((line[21:22].strip(), int(line[22:26].strip())))

        imgt_ranges = {"1": range(27, 39), "2": range(56, 66), "3": range(105, 118)}
        for (chain_orig, resno_orig), label in cdr_old_labels.items():
            if (chain_orig, resno_orig) not in original_ca:
                continue
            idx = original_ca.index((chain_orig, resno_orig))
            _, new_resno = converted_ca[idx]
            assert new_resno in imgt_ranges[label[1]], (
                f"CDR{label} old-resno {resno_orig} mapped to {new_resno}, "
                f"outside IMGT range {imgt_ranges[label[1]].start}-{imgt_ranges[label[1]].stop - 1}"
            )

    def test_cdr_regions_mapped_correctly(self):
        """CDR regions from HLT REMARKs should map to AntiFold IMGT regions."""
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        converter = HLTtoAntiFoldConverter()
        regions = converter.infer_antifold_regions(input_pdb)

        # Should return at least CDR regions
        assert isinstance(regions, list)
        assert len(regions) > 0
        # Valid AntiFold region names
        valid_regions = {"CDR1", "CDR2", "CDR3", "CDRH1", "CDRH2", "CDRH3",
                         "CDRL1", "CDRL2", "CDRL3", "FW1", "FW2", "FW3", "FW4",
                         "FWH1", "FWH2", "FWH3", "FWH4", "FWL1", "FWL2", "FWL3", "FWL4"}
        for r in regions:
            assert r in valid_regions, f"Invalid region name: {r}"

    def test_output_is_valid_pdb_format(self):
        """Output lines should be valid PDB ATOM/TER records."""
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        converter = HLTtoAntiFoldConverter()
        output_lines = converter.convert(input_pdb, heavy_chain="H", light_chain="L")

        for line in output_lines:
            if line.startswith("ATOM"):
                # Standard PDB ATOM line should be at least 54 chars
                assert len(line) >= 54, f"ATOM line too short: {repr(line)}"
                # Record name
                assert line[:4] == "ATOM"
                # Chain ID at position 21
                assert line[21:22].strip() in ("H", "L")


class TestAntiFold_runnerInterface:
    """Test AntiFold_runner has the same interface as ProteinMPNN_runner."""

    def test_runner_has_sequence_optimize_method(self):
        """AntiFold_runner must implement sequence_optimize()."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        assert hasattr(AntiFold_runner, "sequence_optimize")

    def test_runner_has_proteinmpnn_compatible_signature(self):
        """sequence_optimize must accept SampleFeatures and return list of (seq, score)."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        import inspect

        sig = inspect.signature(AntiFold_runner.sequence_optimize)
        params = list(sig.parameters.keys())
        assert "sample_feats" in params, "sequence_optimize must accept sample_feats parameter"

    def test_runner_has_run_model_method(self):
        """AntiFold_runner must implement run_model() for StructManager integration."""
        from rfantibody.antifold.antifold_runner import AntiFold_runner
        assert hasattr(AntiFold_runner, "run_model")
