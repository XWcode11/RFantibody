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

    def test_residue_numbering_starts_at_one(self):
        """Each chain should start numbering from 1 for IMGT compatibility."""
        input_pdb = "test/proteinmpnn/inputs_for_test/ab_des_0.pdb"
        converter = HLTtoAntiFoldConverter()
        output_lines = converter.convert(input_pdb, heavy_chain="H", light_chain="L")

        h_residues = []
        l_residues = []
        for line in output_lines:
            if line.startswith("ATOM"):
                chain = line[21:22].strip()
                resno = int(line[22:26].strip())
                if chain == "H":
                    h_residues.append(resno)
                elif chain == "L":
                    l_residues.append(resno)

        if h_residues:
            assert min(h_residues) == 1, f"H chain numbering should start at 1, got {min(h_residues)}"
        if l_residues:
            assert min(l_residues) == 1, f"L chain numbering should start at 1, got {min(l_residues)}"

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


class TestAntiFoldRunnerInterface:
    """Test AntiFoldRunner has the same interface as ProteinMPNN_runner."""

    def test_runner_has_sequence_optimize_method(self):
        """AntiFoldRunner must implement sequence_optimize()."""
        from rfantibody.antifold.antifold_runner import AntiFoldRunner
        assert hasattr(AntiFoldRunner, "sequence_optimize")

    def test_runner_has_proteinmpnn_compatible_signature(self):
        """sequence_optimize must accept SampleFeatures and return list of (seq, score)."""
        from rfantibody.antifold.antifold_runner import AntiFoldRunner
        import inspect

        sig = inspect.signature(AntiFoldRunner.sequence_optimize)
        params = list(sig.parameters.keys())
        assert "sample_feats" in params, "sequence_optimize must accept sample_feats parameter"

    def test_runner_has_run_model_method(self):
        """AntiFoldRunner must implement run_model() for StructManager integration."""
        from rfantibody.antifold.antifold_runner import AntiFoldRunner
        assert hasattr(AntiFoldRunner, "run_model")
