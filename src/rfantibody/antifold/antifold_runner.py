"""
AntiFold model runner compatible with RFantibody's ProteinMPNN interface.

Wraps AntiFold inverse folding model for antibody CDR sequence design.
"""

import os
import re
import sys
import tempfile

import torch

from rfantibody.antifold.imgt_converter import HLTtoAntiFoldConverter
from rfantibody.config import PathConfig


class AntiFoldRunner:
    """AntiFold sequence design runner with ProteinMPNN-compatible interface."""

    def __init__(self, args, struct_manager):
        """Initialize AntiFold model.

        Args:
            args: Argument namespace with model parameters
            struct_manager: StructManager for I/O handling
        """
        self.struct_manager = struct_manager
        self.converter = HLTtoAntiFoldConverter()

        if torch.cuda.is_available():
            print("Found GPU, will run AntiFold on GPU")
            self.device = "cuda:0"
        else:
            print("No GPU found, running AntiFold on CPU")
            self.device = "cpu"

        self.temperature = getattr(args, "temperature", 0.2)
        self.seqs_per_struct = getattr(args, "seqs_per_struct", 1)
        self.loop_string = getattr(args, "loop_string", "H1,H2,H3,L1,L2,L3")

        # Lazy-load model on first use
        self._model = None

    @property
    def model(self):
        """Lazy-load AntiFold model on first access."""
        if self._model is None:
            self._model = self._load_model()
        return self._model

    def _load_model(self):
        """Load AntiFold model with weights.

        Returns:
            Loaded AntiFold ESM-IF1 model
        """
        # Ensure antifold package is importable
        antifold_repo = "/tmp/antifold_repo"
        if os.path.exists(antifold_repo) and antifold_repo not in sys.path:
            sys.path.insert(0, antifold_repo)

        from antifold.antiscripts import load_IF1_checkpoint
        from antifold.esm.pretrained import _load_IF1_local

        weight_path = PathConfig.get_weight_path("antifold")
        if not weight_path.exists():
            raise FileNotFoundError(
                f"AntiFold weights not found at {weight_path}. "
                "Please download antifold_model.pt to the weights directory."
            )

        print(f"Loading AntiFold model from {weight_path}...")
        model, _ = _load_IF1_local()
        model = load_IF1_checkpoint(model, str(weight_path))
        model = model.eval()
        model = model.to(self.device)
        print(f"AntiFold model loaded on {self.device}")
        return model

    def _parse_score(self, description: str) -> float:
        """Parse score from FASTA description string.

        Args:
            description: SeqRecord description from AntiFold

        Returns:
            Parsed score (global_score if available, else score)
        """
        # Try to extract global_score first, fallback to score
        match = re.search(r"global_score=([\d.]+)", description)
        if match:
            return float(match.group(1))
        match = re.search(r"score=([\d.]+)", description)
        if match:
            return float(match.group(1))
        return 0.0

    def _run_antifold(self, pdb_path: str, regions: list[str]) -> list[tuple[str, float]]:
        """Call AntiFold model to generate sequences.

        Args:
            pdb_path: Path to IMGT-numbered PDB input
            regions: List of AntiFold region names to design

        Returns:
            List of (sequence_str, score_float) tuples
        """
        import pandas as pd
        from antifold.antiscripts import (
            get_pdbs_logits,
            sample_from_df_logits_HL,
            seed_everything,
        )

        pdb_dir = os.path.dirname(pdb_path)
        pdb_basename = os.path.splitext(os.path.basename(pdb_path))[0]

        # Build DataFrame for AntiFold
        df = pd.DataFrame({
            "pdb": [pdb_basename],
            "Hchain": ["H"],
            "Lchain": ["L"],
        })

        # Set seed for reproducibility
        seed = getattr(self, "seed", 42)
        seed_everything(seed)

        # Get logits from model
        df_logits_list = get_pdbs_logits(
            model=self.model,
            pdbs_csv_or_dataframe=df,
            pdb_dir=pdb_dir,
            batch_size=1,
            seed=seed,
        )

        if not df_logits_list:
            raise RuntimeError("AntiFold returned no logits")

        df_logits = df_logits_list[0]

        # Sample sequences
        fasta_dict = sample_from_df_logits_HL(
            df_logits,
            sample_n=self.seqs_per_struct,
            sampling_temp=self.temperature,
            regions_to_mutate=regions,
            seed=seed,
            verbose=False,
        )

        # Parse results: skip original sequence, keep sampled ones
        sequences = []
        for seq_id, seq_record in fasta_dict.items():
            # Original sequence has key like "{name}_HL", skip it
            if "__" not in seq_id:
                continue

            seq = str(seq_record.seq)
            score = self._parse_score(seq_record.description)
            sequences.append((seq, score))

        return sequences

    def sequence_optimize(self, sample_feats) -> list[tuple[str, float]]:
        """Run AntiFold sequence optimization on a pose.

        Args:
            sample_feats: SampleFeatures object with pose and design info

        Returns:
            List of (sequence_str, score_float) tuples
        """
        # Save pose to temporary PDB
        temp_pdb = os.path.join(tempfile.gettempdir(), f"{sample_feats.tag}_antifold_input.pdb")
        sample_feats.pose.dump_pdb(temp_pdb)

        try:
            # Convert HLT to AntiFold-compatible format (H/L only, renumbered)
            antifold_lines = self.converter.convert(
                temp_pdb,
                heavy_chain="H",
                light_chain="L",
            )

            # Write converted PDB to temp file
            converted_pdb = os.path.join(
                tempfile.gettempdir(),
                f"{sample_feats.tag}_antifold_converted.pdb",
            )
            with open(converted_pdb, "w") as f:
                f.writelines(antifold_lines)

            # Infer design regions from loop_string
            regions = self.converter.infer_antifold_regions(temp_pdb)

            # Call AntiFold (or mock in tests)
            sequences = self._run_antifold(converted_pdb, regions)

            # Clean up AntiFold FASTA format (Hseq/Lseq -> continuous sequence)
            cleaned_sequences = []
            for seq, score in sequences:
                # Remove chain separator '/' if present
                cleaned_seq = seq.replace("/", "")
                cleaned_sequences.append((cleaned_seq, score))

            return cleaned_sequences

        finally:
            # Cleanup temp files
            if os.path.exists(temp_pdb):
                os.remove(temp_pdb)
            converted_pdb = os.path.join(
                tempfile.gettempdir(),
                f"{sample_feats.tag}_antifold_converted.pdb",
            )
            if os.path.exists(converted_pdb):
                os.remove(converted_pdb)

    def run_model(self, tag, args):
        """Run full AntiFold design workflow for one structure.

        Args:
            tag: Structure identifier
            args: Arguments with design parameters
        """
        import time

        t0 = time.time()
        print(f"Attempting pose: {tag}")

        # Load the pose
        pose = self.struct_manager.load_pose(tag)

        # Initialize the features
        from rfantibody.proteinmpnn.sample_features import SampleFeatures

        sample_feats = SampleFeatures(pose, tag)

        # Parse the loop string and determine which residues should be designed
        sample_feats.loop_string2fixed_res(args.loop_string)

        # Run AntiFold sequence optimization
        seqs_scores = self.sequence_optimize(sample_feats)

        # Iterate through each seq score pair and thread the sequence onto the pose
        prefix = f"{sample_feats.tag}_dldesign"
        for idx, (seq, _) in enumerate(seqs_scores):
            sample_feats.thread_mpnn_seq(seq)

            outtag = f"{prefix}_{idx}"
            self.struct_manager.dump_pose(sample_feats.pose, outtag)

        seconds = int(time.time() - t0)
        print(f"Struct: {tag} reported success in {seconds} seconds")
