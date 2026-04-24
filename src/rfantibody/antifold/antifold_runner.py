"""
AntiFold model runner compatible with RFantibody's ProteinMPNN interface.

Wraps AntiFold inverse folding model for antibody CDR sequence design.
"""

import os
import tempfile

import torch

from rfantibody.antifold.imgt_converter import HLTtoAntiFoldConverter
from rfantibody.proteinmpnn.util_protein_mpnn import aa_1_3


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

        # TODO: load actual AntiFold model weights when available
        self.model = None
        self.model_path = getattr(args, "antifold_checkpoint_path", "")

    def _run_antifold(self, pdb_path: str, regions: list[str]) -> list[tuple[str, float]]:
        """Call AntiFold model to generate sequences.

        This is a placeholder that should be replaced with actual
        AntiFold inference code once the model is installed.

        Args:
            pdb_path: Path to IMGT-numbered PDB input
            regions: List of AntiFold region names to design

        Returns:
            List of (sequence_str, score_float) tuples
        """
        # Placeholder: return original sequence with dummy score
        # Real implementation will call AntiFold's model.sample()
        raise NotImplementedError(
            "AntiFold model not loaded. "
            "Install AntiFold and implement this method."
        )

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
