import argparse
import os
import random
import sys
import time

import numpy as np
import torch

import rfantibody.proteinmpnn.util_protein_mpnn as mpnn_util
from rfantibody.proteinmpnn.sample_features import SampleFeatures
from rfantibody.proteinmpnn.struct_manager import StructManager

#################################
# Parse Arguments
#################################

parser = argparse.ArgumentParser()

# I/O Arguments
parser.add_argument("-pdbdir", type=str, default="", help='The name of a directory of pdbs to run through the model')
parser.add_argument("-quiver", type=str, default="", help='The name of a quiver file to run this metric on.')

parser.add_argument("-outquiver", type=str, default="",
                    help="The name of the quiver file to which output structs will be written")
parser.add_argument("-outpdbdir", type=str, default="outputs",
                    help='The directory to which the output PDB files will be written')
parser.add_argument("-runlist", type=str, default='',
                    help="The path of a list of pdb tags to run (default: ''; Run all PDBs")
parser.add_argument("-checkpoint_name", type=str, default='check.point',
                    help="The name of a file where tags which have finished will be written (default: check.point)")
parser.add_argument("-debug", action="store_true", default=False,
                    help='When active, errors will cause the script to crash and the error message ' + \
                         'to be printed out (default: False)')
parser.add_argument("-deterministic", action="store_true", default=False,
                    help='Enable deterministic mode by setting fixed seeds and deterministic algorithms (default: False)')

# Design Arguments
parser.add_argument("-loop_string", type=str, default='H1,H2,H3,L1,L2,L3',
                    help='The list of loops which you wish to design')
parser.add_argument("-seqs_per_struct", type=int, default="1",
                    help="The number of sequences to generate for each structure (default: 1)")

# ProteinMPNN Specific Arguments
default_ckpt = os.path.join( os.path.dirname(__file__), '/home/weights/ProteinMPNN_v48_noise_0.2.pt')
parser.add_argument("-checkpoint_path", type=str, default=default_ckpt)
parser.add_argument("-temperature", type=float, default=0.000001, help='An a3m file containing the MSA of your target')
parser.add_argument("-augment_eps", type=float, default=0,
                    help='The variance of random noise to add to the atomic coordinates (default 0)')
parser.add_argument("-protein_features", type=str, default='full',
                    help='What type of protein features to input to ProteinMPNN (default: full)')
parser.add_argument("-omit_AAs", type=str, default='CX',
                    help='A string of all residue types (one letter case-insensitive) that you would not like to ' + \
                         'use for design. Letters not corresponding to residue types will be ignored')
parser.add_argument("-num_connections", type=int, default=48,
                    help='Number of neighbors each residue is connected to, default 48, higher number leads to ' + \
                         'better interface design but will cost more to run the model.')
parser.add_argument("-allow_x", action="store_true", default=False,
                    help='Allow X (unknown) residues in output. Useful for debugging to see exactly which positions were designed.')
parser.add_argument("-backend", type=str, default="proteinmpnn", choices=["proteinmpnn", "antifold"],
                    help='Sequence design backend to use (default: proteinmpnn)')

def _configure_determinism(args) -> None:
    if not args.deterministic:
        return

    print("Setting up deterministic mode with fixed seeds")
    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)
    torch.cuda.manual_seed_all(42)

    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def _normalize_input_tag(tag: str) -> str:
    """Normalize an input identifier to the basename used by output tags."""
    return os.path.basename(tag).split('.')[0]


def _snapshot_output_state(struct_manager, tag: str):
    """Capture the persisted outputs currently associated with one input tag."""
    prefix = f"{_normalize_input_tag(tag)}_dldesign_"

    if struct_manager.output_pdb:
        if not os.path.isdir(struct_manager.outpdbdir):
            return {}

        snapshot = {}
        for name in os.listdir(struct_manager.outpdbdir):
            if not (name.startswith(prefix) and name.endswith('.pdb')):
                continue

            path = os.path.join(struct_manager.outpdbdir, name)
            if not os.path.isfile(path):
                continue

            stat_result = os.stat(path)
            snapshot[name] = (stat_result.st_mtime_ns, stat_result.st_size)

        return snapshot

    if struct_manager.output_quiver:
        return tuple(sorted(tag_name for tag_name in struct_manager.outquiver.get_tags() if tag_name.startswith(prefix)))

    return ()


def _persisted_outputs_changed(before_snapshot, after_snapshot) -> bool:
    """Return whether a backend call created or updated persisted outputs."""
    return before_snapshot != after_snapshot and bool(after_snapshot)


def _run_design_loop(args, struct_manager, runner) -> None:
    """Run the selected sequence-design backend across all pending structures."""
    attempted_structs = 0
    successful_structs = 0
    failed_structs = 0

    for pdb in struct_manager.iterate():
        attempted_structs += 1
        output_snapshot_before = _snapshot_output_state(struct_manager, pdb)

        if args.debug:
            runner.run_model(pdb, args)
            output_snapshot_after = _snapshot_output_state(struct_manager, pdb)
            if not _persisted_outputs_changed(output_snapshot_before, output_snapshot_after):
                raise RuntimeError(f"{args.backend} produced no persisted outputs for {pdb}")
            successful_structs += 1
            struct_manager.record_checkpoint(pdb)
            continue

        t0 = time.time()

        try:
            runner.run_model(pdb, args)
            output_snapshot_after = _snapshot_output_state(struct_manager, pdb)
            if not _persisted_outputs_changed(output_snapshot_before, output_snapshot_after):
                raise RuntimeError(f"{args.backend} produced no persisted outputs for {pdb}")

        except KeyboardInterrupt:
            sys.exit("Script killed by Control+C, exiting")

        except Exception:
            failed_structs += 1
            seconds = int(time.time() - t0)
            print(f"Struct with tag {pdb} failed in {seconds} seconds with error: {sys.exc_info()[0]}")
            continue

        successful_structs += 1
        struct_manager.record_checkpoint(pdb)

    if attempted_structs > 0 and successful_structs == 0:
        sys.exit(f"{args.backend} produced zero persisted outputs across {attempted_structs} attempted structures")

    if failed_structs:
        print(f"Completed with {failed_structs} failed structure(s) that were left uncheckpointed")

class ProteinMPNN_runner():
    '''
    This class is designed to run the ProteinMPNN model on a single input. This class handles the loading of the model,
    the loading of the input data, the running of the model, and the processing of the output
    '''

    def __init__(self, args, struct_manager):
        self.struct_manager = struct_manager

        if torch.cuda.is_available():
            print('Found GPU will run ProteinMPNN on GPU')
            self.device = "cuda:0"
        else:
            print('No GPU found, running ProteinMPNN on CPU')
            self.device = "cpu"

        self.mpnn_model = mpnn_util.init_seq_optimize_model(
            self.device,
            hidden_dim=128,
            num_layers = 3,
            backbone_noise = args.augment_eps,
            num_connections = args.num_connections,
            checkpoint_path = args.checkpoint_path
        )

        self.temperature = args.temperature
        self.seqs_per_struct = args.seqs_per_struct
        self.omit_AAs = [ letter for letter in args.omit_AAs.upper() if letter in list("ARNDCQEGHILKMFPSTWYVX") ]
        self.allow_x = args.allow_x

    def sequence_optimize(self, sample_feats: SampleFeatures) -> list[tuple[str, float]]:
        t0 = time.time()

        # ProteinMPNN backend is designed for antibody-target interface design;
        # it uses the last chain as visible context (target) and redesigns the rest.
        # Antibody-only inputs (no T chain) should use the AntiFold backend instead.
        assert 'T' in sample_feats.chains, (
            f"ProteinMPNN backend requires a target chain 'T'. Got chains {sample_feats.chains}. "
            "For antibody-only inputs, use '-backend antifold'."
        )

        # Once we have figured out pose I/O without Rosetta this will be easy to swap in
        pdbfile = 'temp.pdb'
        sample_feats.pose.dump_pdb(pdbfile)

        feature_dict = mpnn_util.generate_seqopt_features(pdbfile, sample_feats.chains)

        os.remove(pdbfile)

        arg_dict = mpnn_util.set_default_args(self.seqs_per_struct, omit_AAs=self.omit_AAs, allow_x=self.allow_x)
        arg_dict['temperature'] = self.temperature

        masked_chains = sample_feats.chains[:-1]
        visible_chains = [sample_feats.chains[-1]]

        fixed_positions_dict = {pdbfile[:-len('.pdb')]: sample_feats.fixed_res}

        sequences = mpnn_util.generate_sequences(
            self.mpnn_model,
            self.device,
            feature_dict,
            arg_dict,
            masked_chains,
            visible_chains,
            fixed_positions_dict=fixed_positions_dict
        )
        
        print( f"MPNN generated {len(sequences)} sequences in {int( time.time() - t0 )} seconds" ) 

        print(f'sequence_optimize: {sequences}')

        return sequences

    def proteinmpnn(self, sample_feats: SampleFeatures) -> None:
        '''
        Run MPNN sequence optimization on the pose
        '''
        seqs_scores = self.sequence_optimize(sample_feats)

        # Iterate though each seq score pair and thread the sequence onto the pose
        # Then write each pose to a pdb file
        prefix = f"{sample_feats.tag}_dldesign"
        for idx, (seq, _) in enumerate(seqs_scores): 
            sample_feats.thread_mpnn_seq(seq)

            outtag = f"{prefix}_{idx}"

            self.struct_manager.dump_pose(sample_feats.pose, outtag)

    def run_model(self, tag, args):
        t0 = time.time()

        print(f"Attempting pose: {tag}")
        
        # Load the pose 
        pose = self.struct_manager.load_pose(tag)

        # Initialize the features
        sample_feats = SampleFeatures(pose, tag)

        # Parse the loop string and determine which residues should be designed
        sample_feats.loop_string2fixed_res(args.loop_string)

        self.proteinmpnn(sample_feats)

        seconds = int(time.time() - t0)

        print(f"Struct: {tag} reported success in {seconds} seconds")


####################
####### Main #######
####################

def main(argv: list[str] | None = None) -> None:
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    _configure_determinism(args)

    struct_manager = StructManager(args)

    if args.backend == "antifold":
        from rfantibody.antifold.antifold_runner import AntiFold_runner

        runner = AntiFold_runner(args, struct_manager)
        print("Using AntiFold backend for sequence design")
    else:
        runner = ProteinMPNN_runner(args, struct_manager)
        print("Using ProteinMPNN backend for sequence design")

    _run_design_loop(args, struct_manager, runner)


if __name__ == "__main__":
    main()
    



