
antigen_hlt_pdb=$1
antibody_hlt_pdb=$2
tag=$3
# If tag is not provided, use today's date as the tag
if [ -z "$tag" ]; then
    tag=$(date +%Y%m%d%H%M)
fi
# get filenames and remove the extension
antigen_filename=$(basename "$antigen_hlt_pdb" .pdb)

antibody_filename=$(basename "$antibody_hlt_pdb" .pdb)

output_name=${antigen_filename}-${antibody_filename}


# create a directory to store the rfdiffusion results
rfdiffusion_output_dir=/home/output_data/rfdiffusion/$output_name-${tag}
mkdir -p $rfdiffusion_output_dir

# create a directory to store the proteinmpnn results
proteinmpnn_output_dir=/home/output_data/proteinmpnn/$output_name-${tag}
mkdir -p $proteinmpnn_output_dir

# create a directory to store the rf2 results
rf2_output_dir=/home/output_data/rf2/$output_name-${tag}
mkdir -p $rf2_output_dir

output_prefix=$rfdiffusion_output_dir/${output_name}-${tag}

poetry run python  /home/src/rfantibody/rfdiffusion/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=$antigen_hlt_pdb  \
    antibody.framework_pdb=$antibody_hlt_pdb \
    inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \
    'ppi.hotspot_res=[T103,T104,T105,T106,T107,T1,T48,T95,T8]' \
    'antibody.design_loops=[L1:,L2:,L3:,H1:,H2:,H3:]' \
    inference.num_designs=20 \
    inference.output_prefix=$output_prefix



poetry run python /home/scripts/proteinmpnn_interface_design.py \
    -pdbdir $rfdiffusion_output_dir \
    -outpdbdir $proteinmpnn_output_dir



poetry run python /home/scripts/rf2_predict.py \
    input.pdb_dir=$proteinmpnn_output_dir \
    output.pdb_dir=$rf2_output_dir 

# extract the scores from the rf2 results
poetry run python extract_pdb_scores.py $rf2_output_dir --output $rf2_output_dir/${output_name}-${tag}-scores.csv 

echo $output_name "finished!"