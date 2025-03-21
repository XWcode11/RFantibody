
antigen_hlt_pdb=$1
antibody_hlt_pdb=$2

# get filenames and remove the extension
antigen_filename=$(basename "$antigen_hlt_pdb")
antigen_filename="${antigen_filename%.*}"

antibody_filename=$(basename "$antibody_hlt_pdb")
antibody_filename="${antibody_filename%.*}"

output_name=${antigen_filename}_${antibody_filename}


# create a directory to store the rfdiffusion results
rfdiffusion_output_dir=/home/output_data/rfdiffusion/$output_name
mkdir -p $rfdiffusion_output_dir

# create a directory to store the proteinmpnn results
proteinmpnn_output_dir=/home/output_data/proteinmpnn/$output_name
mkdir -p $proteinmpnn_output_dir

# create a directory to store the rf2 results
rf2_output_dir=/home/output_data/rf2/$output_name
mkdir -p $rf2_output_dir

output_prefix=$rfdiffusion_output_dir/$output_name

poetry run python  /home/src/rfantibody/rfdiffusion/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=$antigen_hlt_pdb  \
    antibody.framework_pdb=$antibody_hlt_pdb \
    inference.ckpt_override_path=/home/weights/RFdiffusion_Ab.pt \
    'ppi.hotspot_res=[T305,T456]' \
    'antibody.design_loops=[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]' \
    inference.num_designs=20 \
    inference.output_prefix=$output_prefix



poetry run python /home/scripts/proteinmpnn_interface_design.py \
    -pdbdir $rfdiffusion_output_dir \
    -outpdbdir $proteinmpnn_output_dir



poetry run python /home/scripts/rf2_predict.py \
    input.pdb_dir=$proteinmpnn_output_dir \
    output.pdb_dir=$rf2_output_dir 

# extract the scores from the rf2 results
poetry run python extract_pdb_scores.py $rf2_output_dir --output $rf2_output_dir/${output_name}_scores.csv 

echo $output_name "finished!"