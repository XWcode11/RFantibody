
# Iterate through all pdb files in the rf2 output directory and run prodigy on each
for pdb_file in output_data/rf2/hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1/*.pdb; do
    # 获取文件名
    filename=$(basename "$pdb_file")
    # 获取文件名的前缀（去掉扩展名）
    prefix="${filename%.*}"
    prodigy "$pdb_file" --selection H,L T > output_data/prodigy/${prefix}.txt
done




