#!/bin/bash
# 运行 rosetta 的脚本
ROSETTA_PATH=/home/xw/rosetta/rosetta.binary.linux.release-371/main/source/bin
subdir=hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1
PDB_DIR=./output_data/rf2/${subdir}
rosetta_outdir=./output_data/rosetta/${subdir}

#create output directory
mkdir -p $rosetta_outdir

# 迭代所有 pdb 文件并运行 rosetta
for pdb_file in ${PDB_DIR}/*.pdb; do
    # 获取文件名
    filename=$(basename "$pdb_file")
    # 获取文件名的前缀（去掉扩展名）
    prefix="${filename%.*}"
    # 运行 rosetta

    echo $prefix
    # ${ROSETTA_PATH}/idealize_jd2.static.linuxgccrelease -in:file:s $pdb_file \
    #     -out:path:pdb $rosetta_outdir -score:weights ref2015 -out:pdb  -overwrite 
    ${ROSETTA_PATH}/relax.static.linuxgccrelease -in:file:s $pdb_file \
        -out:pdb -out:path:pdb $rosetta_outdir -relax:thorough -overwrite
    ${ROSETTA_PATH}/InterfaceAnalyzer.static.linuxgccrelease -s $rosetta_outdir/${prefix}_0001.pdb \
        -interface LH_T -pack_separated  -out:file:score_only $rosetta_outdir/${subdir}_relax.csv
done
