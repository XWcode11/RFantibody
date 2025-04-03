#!/bin/bash
# 运行 rosetta 的脚本
ROSETTA_PATH=/home/xw/rosetta/rosetta.binary.linux.release-371/main/source/bin
subdir=$1
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
        -interface LH_T -pack_separated  -out:file:score_only $rosetta_outdir/${subdir}_relax.out
done

echo "转换Rosetta结果为标准CSV格式..."
for relax_file in ${rosetta_outdir}/*_relax.out; do
    if [ -f "$relax_file" ]; then
        # 使用相同文件名，仅更改扩展名
        output_file="${relax_file/.out/.csv}"
        
        # 处理数据：跳过SEQUENCE行，只处理SCORE开头的行，转换空格为逗号
        grep "^SCORE:" "$relax_file" | tr -s ' ' | tr ' ' ',' > "$output_file"
        
        echo "已转换: $(basename $relax_file) → $(basename $output_file)"
    fi
done

echo "处理完成。结果保存在: $rosetta_outdir"
