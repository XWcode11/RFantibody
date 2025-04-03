#!/bin/bash
# 一站式抗体设计与评估脚本

# 检查命令行参数
if [ "$#" -lt 2 ]; then
    echo "用法: $0 <抗原PDB文件> <抗体PDB文件> [可选:标签]"
    exit 1
fi

# 获取参数
ANTIGEN_PDB=$1  #相对路径
ANTIBODY_PDB=$2 #相对路径
TAG=$3

# 如果没有提供TAG，使用当前日期
if [ -z "$TAG" ]; then
    TAG=$(date +%Y%m%d-%H%M)
fi

# 获取文件名（不带路径和扩展名）
ANTIGEN_NAME=$(basename "$ANTIGEN_PDB" .pdb)
ANTIBODY_NAME=$(basename "$ANTIBODY_PDB" .pdb)
OUTPUT_NAME="${ANTIGEN_NAME}-${ANTIBODY_NAME}-${TAG}"

echo "========================================================"
echo "开始抗体设计流程"
echo "抗原: $ANTIGEN_PDB"
echo "抗体: $ANTIBODY_PDB"
echo "标签: $TAG"
echo "输出名称: $OUTPUT_NAME"
echo "========================================================"

# 1. 启动Docker容器
echo "步骤1: 启动Docker容器..."
docker start rfantibody_xw

# 2. 在容器内运行pipeline脚本
echo "步骤2: 运行抗体设计pipeline..."
docker exec rfantibody_xw bash -c "cd /home && ./run_pipeline.sh $ANTIGEN_PDB $ANTIBODY_PDB $TAG"

# 3. 准备Rosetta输入
echo "步骤3: 准备Rosetta结构优化..."
# 更新run_rosetta.sh中的子目录路径
subdir="${ANTIGEN_NAME}-${ANTIBODY_NAME}-${TAG}"

# 4. 运行Rosetta分析
echo "步骤4: 运行Rosetta结构优化和分析..."
./run_rosetta.sh $subdir

echo "========================================================"
echo "抗体设计和评估流程完成!"
echo "RFdiffusion结果: output_data/rfdiffusion/$subdir"
echo "ProteinMPNN结果: output_data/proteinmpnn/$subdir"
echo "RF2结果: output_data/rf2/$subdir"
echo "Rosetta结果: output_data/rosetta/$subdir" 
echo "========================================================"