#!/bin/bash

# 脚本名称: run_gromacs_fep.sh
# 目的: 使用 GROMACS 和 FEP 方法计算抗原-抗体结合能

# 设置工作目录和文件路径
WORK_DIR="/mnt/data1/xiongw/Projects/active/RFantibody/output_data/gromacs"
INPUT_PDB="/mnt/data1/xiongw/Projects/active/RFantibody/output_data/rf2/hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1/hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1_0_dldesign_0_best.pdb"
BASE_NAME="hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1_0_dldesign_0_best"
MDP_DIR="/mnt/data1/xiongw/Projects/active/RFantibody/output_data/gromacs/MDP_DIR"

# 创建工作目录（如果不存在）
# mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit

# 检查 GROMACS 是否可用
if ! command -v gmx &> /dev/null; then
    echo "错误: GROMACS 未安装或未正确配置环境变量"
    exit 1
fi

# 1. 生成拓扑文件和初始结构
echo "步骤 1: 生成拓扑文件和初始结构..."
gmx pdb2gmx -f "$INPUT_PDB" \
    -o "${BASE_NAME}.gro" \
    -water spce \
    -ff amber99sb \
    -ignh || { echo "pdb2gmx 失败"; exit 1; }
[ -f "topol.top" ] || { echo "topol.top 未生成"; exit 1; }

# 创建索引文件，定义抗原和抗体组（H 和 L 为抗体，T 为抗原）
echo "生成索引文件..."
echo -e "chain H L\nchain T\nq" | gmx make_ndx -f "${BASE_NAME}.gro" -o index.ndx || { echo "make_ndx 失败"; exit 1; }

# 2. 定义模拟盒子
echo "步骤 2: 定义模拟盒子..."
gmx editconf -f "${BASE_NAME}.gro" \
    -o "${BASE_NAME}-box.gro" \
    -d 1.2 \
    -bt cubic || { echo "editconf 失败"; exit 1; }


# 3. 添加溶剂（水分子）
echo "步骤 3: 添加溶剂..."
gmx solvate -cp "${BASE_NAME}-box.gro" \
    -cs spc216.gro \
    -o "${BASE_NAME}_solv.gro" \
    -p topol.top || { echo "solvate 失败"; exit 1; }

# 4. 添加离子以中和电荷
echo "步骤 4: 添加离子..."
gmx grompp -f "${MDP_DIR}/ions.mdp" \
    -c "${BASE_NAME}_solv.gro" \
    -p topol.top \
    -maxwarn 1 \
    -o "${BASE_NAME}-ions.tpr" || { echo "grompp for ions 失败"; exit 1; }
gmx genion -s "${BASE_NAME}-ions.tpr" \
    -o "${BASE_NAME}_ions.gro" \
    -p topol.top \
    -pname NA \
    -nname CL \
    -neutral << EOF
SOL
EOF
[ $? -eq 0 ] || { echo "genion 失败"; exit 1; }

# 5. 能量最小化
echo "步骤 5: 能量最小化..."
gmx grompp -f "${MDP_DIR}/minim.mdp" \
    -c "${BASE_NAME}_ions.gro" \
    -p topol.top \
    -o em.tpr || { echo "grompp for minimization 失败"; exit 1; }
gmx mdrun -v -deffnm em || { echo "mdrun for minimization 失败"; exit 1; }

# 6. NVT 平衡
echo "步骤 6: NVT 平衡..."
gmx grompp -f "${MDP_DIR}/nvt.mdp" \
    -c em.gro \
    -p topol.top \
    -n index.ndx \
    -o nvt.tpr || { echo "grompp for NVT 失败"; exit 1; }
gmx mdrun -v -deffnm nvt || { echo "mdrun for NVT 失败"; exit 1; }

# 7. NPT 平衡
echo "步骤 7: NPT 平衡..."
gmx grompp -f "${MDP_DIR}/npt.mdp" \
    -c nvt.gro \
    -p topol.top \
    -n index.ndx \
    -o npt.tpr || { echo "grompp for NPT 失败"; exit 1; }
gmx mdrun -v -deffnm npt || { echo "mdrun for NPT 失败"; exit 1; }

# 8. 生产性 MD 模拟
echo "步骤 8: 生产性 MD 模拟..."
gmx grompp -f "${MDP_DIR}/md.mdp" \
    -c npt.gro \
    -p topol.top \
    -n index.ndx \
    -o md.tpr || { echo "grompp for MD 失败"; exit 1; }
gmx mdrun -v -deffnm md || { echo "mdrun for MD 失败"; exit 1; }

# 生成 FEP 拓扑文件
echo "生成 FEP 拓扑文件..."
cp topol_Protein_chain_T.itp topol_Protein_chain_T_fep.itp
sed -i '/\[ atoms \]/,/^\s*$/ s/\(\s*[0-9]\+\s\+[A-Za-z0-9]\+\s\+[0-9]\+\s\+[A-Za-z]\+\s\+[A-Za-z0-9]\+\s\+[0-9]\+\s\+[-0-9.]\+\s\+[0-9.]\+\)/\1 \2 0.0 \3/' topol_Protein_chain_T_fep.itp
sed -i '/#include "topol_Protein_chain_T.itp"/a #include "topol_Protein_chain_T_fep.itp"' topol.top

# 9. FEP 准备：创建多个 λ 窗口的目录和拓扑
echo "步骤 9: 准备 FEP 模拟..."
mkdir -p fep_windows
cd fep_windows || exit

# 定义 λ 值（21 个窗口）
LAMBDA_VALUES=(0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00)
for i in "${!LAMBDA_VALUES[@]}"; do
    L_DIR="lambda_${LAMBDA_VALUES[$i]}"
    mkdir -p "$L_DIR"
    cd "$L_DIR" || exit
    
    # 复制初始结构和拓扑
    cp "../../md.gro" "./${BASE_NAME}_start.gro"
    cp "../../topol.top" "./topol.top"
    cp "../../topol_Protein_chain_H.itp" "./"
    cp "../../topol_Protein_chain_L.itp" "./"
    cp "../../topol_Protein_chain_T.itp" "./"
    cp "../../topol_Protein_chain_T_fep.itp" "./"
    cp "../../index.ndx" "./index.ndx"
    cp "../../${MDP_DIR}/fep.mdp" "./fep.mdp"
    
    # 生成 .tpr 文件，指定当前 λ 状态
    gmx grompp -f fep.mdp \
        -c "${BASE_NAME}_start.gro" \
        -p topol.top \
        -n index.ndx \
        -o fep.tpr \
        -maxwarn 1 \
        -l "$i" || { echo "grompp for λ=${LAMBDA_VALUES[$i]} 失败"; exit 1; }
    
    # 运行 FEP 模拟（后台运行）
    gmx mdrun -v -deffnm fep &
    cd .. || exit
done
wait  # 等待所有后台任务完成

# 10. 分析 FEP 结果
echo "步骤 10: 分析 FEP 结果..."
cd "$WORK_DIR" || exit
gmx bar -f fep_windows/lambda_*/fep.xvg -o bar.xvg -b 1000 || { echo "gmx bar 失败"; exit 1; }

echo "完成！请查看 bar.xvg 文件以获取结合自由能结果。"