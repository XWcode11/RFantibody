# 目标PDB是抗原抗体复合物，抗原为链T，抗体为链H和L。使用GROMACS进行分子动力学模拟和FEP计算。以便
# 进行结合自由能计算。

# ========== 设置工作目录和文件路径 ==========
WORK_DIR="/mnt/data1/xiongw/Projects/active/RFantibody/output_data/gromacs"
INPUT_PDB="/mnt/data1/xiongw/Projects/active/RFantibody/output_data/rf2/hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1/hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1_0_dldesign_0_best.pdb"
BASE_NAME="hPDL1_N_hlt-1C2_Chothia_hlt-20250327-1_0_dldesign_0_best"
MDP_DIR="/mnt/data1/xiongw/Projects/active/RFantibody/output_data/gromacs/MDP_DIR"


# ========== 清理中间文件 (新增) ==========
# rm -f *.gro *.tpr *.top *.edr *.log *.cpt index.ndx *.xvg 

# ========== 控制每个步骤是否执行的开关 ==========
DO_PDB2GMX=true       # 生成拓扑文件
DO_MAKE_NDX=true      # 创建索引文件
DO_EDITCONF=true      # 定义模拟盒子
DO_SOLVATE=true       # 添加溶剂
DO_IONS=true          # 添加离子
DO_EM=true            # 能量最小化
DO_NVT=true           # NVT平衡
DO_NPT=true          # NPT平衡
DO_MD=true           # 生产性MD模拟
DO_FEP_PREP=true     # 准备FEP拓扑
DO_FEP_WINDOWS=true  # 创建FEP窗口
DO_FEP_ANALYSIS=false # 分析FEP结果

# ========== 日志设置 ==========
LOG_FILE="${WORK_DIR}/gromacs_fep_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "开始运行时间: $(date)"

# ========== 基础准备 ==========
cd "$WORK_DIR" || { echo "错误: 无法进入工作目录"; exit 1; }

# 检查 GROMACS 是否可用
if ! command -v gmx &> /dev/null; then
    echo "错误: GROMACS 未安装或未正确配置环境变量"
    exit 1
fi

# ========== 步骤1: 生成拓扑文件和初始结构 ==========
if $DO_PDB2GMX; then
    echo "步骤 1: 生成拓扑文件和初始结构..."
    gmx pdb2gmx -f "$INPUT_PDB" \
        -o "${BASE_NAME}.gro" \
        -water spce \
        -ff amber99sb \
        -ignh || { echo "pdb2gmx 失败"; exit 1; }
    [ -f "topol.top" ] || { echo "topol.top 未生成"; exit 1; }
    echo "✓ 步骤1已完成"
else
    echo "步骤 1: 跳过拓扑文件生成"
fi


# ========== 步骤3: 定义模拟盒子 ==========
if $DO_EDITCONF; then
    echo "步骤 3: 定义模拟盒子..."
    gmx editconf -f "${BASE_NAME}.gro" \
        -o "${BASE_NAME}-box.gro" \
        -d 1.2 \
        -bt cubic || { echo "editconf 失败"; exit 1; }
    echo "✓ 步骤3已完成"
else
    echo "步骤 3: 跳过模拟盒子定义"
fi

# ========== 步骤4: 添加溶剂 ==========
if $DO_SOLVATE; then
    echo "步骤 4: 添加溶剂..."
    # 备份原始拓扑文件以确保正确恢复
    cp topol.top topol.top.bak

    gmx solvate -cp "${BASE_NAME}-box.gro" \
        -cs spc216.gro \
        -o "${BASE_NAME}_solv.gro" \
        -p topol.top || { echo "solvate 失败"; exit 1; }

    # 验证溶剂化后的拓扑与结构匹配
    SOLV_ATOMS=$(grep -c "ATOM" "${BASE_NAME}_solv.gro")
    echo "溶剂化后结构中的原子数: $SOLV_ATOMS"

    echo "✓ 步骤4已完成"
else
    echo "步骤 4: 跳过溶剂添加"
fi

# ========== 步骤5: 添加离子 ==========
if $DO_IONS; then
    echo "步骤 5: 添加离子..."

    # 先检查结构文件和拓扑文件是否匹配
    if [ -f "${BASE_NAME}_solv.gro" ] && [ -f "topol.top" ]; then
        # 获取拓扑文件中的原子数（估计值）
        TOP_ATOMS=$(grep -A 1000 "\[ molecules \]" topol.top | grep -v "\[" | awk '{sum+=$2} END {print sum}')

        # 获取结构文件中的原子数
        GRO_ATOMS=$(grep -c "ATOM" "${BASE_NAME}_solv.gro" || wc -l < "${BASE_NAME}_solv.gro")

        echo "拓扑文件中的预估原子数: $TOP_ATOMS"
        echo "结构文件中的原子数: $GRO_ATOMS"

        # 如果不匹配，尝试修复
        if [ "$TOP_ATOMS" != "$GRO_ATOMS" ]; then
            echo "警告: 拓扑文件与结构文件不匹配，尝试修复..."

            # 恢复备份的拓扑文件
            if [ -f "topol.top.bak" ]; then
                cp topol.top.bak topol.top
                echo "已恢复原始拓扑文件"
            fi

            # 重新生成溶剂化系统
            echo "重新进行溶剂化..."
            gmx solvate -cp "${BASE_NAME}-box.gro" \
                -cs spc216.gro \
                -o "${BASE_NAME}_solv.gro" \
                -p topol.top || { echo "溶剂化失败"; exit 1; }
        fi
    else
        echo "错误: 缺少溶剂化结构文件或拓扑文件"
        exit 1
    fi

    # 使用改进的离子添加步骤，增加maxwarn
    gmx grompp -f "${MDP_DIR}/ions.mdp" \
        -c "${BASE_NAME}_solv.gro" \
        -p topol.top \
        -o "${BASE_NAME}-ions.tpr" \
        -maxwarn 2 || { echo "grompp for ions 失败"; exit 1; }

    echo -e "SOL" | gmx genion -s "${BASE_NAME}-ions.tpr" \
        -o "${BASE_NAME}_ions.gro" \
        -p topol.top \
        -pname NA \
        -nname CL \
        -neutral || { echo "genion 失败"; exit 1; }

    echo "✓ 步骤5已完成"
else
    echo "步骤 5: 跳过离子添加"
fi

# ========== 步骤6: 能量最小化 ==========
if $DO_EM; then
    echo "步骤 6: 能量最小化..."
    # 检查文件是否存在
    if [ ! -f "${BASE_NAME}_ions.gro" ]; then
        echo "错误: 离子化结构文件不存在"
        exit 1
    fi

    gmx grompp -f "${MDP_DIR}/minim.mdp" \
        -c "${BASE_NAME}_ions.gro" \
        -p topol.top \
        -o em.tpr \
        -maxwarn 3 || { echo "grompp for minimization 失败"; exit 1; }

    gmx mdrun -v -deffnm em || { echo "mdrun for minimization 失败"; exit 1; }

    # 验证能量最小化结果
    if [ -f "em.gro" ] && [ -f "em.edr" ]; then
        echo "能量最小化完成，生成了em.gro和em.edr文件"
    else
        echo "警告: 能量最小化可能未正确完成"
    fi

    echo "✓ 步骤6已完成"
else
    echo "步骤 6: 跳过能量最小化"
fi

# ========== 步骤7: NVT平衡 ==========
if $DO_NVT; then
    echo "步骤 7: NVT平衡..."
    if [ ! -f "em.gro" ]; then
        echo "错误: 能量最小化结构文件不存在"
        exit 1;
    fi

    if $DO_MAKE_NDX; then
        echo " 生成索引文件..."
        # 直接使用系统自动生成的默认组，不做特殊分组
        # 0号组"System"默认包含所有原子，完全满足nvt.mdp和npt.mdp的需求
        echo -e "q" | gmx make_ndx -f em.gro -o index.ndx || { echo "make_ndx 失败"; exit 1; }
        
        # 显示可用的组
        echo "索引组已生成，可用组:"
        gmx make_ndx -f em.gro -n index.ndx -o index.ndx < /dev/null 2>&1 | grep -A 20 "Group"
        echo "✓ 步骤2已完成"
    else
        echo "跳过索引文件生成"
    fi
    
    # 创建位置约束文件
    echo "为NVT平衡创建位置约束文件..."
    echo -e "Protein\nSystem" | gmx genrestr -f em.gro -o posre.itp -fc 1000 1000 1000 || { echo "genrestr 失败"; exit 1; }
    
    # 检查或创建#ifdef POSRES语句
    for itp in topol_Protein_*.itp; do
        if ! grep -q "#ifdef POSRES" "$itp"; then
            # 将位置约束添加到拓扑文件中
            posre_section="\n; Position restraint for NVT equilibration\n#ifdef POSRES\n#include \"posre.itp\"\n#endif\n"
            sed -i "/\[ moleculetype \]/i ${posre_section}" "$itp"
            echo "已添加位置约束到 $itp"
        fi
    done

    echo "运行一个非常温和的初始NVT平衡..."
    gmx grompp -f "${MDP_DIR}/nvt.mdp" \
        -c em.gro \
        -p topol.top \
        -n index.ndx \
        -r em.gro \
        -o nvt.tpr \
        -maxwarn 3 || { echo "grompp for NVT 失败"; exit 1; }

    # 使用更稳定的CPU模式和额外的稳定性选项 - 移除noconfout选项以保存输出
    gmx mdrun -v -deffnm nvt \
        -nb cpu \
        -ntmpi 1 \
        -ntomp 4 \
        -rdd 1.4 \
        -dlb yes || { echo "mdrun for NVT 失败"; exit 1; }
    
    # 确认输出文件存在
    if [ ! -f "nvt.gro" ]; then
        echo "错误: NVT平衡没有生成nvt.gro文件"
        exit 1
    fi
        
    echo "✓ 步骤7已完成"
else
    echo "步骤 7: 跳过NVT平衡"
fi

# ========== 步骤8: NPT平衡 ==========
if $DO_NPT; then
    echo "步骤 8: NPT平衡..."
    
    # 检查NVT输出文件是否存在
    if [ ! -f "nvt.gro" ]; then
        echo "错误: 找不到nvt.gro文件，无法进行NPT平衡"
        exit 1
    fi
    
    # 为NPT平衡更新位置约束强度（如果存在posre.itp文件）
    if [ -f "posre.itp" ]; then
        echo "为NPT减弱位置约束强度..."
        cp posre.itp posre.itp.nvt
        # 将力常数降低到原来的一半
        awk '{if (NF>=7 && $1 !~ /^\;/) {print $1,$2,$3,$4/2,$5/2,$6/2} else {print $0}}' posre.itp.nvt > posre.itp
    fi

    # 使用改进的NPT设置
    gmx grompp -f "${MDP_DIR}/npt.mdp" \
        -c nvt.gro \
        -p topol.top \
        -n index.ndx \
        -r nvt.gro \
        -o npt.tpr \
        -maxwarn 2 || { echo "grompp for NPT 失败"; exit 1; }
    
    # 使用与NVT相同的稳定性设置
    gmx mdrun -v -deffnm npt \
        -nb cpu \
        -ntmpi 1 \
        -ntomp 4 \
        -rdd 1.4 \
        -dlb yes || { echo "mdrun for NPT 失败"; exit 1; }
    
    echo "✓ 步骤8已完成"
else
    echo "步骤 8: 跳过NPT平衡"
fi

# ========== 步骤9: 生产性MD模拟 ==========
if $DO_MD; then
    echo "步骤 9: 生产性MD模拟..."
    
    # 检查NPT输出文件是否存在
    if [ ! -f "npt.gro" ]; then
        echo "错误: 找不到npt.gro文件，无法进行MD模拟"
        exit 1
    fi
    
    
    # 显示索引文件中的组以便诊断
    echo "索引文件中的组:"
    gmx make_ndx -f npt.gro -n index.ndx -o index.ndx < /dev/null 2>&1 | grep -A 20 "Group"

    # 移除位置约束（如果存在于拓扑文件中）
    for itp in topol_Protein_*.itp; do
        if grep -q "#ifdef POSRES" "$itp"; then
            echo "从$itp中移除位置约束标记..."
            sed -i '/#ifdef POSRES/,/#endif/d' "$itp"
        fi
    done
    
    # 使用改进的MD设置
    gmx grompp -f "${MDP_DIR}/md.mdp" \
        -c npt.gro \
        -t npt.cpt \
        -p topol.top \
        -n index.ndx \
        -o md.tpr \
        -maxwarn 3 || { echo "grompp for MD 失败"; exit 1; }
    
    # 使用与NPT相同的稳定性设置
    gmx mdrun -v -deffnm md \
        -nb cpu \
        -ntmpi 1 \
        -ntomp 4 \
        -dlb yes || { echo "mdrun for MD 失败"; exit 1; }
    
    echo "✓ 步骤9已完成"
else
    echo "步骤 9: 跳过生产性MD模拟"
fi
# 计算抗原-抗体结合能