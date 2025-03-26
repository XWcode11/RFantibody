import argparse
import glob
import os
import pandas as pd

def extract_sequence(lines):
    """从PDB文件中提取重链(H)和轻链(L)的氨基酸序列"""
    # 氨基酸三字母代码到单字母代码的映射
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    chains = {}  # 存储所有链的残基信息
    

    for line in lines:
        if line.startswith("ATOM"):
            # 只考虑CA原子以避免重复残基
            if line[12:16].strip() == "CA":
                chain_id = line[21:22].strip()
                residue_name = line[17:20].strip()
                residue_number = int(line[22:26].strip())
                
                if residue_name in aa_map:
                    amino_acid = aa_map[residue_name]
                    
                    if chain_id not in chains:
                        chains[chain_id] = {}
                    
                    # 只记录新的残基，避免重复
                    if residue_number not in chains[chain_id]:
                        chains[chain_id][residue_number] = amino_acid
    
    # 提取H和L链序列
    h_sequence = ""
    l_sequence = ""
    
    # 如果明确标记了H和L链
    if "H" in chains:
        h_chain_residues = chains["H"]
        h_sequence = ''.join([h_chain_residues[num] for num in sorted(h_chain_residues.keys())])
    
    if "L" in chains:
        l_chain_residues = chains["L"]
        l_sequence = ''.join([l_chain_residues[num] for num in sorted(l_chain_residues.keys())])
    
    return h_sequence, l_sequence

def extract_metrics(file):
    with open(file) as f:
        lines = f.readlines()
        
    # 提取序列
    h_sequence, l_sequence = extract_sequence(lines)
    
    # 提取metrics
    metrics_lines = lines[-13:]
    records = []
    for line in metrics_lines:
        parts = line.strip().split(':')
        if len(parts) < 2:
            continue
        key = parts[0].strip()
        try:
            value = float(parts[1].strip())
        except ValueError:
            continue
        records.append({
            "file_name": os.path.basename(file), 
            "metrics": key, 
            "value": value,
            "heavy_chain": h_sequence,
            "light_chain": l_sequence
        })
    
    return records

def main():
    parser = argparse.ArgumentParser(description="从PDB文件中提取metrics并保存到CSV")
    parser.add_argument("folder", help="包含PDB文件的文件夹路径")
    parser.add_argument("--output", default="output.csv", help="输出CSV文件路径，默认为output.csv")
    args = parser.parse_args()
    
    pdb_files = glob.glob(os.path.join(args.folder, "*.pdb"))
    all_records = []
    for file in pdb_files:
        all_records.extend(extract_metrics(file))
    
    if all_records:
        df = pd.DataFrame(all_records)
        df.to_csv(args.output, index=False)
        print(f"CSV文件已保存到 {args.output}")
    else:
        print("未找到符合条件的metrics记录。")

def test():
    file = "output_data/rf2/hPDL1_N_hlt_1C2_Chothia_hlt/hPDL1_N_hlt_1C2_Chothia_hlt_0_dldesign_0_best.pdb"
    records = extract_metrics(file)
    print(records)
if __name__ == "__main__":
    main()