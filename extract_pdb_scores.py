import argparse
import glob
import os
import pandas as pd

def extract_metrics(file):
    with open(file) as f:
        lines = f.readlines()
    # 假设 metrics 在文件末尾 13 行内
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
        records.append({"file_name": os.path.basename(file), "metrics": key, "value": value})
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

if __name__ == "__main__":
    main()