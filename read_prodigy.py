import pandas as pd
import glob
import os

file_dir = "output_data/prodigy"
file_list = glob.glob(os.path.join(file_dir, "*.txt"))

# 初始化一个字典来存储数据
data = {"pdb": [],  "Predicted binding affinity (kcal.mol-1)": [], "Predicted dissociation constant (M) at 25.0˚C": []}

for file in file_list:
    filename = os.path.basename(file)
    with open(file, 'r', encoding='utf-8') as f:
        # 读取最后两行并去掉换行符
        last_two_lines = [line.strip() for line in f.readlines()[-2:]]
        # 提取列名和值
        column1, value1 = last_two_lines[0].split(":")[0][5:], last_two_lines[0].split(":")[1].strip()
        column2, value2 = last_two_lines[1].split(":")[0][5:], last_two_lines[1].split(":")[1].strip()
        # 添加到字典
        data["pdb"].append(filename)
        data["Predicted binding affinity (kcal.mol-1)"].append(value1)
        data["Predicted dissociation constant (M) at 25.0˚C"].append(value2)

# 转换为 DataFrame
df = pd.DataFrame(data)
df.to_csv("output_data/prodigy/prodigy.csv", index=False)