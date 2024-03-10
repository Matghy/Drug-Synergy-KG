import pandas as pd
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate

# 初始化sigFactory，调整bins定义
fdefName = 'MinimalFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)
sigFactory = SigFactory(featFactory, minPointCount=2, maxPointCount=3,trianglePruneBins=False)
sigFactory.SetBins([(0,2), (2,5), (5,8)])  # 调整了这里的距离区间
sigFactory.Init()

# 读取CSV文件中的SMILES列
df = pd.read_csv('drug_init_smiles.csv')  # 确保CSV文件位于脚本的同一目录下，或者提供完整路径
smiles_list = df['SMILES'].tolist()

# 初始化存储有效fingerprint的列表
valid_fps = []

# 为每个SMILES字符串生成fingerprint
for i, smiles in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(smiles)
    if mol:  # 确保分子对象有效
        try:
            fp = Generate.Gen2DFingerprint(mol, sigFactory)
            # 将fingerprint转换为数字串，并确保开头的零被保留
            fp_bits = [int(bit) for bit in list(fp.ToBitString())]
            valid_fps.append(fp_bits)
        except Exception as e:
            # 遇到有问题的分子时跳过，并可选择打印错误信息
            print(f"Skipping molecule at index {i} due to error: {e}")
            valid_fps.append([None] * sigFactory.GetSigSize())  # 使用None或适当的值填充，保持行索引对齐
    else:
        print(f"Invalid SMILES at index {i}: {smiles}")
        valid_fps.append([None] * sigFactory.GetSigSize())  # 无效的SMILES，使用None填充

# 将数据转换为DataFrame
fp_df = pd.DataFrame(valid_fps)

# 将原始SMILES数据与fingerprints合并
result_df = pd.concat([df, fp_df], axis=1)

# 保存到新的CSV文件
result_df.to_csv('drug_init_embs.csv', index=False)
