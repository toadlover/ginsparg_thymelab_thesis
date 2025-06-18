import pandas as pd
import numpy as np
#import seaborn as sns
#import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

print("hi")

df = pd.read_csv('drug_summary_data.csv')

print("initial df")
print(df)

# Step 1: Preprocess
# Make sure the SMILES strings are valid
df['mol'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
df = df[df['mol'].notnull()].copy()  # drop invalid molecules

# Step 2: Calculate Morgan Fingerprints
def get_morgan_fp(mol, radius=2, nBits=1024):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    arr = np.zeros((nBits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

df['fingerprint'] = df['mol'].apply(get_morgan_fp)

print(df)

df.to_csv("with_fingerprints.csv", index=False)
