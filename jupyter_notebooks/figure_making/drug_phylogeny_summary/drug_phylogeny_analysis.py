#script to process making a drug phylogeny analysis by fingerprints
#written with the help of chatgpt to assist with logic and syntax for assembling the phylogeny

#package imports
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# Load your CSV
df = pd.read_csv('drug_summary_data.csv')

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

# Step 3: Calculate SSMD vs DMSO group
# Assumes "group" column contains a label like 'DMSO', 'Drug A', etc.
# We compute SSMD for each non-DMSO group vs DMSO for each group

def ssmd(x, y):
    return (np.mean(x) - np.mean(y)) / np.sqrt(np.var(x) + np.var(y))

# Collect SSMDs per group
dmso_vals = df[df['group'] == 'DMSO']['normalized_rlu_10uM']
ssmd_results = []

for grp in df['group'].unique():
    if grp == 'DMSO':
        continue
    grp_vals = df[df['group'] == grp]['normalized_rlu_10uM']
    ssmd_val = ssmd(grp_vals, dmso_vals)
    ssmd_results.append({'group': grp, 'SSMD_vs_DMSO': ssmd_val})

ssmd_df = pd.DataFrame(ssmd_results)

# Step 4: Merge SSMD results with main dataframe
df = df.merge(ssmd_df, on='group', how='left')

# Step 5: Plot SSMD Heatmap (1 value per group)
plt.figure(figsize=(10, 6))
ssmd_matrix = ssmd_df.pivot_table(index='group', values='SSMD_vs_DMSO')
sns.heatmap(ssmd_matrix, annot=True, cmap='vlag', center=0, linewidths=0.5)
plt.title('SSMD vs DMSO per Group')
plt.xlabel('SSMD')
plt.ylabel('Drug Group')
plt.tight_layout()
plt.savefig('my_plot.png', dpi=300, bbox_inches='tight')
#plt.show()