# Histogram of each group 1 ligand's highest Tanimoto similarity to any other group 1 ligand.
# Fingerprints: Morgan, radius 2, 1024 bits (same as drug_phylogeny_analysis.py).

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

INPUT_CSV = "drug_summary_data.csv"
OUT_PNG = "group1_max_tanimoto_histogram.png"
OUT_CSV = "group1_max_tanimoto_per_ligand.csv"
GROUP = 1
EXCLUDE = ["Z1312379811", "PV-000555062637", "PV-003960018814"]
EXPECTED_N = 28

df = pd.read_csv(INPUT_CSV)
df = df[(df["group"] == GROUP) & (~df["enamine_name"].isin(EXCLUDE))].reset_index(drop=True)
assert len(df) == EXPECTED_N, f"Expected {EXPECTED_N} ligands, got {len(df)}"

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
fps = []
for name, smi in zip(df["enamine_name"], df["smiles"]):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {name}: {smi}")
    fps.append(gen.GetFingerprint(mol))

max_sim, best_match = [], []
for i, fp in enumerate(fps):
    sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
    sims[i] = -1  # exclude self
    j = int(sims.argmax())
    max_sim.append(sims[j])
    best_match.append(df["enamine_name"][j])

df["max_tanimoto"] = max_sim
df["most_similar_ligand"] = best_match
df[["enamine_name", "smiles", "max_tanimoto", "most_similar_ligand"]].to_csv(OUT_CSV, index=False)
print(df[["enamine_name", "max_tanimoto", "most_similar_ligand"]].to_string())

bins = np.linspace(0, 1, 11)
counts, _ = np.histogram(max_sim, bins=bins)
print("Bin counts:", dict(zip([f"{b:.1f}-{b+0.1:.1f}" for b in bins[:-1]], counts)))

fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(max_sim, bins=bins, edgecolor="black", color="steelblue")
ax.set_xticks(bins)
ax.set_xlabel("Highest Tanimoto similarity to another ligand")
ax.set_ylabel("Number of ligands")
ax.set_title(f"Group {GROUP} nearest-neighbor Tanimoto similarity (n={len(df)})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
