import numpy as np
from scipy import stats



dmso = np.array([262.5238,109.9048,108.0952,116.1905,102.2024,342.1429,

                 87.6548,130.4167,307.0357,219.6905,107.3929,192.4762])



drug = np.array([128.4881,41.4167,36.8690,118.5119,92.2262,76.9167,

                 34.7619,113.0476,67.3333,137.8214,180.0,164.6429])



print("Mann–Whitney U:", stats.mannwhitneyu(dmso, drug, alternative='two-sided'))
print("Welch t-test:", stats.ttest_ind(dmso, drug, equal_var=False))



print("\nNormality tests (Shapiro-Wilk):")
print("DMSO:", stats.shapiro(dmso))
print("Drug:", stats.shapiro(drug))
