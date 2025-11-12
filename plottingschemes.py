import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

#trail-and-error manual fine-tuning of z values to find outliers in the trend between solubility and molecular weight

df = pd.read_csv("https://raw.githubusercontent.com/schwallergroup/ai4chem_course/main/notebooks/01%20-%20Basics/data/delaney-processed.csv" )

#Including only the columns that are relevant
dfsol = df[["Compound ID", "measured log solubility in mols per litre"]]
dfsmi = df[["Compound ID", "smiles"]]
dfmw = df[["Compound ID", "Molecular Weight"]]

merged = pd.merge(dfsol, dfsmi, on="Compound ID")
merged = pd.merge(merged, dfmw, on="Compound ID")

sns.regplot(x='Molecular Weight', y='measured log solubility in mols per litre', data=df)
plt.show()

# Compute z-scores for Molar Solubility
z_scores = np.abs(stats.zscore(merged['measured log solubility in mols per litre']))

# Common threshold: 3 standard deviations
# heavier outliers (higher MW)
# outliers = merged[(z_scores > 0.3) & (merged['Molecular Weight'] > 500)]

outliers = merged[(z_scores > 1.8) & (merged['measured log solubility in mols per litre'] > -4)]

sorted_outliers=outliers.sort_values("Molecular Weight", ascending=False)

print(sorted_outliers)
