# drug_discovery_capstone_full.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

plots_folder = "plots"
os.makedirs(plots_folder, exist_ok=True)

#1.Setup Paths and Folders

data_folder = "data"
results_folder = "results"
os.makedirs(results_folder, exist_ok=True)

chembl_file = os.path.join(data_folder, "chembl_raw.csv")
pubchem_file = os.path.join(data_folder, "pubchem_raw.csv")

#2. Load and Clean ChEMBL CSV
chembl = pd.read_csv(chembl_file, sep=';', on_bad_lines="skip")  #skips malformed rows

# Clean column names
chembl.columns = chembl.columns.str.strip().str.replace('"','').str.replace(' ', '_').str.lower() #remove spaces and quotes

# Keep necessary columns
keep_cols = ['molecule_chembl_id', 'standard_value', 'standard_units']
chembl = chembl[[col for col in keep_cols if col in chembl.columns]]

# Rename columns
chembl.rename(columns={
    'molecule_chembl_id':'compound_id',
    'standard_value':'IC50_nM',
    'standard_units':'units'
}, inplace=True)

# Filter units if exists to keep measurements in nM only
if 'units' in chembl.columns:
    chembl['units'] = chembl['units'].astype(str)
    chembl = chembl[chembl['units'].str.upper() == 'NM']

# Convert IC50 to numeric
chembl['IC50_nM'] = pd.to_numeric(chembl['IC50_nM'], errors='coerce')
chembl = chembl.dropna(subset=['IC50_nM'])

# Convert to uM and pIC50
chembl['IC50_uM'] = chembl['IC50_nM'] / 1000
chembl['pIC50'] = -np.log10(chembl['IC50_uM'])
chembl['source'] = "ChEMBL"

chembl.to_csv(os.path.join(results_folder,"chembl_clean.csv"), index=False)
print("✅ ChEMBL cleaned. Top 5 compounds:")
print(chembl.sort_values('pIC50', ascending=False).head())

#3. Load and Clean PubChem CSV
pubchem = pd.read_csv(pubchem_file, sep=None, engine='python', on_bad_lines="skip")  #allows automatic sepearator detection and skip bad rows

# Clean column names
pubchem.columns = pubchem.columns.str.strip().str.replace('"','').str.replace(' ', '_').str.lower()

# Detect columns to try to automatically find the right columns
compound_col = [c for c in pubchem.columns if 'cid' in c.lower() or 'compound' in c.lower()]
ic50_col = [c for c in pubchem.columns if 'ic50' in c.lower() or 'value' in c.lower()]

if not compound_col or not ic50_col:
    raise Exception("Cannot detect PubChem compound/IC50 columns. Columns found:", pubchem.columns.tolist())

compound_col = compound_col[0]
ic50_col = ic50_col[0]

pubchem = pubchem[[compound_col, ic50_col]]
pubchem.rename(columns={compound_col:'compound_id', ic50_col:'IC50_nM'}, inplace=True)

# Convert IC50 to numeric
pubchem['IC50_nM'] = pd.to_numeric(pubchem['IC50_nM'], errors='coerce')
pubchem = pubchem.dropna(subset=['IC50_nM'])

# Convert to uM and pIC50
pubchem['IC50_uM'] = pubchem['IC50_nM'] / 1000
pubchem['pIC50'] = -np.log10(pubchem['IC50_uM'])
pubchem['source'] = "PubChem"

pubchem.to_csv(os.path.join(results_folder,"pubchem_clean.csv"), index=False)
print("✅ PubChem cleaned. Top 5 compounds:")
print(pubchem.sort_values('pIC50', ascending=False).head())

# ---------------------------
# 4️⃣ Combine Datasets
# ---------------------------
combined = pd.concat([
    chembl[['compound_id','IC50_nM','IC50_uM','pIC50','source']],
    pubchem[['compound_id','IC50_nM','IC50_uM','pIC50','source']]
], ignore_index=True)

combined.to_csv(os.path.join(results_folder,"combined_hits.csv"), index=False)
print("✅ Combined ChEMBL + PubChem dataset saved.")

# ---------------------------
# 5️⃣ Select Top Candidates
# ---------------------------
top_candidates = combined.sort_values('pIC50', ascending=False).head(20) #picks top 20 compounds with highest pIC50

# Add interpretation
def interpret_pIC50(p):
    if p >= 7: return "Very potent"
    elif p >= 6: return "Potent"
    elif p >= 5: return "Moderate"
    else: return "Weak"

top_candidates['Interpretation'] = top_candidates['pIC50'].apply(interpret_pIC50)
top_candidates.to_csv(os.path.join(results_folder,"top_candidates.csv"), index=False)

print("✅ Top candidates saved to 'results/top_candidates.csv'")
print(top_candidates[['compound_id','source','IC50_nM','pIC50','Interpretation']])


# 1️⃣ Histogram of pIC50 (all compounds)
# ---------------------------
plt.figure(figsize=(8,5))
plt.hist(combined['pIC50'], bins=30, color='teal', edgecolor='black')
plt.title("Distribution of pIC50 Values")
plt.xlabel("pIC50")
plt.ylabel("Number of Compounds")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(plots_folder,"pIC50_distribution.png"))
plt.show()

# ---------------------------
# 2️⃣ Top 20 Compounds Bar Plot
# ---------------------------



# Clean compound IDs for plotting
top_candidates['compound_id'] = top_candidates['compound_id'].astype(str).str.strip().str.replace('"','')

plt.figure(figsize=(10,6))
plt.barh(top_candidates['compound_id'], top_candidates['pIC50'], color='darkorange')
plt.xlabel("pIC50")
plt.ylabel("Compound ID")
plt.title("Top 20 Compounds by pIC50")
plt.gca().invert_yaxis()  # highest pIC50 on top
plt.tight_layout()
plt.savefig(os.path.join(plots_folder,"top20_compounds.png"))
plt.show()



# ---------------------------
# 3️⃣ Source Comparison Boxplot
# ---------------------------
plt.figure(figsize=(7,5))
sources = combined['source'].unique()
data_to_plot = [combined[combined['source']==s]['pIC50'] for s in sources]
plt.boxplot(data_to_plot, labels=sources)
plt.ylabel("pIC50")
plt.title("pIC50 Distribution by Source")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(plots_folder,"source_comparison.png"))
plt.show()


