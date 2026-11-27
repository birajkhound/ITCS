# Article title: Identification of Thermal Class-Specific Genomic and Proteomic Compositional Signatures in Bacteria Using a Weighted Ensemble Feature Ranking Framework, 
# Author Names: Biraj Khound1, Deb K Boruah1, Pranjal Phukan1, Kalyan Sarma1, Prince Das1, Robita Doley2, Alicia Passah1;
# Affiliation: 
# 1 Department of Diagnostic and Interventional Radiology, All India Institute of Medical Sciences, Changsari, Guwahati, 781101, Assam, India.
# 2 Department of Molecular Biology and Biotechnology, Tezpur University, Napaam, Tezpur, 784028, Assam, India.
# E-mail address of the corresponding author: prsbirajaiimsguwahati@gmail.com


import pandas as pd

# ==========================================================
# USER INPUT
# ==========================================================
files = ["D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\Details of Bacteria Considered for Compositional Bias Study.xlsx",
         "D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\AAB_values.xlsx",
         "D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\AAPB_values.xlsx",
         "D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\CPB_values.xlsx",
         "D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\CUB_values.xlsx",
         "D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\DUB_values.xlsx"]
# for i in range(1, 7):
#     path = input(f"Enter path of Excel file {i}: ").strip()
#     files.append(path)

common_col = "zip file name"

output_file = "D:\\TEMPURA based Dataset\\featured dataset preparation\\combined dataset\\complete_featured_dataset.xlsx"

# ==========================================================
# READ FILES
# ==========================================================
dfs = []
for i, file in enumerate(files, start=1):
    df = pd.read_excel(file)

    if common_col not in df.columns:
        raise ValueError(f"{common_col} not found in File {i}")

    # Convert common column to string for safe matching
    df[common_col] = df[common_col].astype(str).str.strip()

    # Rename all columns except common column
    new_cols = {}
    for col in df.columns:
        if col != common_col:
            new_cols[col] = f"File{i}_{col}"

    df = df.rename(columns=new_cols)

    dfs.append(df)

# ==========================================================
# FIND COMMON VALUES IN ALL SIX FILES
# ==========================================================
common_values = set(dfs[0][common_col])

for df in dfs[1:]:
    common_values = common_values.intersection(set(df[common_col]))

# Keep only common rows
filtered_dfs = []
for df in dfs:
    filtered_dfs.append(df[df[common_col].isin(common_values)])

# ==========================================================
# MERGE ALL FILES ON COMMON COLUMN
# ==========================================================
result = filtered_dfs[0]

for df in filtered_dfs[1:]:
    result = pd.merge(result, df, on=common_col, how="inner")

# ==========================================================
# SAVE OUTPUT
# ==========================================================
result.to_excel(output_file, index=False)

print(f"Done! Output saved as {output_file}")