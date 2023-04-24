# Load the standard libraries
import numpy as np
import pandas as pd

# read in excel file 
apple_dfs = pd.read_excel('../data/GSE182822_Matrix_FPKM.xlsx', sheet_name=None)

# concatenate all DataFrames into one and add a 'Sheet' column
apple_df = pd.concat([df.assign(Apple_Type=sheet) for sheet, df in apple_dfs.items()], ignore_index=True)

# Move apple name column to the first position
last_col = apple_df.pop('Apple_Type')
apple_df.insert(0, 'Apple_Type', last_col)

# Split the values in the gene_id column by the '|' character
new_cols = apple_df['GENE_ID'].str.split('|', expand=True)

# Assign the new column names to the DataFrame
new_cols.columns = ['col1', 'Gene_ID', 'col3', 'Gene_Ref', 'col5']
new_cols = new_cols.drop(['col1','col3', 'col5'], axis=1)

apple_df = pd.concat([apple_df.iloc[:, :1], new_cols, apple_df.iloc[:, 1:]], axis=1)

# make new df with only genes shared by all three apple types

# Group dataframe by Gene_ID and count the number of unique Apple_Type values for each group
counts = apple_df.groupby('Gene_ID')['Apple_Type'].nunique()

# Get the Gene_IDs that are present in all three apple types
genes_in_all_types = counts[counts == 3].index.tolist()

# Filter the original dataframe to keep only the rows corresponding to these genes
shared_genes_df = apple_df[apple_df['Gene_ID'].isin(genes_in_all_types)]

# Save the DataFrame to a CSV file
shared_genes_df.to_csv(r'../data/apple_shared_genes.csv', index=False)
