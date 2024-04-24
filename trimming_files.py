import pandas as pd

# Get the paths to the 2 files
file_path1 = input("type in the path to the arenosa file")
file_path2 = input("type in the path to the lyrata file")

# Read in the files
df1 = pd.read_csv(file_path1, sep='\t')
df2 = pd.read_csv(file_path2, sep='\t')

# Extract unique chromosome positions from lyrata file
unique_positions = df2[['CHROM', 'POS']]

# Filter rows of arenosa file to keep only those that are in teh lyrata file
df1_trimmed = df1.merge(unique_positions, on=['CHROM', 'POS'], how='inner')

# Save the trimmed DataFrame to a new file
output_file = 'arenosa_632_trimmed.txt'
df1_trimmed.to_csv(output_file, sep='\t', index=False)
