import pandas as pd

# Specify the columns to keep
columns_to_keep = ['CHROM', 'POS', 'AF', 'AC', 'AN'] 

# Get the file paths
file_path1 = input("type in the path to the trimmed arenosa file")
file_path2 = input("type in the path to the lyrata file")

# Read in the files
df1 = pd.read_csv(file_path1, sep='\t', usecols=columns_to_keep)
df2 = pd.read_csv(file_path2, sep='\t', usecols=columns_to_keep)

# Merge the two DataFrames on CHROM and POS columns
merged_df = pd.merge(df1, df2, on=['CHROM', 'POS'], how='inner', suffixes=('_arenosa', '_lyrata'))

# Calculate the difference between AFs
merged_df['AF_difference'] = abs(merged_df['AF_lyrata'] - merged_df['AF_arenosa'])

# Sort the absolute differences in descending order
sorted_diff = merged_df.sort_values(by='AF_difference', ascending=False)

# Count the number of differences above a certain threshold
threshold = 0.8  # Adjust the threshold as needed
significant_differences = (sorted_diff['AF_difference'] > threshold).sum()

# Print the most different AFs
print("Most different AFs:")
print(sorted_diff.head())

# Print results to a file
sorted_diff.to_csv('significant_differences.csv', index=False)


# Print the total number of differences above the threshold
print("Total significant differences above threshold:", significant_differences)
