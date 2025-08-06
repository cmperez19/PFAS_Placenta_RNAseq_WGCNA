#pip install pandas
import os
import pandas as pd

# Directory containing your text files
directory = '/projects/eversongroup/glowing/transcriptcounts/counts'

# Initialize an empty dataframe to hold all the extracted columns
combined_df = pd.DataFrame()

# Function to extract the column label from the filename
def extract_label(filename):
    parts = filename.split('_')
    if len(parts) > 2:
        return '_'.join(parts[:2])
    return filename

# Loop through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('_FC.txt'):  # Ensure it is a text file
        filepath = os.path.join(directory, filename)
        
        # Read the file into a dataframe
        df = pd.read_csv(filepath, sep='\t', skiprows=1, index_col='Geneid')  # Set 'Geneid' as index
        
        # Extract the last column and set its name
        last_column = df.iloc[:, -1]
        last_column.name = extract_label(filename)
        
        # If combined_df is empty, add the row names (index) first
        if combined_df.empty:
            combined_df = last_column.to_frame()
        else:
            combined_df = combined_df.join(last_column)

# Save the combined dataframe to a CSV file
combined_df.to_csv('unfiltered_glowing_counts.csv')

