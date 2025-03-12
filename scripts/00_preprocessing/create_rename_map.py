import pandas as pd

# Specify the input and output file paths
input_table_path = "sequence_order.txt"  # Replace with the path to your input table
output_map_path = "rename_map.txt"    # Path to the output rename_map.txt

# Read the input table into a DataFrame
# Assuming the input table has two columns: original_filename and new_name separated by a tab
df = pd.read_csv(input_table_path, sep="\t", header=None, names=["original_filename", "new_name"])

# Process the filenames to extract the 8th and 9th attributes
rename_map = []
for _, row in df.iterrows():
    file_name = row["original_filename"]
    new_name = row["new_name"]
    
    # Split the filename into parts and extract the 8th and 9th attributes
    split = file_name.split('_')
    stem_name = f"{split[8]}_{split[9]}"  # Extract the stem name (8th and 9th attributes)
    
    # Append the processed data to the rename_map list
    rename_map.append((stem_name, new_name))

# Write the rename map to the output file
with open(output_map_path, "w") as outfile:
    for stem_name, new_name in rename_map:
        outfile.write(f"{stem_name}\t{new_name}\n")

print(f"Rename map saved to {output_map_path}")

