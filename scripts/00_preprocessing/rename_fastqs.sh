#!/bin/bash

# Default directory (can be overridden by an argument)
file_dir="../fastq"

# If an argument is provided, use it as the directory
if [[ -n $1 ]]; then
  file_dir="$1"
fi

# Read the mapping file (adjust path to where `rename_map.txt` is stored)
rename_map_file="rename_map.txt"

# Read the mapping file line by line
while read -r pattern replacement; do
  # Rename files matching the pattern
  for file in "$file_dir"/*"$pattern"*; do
    if [[ -f $file ]]; then # Check if it's a file
      mv "$file" "${file/$pattern/$replacement}"
    fi
  done
done < "$rename_map_file"

