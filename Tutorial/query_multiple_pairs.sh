#!/bin/bash

# Input file with protein-RNA pairs
input_file="protein_rna_pairs.txt"

# Output file for concatenated results
output_file="concatenated_output.csv"

# Create a temporary file
temp_file=$(mktemp)

# Loop over the protein-RNA pairs in the input file
while IFS=$'\t' read -r protein rna; do
    # Remove leading and trailing whitespace from protein and rna variables
    protein=$(echo "$protein" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    rna=$(echo "$rna" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')

    echo "$protein"
    echo "$rna"

    # Download the protein-RNA pair using curl
    curl "http://scrapid.tartaglialab.com/database?rna=${rna}&protein=${protein}" > "$temp_file"

    # If it's not the first iteration, remove the header from the temporary file
    if [ -f "$output_file" ]; then 
        tail -n +2 "$temp_file" >> "$output_file"
    else
        # Move the temporary file to the output file directly
        mv "$temp_file" "$output_file"
    fi

done < <(tr -d '\r' < "$input_file")

# Remove the temporary file
rm "$temp_file"
