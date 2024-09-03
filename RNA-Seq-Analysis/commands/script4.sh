#!/bin/bash

# Step 1: Extract first and last columns from test.bam.bam.txt and save to output.csv
awk '{print $1 "," $NF}' test_trimmed.txt > output.csv

# Step 2: Loop through remaining txt files and extract last column, adding new column to output.csv
for file in *.txt; do
    if [ "$file" != "test_trimmed.txt" ]; then
        awk '{print $NF}' "$file" | paste -d, output.csv - > output2.csv
        mv output2.csv output.csv
    fi
done

