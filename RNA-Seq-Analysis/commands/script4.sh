#!/bin/bash

# Step 1: Extract first and last columns from test.bam.bam.txt and save to output.csv
# awk '{print $1 "," $NF}' test_trimmed.txt > output.csv

# Step 2: Loop through remaining txt files and extract last column, adding new column to output.csv
# for file in *.txt; do
#     if [ "$file" != "test_trimmed.txt" ]; then
#         awk '{print $NF}' "$file" | paste -d, output.csv - > output2.csv
#         mv output2.csv output.csv
#     fi
# done

# Set the directory
DIR="count_matrix"
# Step 1: Extract the first and last columns from the first file and save to output.csv
FIRST_FILE=$(ls ${DIR}/*.txt | head -n 1)
awk '{print $1 "," $NF}' $FIRST_FILE > ${DIR}/output.csv

# Step 2: Loop through the remaining txt files and extract the last column, adding it to output.csv
for file in ${DIR}/*.txt; do
    if [ "$file" != "$FIRST_FILE" ]; then
        awk '{print $NF}' "$file" | paste -d, ${DIR}/output.csv - > ${DIR}/output2.csv
        mv ${DIR}/output2.csv ${DIR}/output_tmp.csv
    fi
done

# Step 3: Remove folder names from the columns in output.csv
awk -F',' '{
    for(i=1; i<=NF; i++) {
        sub(".*/", "", $i)
    }
    print $0
}' ${DIR}/output_tmp.csv > ${DIR}/output.csv

# Step 4: Remove the temporary file
rm ${DIR}/output_tmp.csv