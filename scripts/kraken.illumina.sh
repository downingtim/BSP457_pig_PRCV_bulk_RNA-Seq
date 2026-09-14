#!/bin/bash

# Directory containing the FASTP files
fastp_files_dir="FASTQS"
base_command="kraken2"
db_path="/home/share/"
threads=3
output_dir_base="KRAKEN_FILES"
output_dir_base2="KRAKEN_REPORT"

# Get a list of sample names by listing the files in the directory and extracting the base names _R1_001.fastq.gz
sample_names=$(ls $fastp_files_dir | grep -E '.fastq.gz$' | sed -E 's/_R([12])_001.fastq.gz$//' | sort | uniq)

# Iterate over sample names and construct the command for each
for sample_name in $sample_names; do
    # Construct output directories
    output_dir="${output_dir_base}/${sample_name}"
    output_dir2="${output_dir_base2}/${sample_name}"
    command="${base_command} --db ${db_path} --threads ${threads} --report ${output_dir2} --use-names --output ${output_dir} ${fastp_files_dir}/${sample_name}_R1_001.fastq.gz "
    echo "${command}"
done
