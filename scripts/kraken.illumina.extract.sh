#!/bin/bash

fastp_files_dir="FASTQS"

# Base command components - Porcine_traheal_ALI_cultures-Animal-4100-6hpi-mock_S26_R1_001.fastq.gz
base_command="python /mnt/lustre/RDS-live/downing/KrakenTools/extract_kraken_reads.py -k "
output_dir_base="KRAKEN_FILES"
output_dir_base2="KRAKEN_REPORT"
sample_names=$(ls $fastp_files_dir | grep -E '_001.fastq.gz$' | sed -E 's/_R1_001.fastq.gz$//' | sort | uniq)

for sample_name in $sample_names; do
    output_dir="${output_dir_base}/${sample_name}"
    output_dir2="${output_dir_base2}/${sample_name}"
    out_file=${sample_name}

    taxon=11146
    command="${base_command} $output_dir_base/${sample_name} -s1 $fastp_files_dir/${sample_name}_R1_001.fastq.gz -o KRAKEN_VALID_FILES/${out_file}.fastq -t $taxon --fastq-output -r $output_dir_base2/${sample_name} --include-children --include-parents "
    echo "${command}"

    taxon=9822 # pig

done
