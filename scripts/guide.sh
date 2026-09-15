#!/bin/bash

# exectution: ../KRAKEN_VALID_FILES/Porcine_traheal_ALI_cultures-Animal-4100-24hpi-infected_S37.fastq
# ls ../KRAKEN_VALID_FILES/* | \
# perl -e 'while(<>){ chomp; $_=~ s/.fastq//g; $_=~s/..\/KRAKEN_VALID_FILES\///g; print "sh guide.sh $_ &> REPORT/$_.txt & \n";}'

sample=$1
fasta="DQ811787.fasta"

set -e
check_command() {
    if [ $? -ne 0 ]; then
        echo "Error: Command failed - $1"
        exit 1
    fi
}

check_file() {
    if [[ ! -f "$1" || ! -s "$1" ]]; then
        echo "Error: File $1 is missing or empty"
        exit 1
    fi
}
mkdir -p BAM_FILES SAM_FILES BAM_FINAL COVERAGE ERROR_FILES FB_VCF_FILES BCF_VCF_FILES CONSENSUS REPORT DEPTH

check_file "$fasta"
ls -lt $fasta 
ls -lt      ../KRAKEN_VALID_FILES/${sample}.fastq
check_file "../KRAKEN_VALID_FILES/${sample}.fastq"
echo "Processing sample: $sample"

# Mapping Illumina reads with validation
echo "Mapping Illumina reads..."
echo "minimap2 -ax sr $fasta  ../KRAKEN_VALID_FILES/${sample}.fastq  > SAM_FILES/${sample}.${fasta}.illumina.sam" 
minimap2 -ax sr $fasta    "../KRAKEN_VALID_FILES/${sample}.fastq"  \
	 > "SAM_FILES/${sample}.${fasta}.illumina.sam" 2> "ERROR_FILES/${sample}.illumina.sam.errors.txt"
check_command "Illumina mapping"
check_file "SAM_FILES/${sample}.${fasta}.illumina.sam"
ls -lt "SAM_FILES/${sample}.${fasta}.illumina.sam"

samtools view -bS "SAM_FILES/${sample}.${fasta}.illumina.sam" > "BAM_FILES/${sample}.${fasta}.illumina.bam"
check_file "BAM_FILES/${sample}.${fasta}.illumina.bam"
rm -rf "SAM_FILES/${sample}.${fasta}.illumina.sam" 
# Sort by queryname first
echo "Sorting Illumina BAM by queryname..."
samtools sort -n -o "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam" "BAM_FILES/${sample}.${fasta}.illumina.bam"
check_file "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam"

# Fix mate information
echo "Fixing mate information for Illumina reads..."
samtools fixmate -m "BAM_FILES/${sample}.${fasta}.illumina.queryname.bam" "BAM_FILES/${sample}.${fasta}.illumina.fixmate.bam"
check_file "BAM_FILES/${sample}.${fasta}.illumina.fixmate.bam"

# Sort by coordinate for markdup
echo "Sorting Illumina BAM by coordinate..."
samtools sort -o "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam" "BAM_FILES/${sample}.${fasta}.illumina.fixmate.bam"
check_file "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam"

# Mark duplicates
echo "Marking duplicates in Illumina reads..."
samtools markdup -r -s "BAM_FILES/${sample}.${fasta}.illumina.coordsort.bam" "BAM_FINAL/${sample}.${fasta}.merged.bam"
check_file "BAM_FINAL/${sample}.${fasta}.merged.bam"
samtools index "BAM_FINAL/${sample}.${fasta}.merged.bam"
samtools flagstat "BAM_FINAL/${sample}.${fasta}.merged.bam" > "BAM_FINAL/${sample}.${fasta}.flagstat"
samtools coverage "BAM_FINAL/${sample}.${fasta}.merged.bam" > "COVERAGE/${sample}.${fasta}.coverage.txt"
samtools depth "BAM_FINAL/${sample}.${fasta}.merged.bam" > "DEPTH/${sample}.${fasta}.depth.txt"

# Variant calling with freebayes
echo "Calling variants with freebayes..."
/mnt/lustre/RDS-live/downing/freebayes/build/freebayes -f $fasta -F 0.01 -p 1   --min-alternate-count 1 --min-alternate-fraction 0.001 \
	  "BAM_FINAL/${sample}.${fasta}.merged.bam" > "FB_VCF_FILES/${sample}.${fasta}.fb.vcf" \
	  2> "ERROR_FILES/${sample}.${fasta}.fb.errors.txt"
bgzip -f "FB_VCF_FILES/${sample}.${fasta}.fb.vcf"
tabix -f -p vcf "FB_VCF_FILES/${sample}.${fasta}.fb.vcf.gz"

# Normalize variants
bcftools norm -c w -f $fasta -m-both -Oz -o "FB_VCF_FILES/${sample}.${fasta}.norm.fb.vcf.gz" \
	 "FB_VCF_FILES/${sample}.${fasta}.fb.vcf.gz"
tabix -f -p vcf "FB_VCF_FILES/${sample}.${fasta}.norm.fb.vcf.gz"

# BCFtools variant calling
echo "Calling variants with BCFtools..."
bcftools mpileup -A -Ob -f $fasta "BAM_FINAL/${sample}.${fasta}.merged.bam" |\
    bcftools call -cvO v --ploidy 1 -o "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf"
bgzip -f "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf"
tabix -f -p vcf "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf.gz"

# Normalize BCF variants
bcftools norm -c w -f $fasta -m-both -Oz -o "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz" \
	 "BCF_VCF_FILES/${sample}.${fasta}.bcf.vcf.gz"
tabix -f -p vcf "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz"

# Generate consensus
echo "Generating consensus sequences..."
bcftools consensus -H LA -f $fasta     "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz" \
    -o "CONSENSUS/${sample}.${fasta}.LA.fasta"
bcftools consensus -H LR -f $fasta     "BCF_VCF_FILES/${sample}.${fasta}.norm.bcf.vcf.gz" \
    -o "CONSENSUS/${sample}.${fasta}.LR.fasta"
