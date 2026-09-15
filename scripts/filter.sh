#!/bin/bash

set -euo pipefail

FASTA="DQ811787.fasta"

GENOME_LENGTH=$(grep -v '^>' "$FASTA" | tr -d '\n' | wc -c)

START=500
END=$((GENOME_LENGTH-500))

mkdir -p FILTERED_VCF SNP_TABLES

for FB in FB_VCF_FILES/*.norm.fb.vcf.gz
do
    SAMPLE=$(basename "$FB" .norm.fb.vcf.gz)
    BCF="BCF_VCF_FILES/${SAMPLE}.norm.bcf.vcf.gz"

    echo "Processing ${SAMPLE}"

    ####################################################################
    # FREEBAYES FILTERS
    ####################################################################

    bcftools view \
        -v snps \
        -i 'QUAL>=20 && INFO/DP>=5 && INFO/AO>=5 && INFO/AF>0.05' \
        "$FB" |
    awk -v start="$START" -v end="$END" '
        /^#/ {print; next}
        $2 > start && $2 < end
    ' | bgzip -c \
    > FILTERED_VCF/${SAMPLE}.freebayes.filtered.vcf.gz

    tabix -f -p vcf \
        FILTERED_VCF/${SAMPLE}.freebayes.filtered.vcf.gz

    ####################################################################
    # BCFTOOLS FILTERS
    ####################################################################

    (
    bcftools view -h "$BCF"

    bcftools view -H "$BCF" | \
    awk -F'\t' '
    {
        qual=$6

        split($8,info,";")

        dp=0
        altreads=0

        for(i in info)
        {
            if(info[i] ~ /^DP=/)
            {
                split(info[i],a,"=")
                dp=a[2]
            }

            if(info[i] ~ /^DP4=/)
            {
                split(info[i],a,"=")
                split(a[2],d,",")

                altreads=d[3]+d[4]
            }
        }

        if(length($4)==1 &&
           length($5)==1 &&
           qual>=20 &&
           dp>=5 &&
           altreads>=5 &&
           $2>500 &&
           $2<'$END')
        {
            print
        }
    }'
    ) | bgzip -c \
    > FILTERED_VCF/${SAMPLE}.bcftools.filtered.vcf.gz

    tabix -f -p vcf \
        FILTERED_VCF/${SAMPLE}.bcftools.filtered.vcf.gz

    ####################################################################
    # UNION OF BOTH CALLERS
    ####################################################################

    (
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\tFB\n' \
        FILTERED_VCF/${SAMPLE}.freebayes.filtered.vcf.gz

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\tBCF\n' \
        FILTERED_VCF/${SAMPLE}.bcftools.filtered.vcf.gz
    ) | \
    sort -k1,1 -k2,2n -u \
    > SNP_TABLES/${SAMPLE}.merged_snps.tsv

    echo "${SAMPLE}"
    echo -n "  FreeBayes SNPs : "
    bcftools view -H FILTERED_VCF/${SAMPLE}.freebayes.filtered.vcf.gz | wc -l

    echo -n "  BCFtools SNPs  : "
    bcftools view -H FILTERED_VCF/${SAMPLE}.bcftools.filtered.vcf.gz | wc -l

    echo -n "  Union SNPs     : "
    wc -l < SNP_TABLES/${SAMPLE}.merged_snps.tsv
done
