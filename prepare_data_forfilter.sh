#!/bin/bash

# extract SNPs with minor allele count of at least 3
bcftools view -m2 -M2 -v snps -e 'MAC < 3' wsct.merged.var.vcf.gz -Oz -o wsct.snps.var.mac3.vcf.gz

# calculate the percentage of missing sites per individual
vcftools --gzvcf wsct.snps.var.mac3.vcf.gz --missing-indv --out tmp

# extract individuals with less than 25% missing data, then re-filter as before
# allow only sites with no more than 25% of all individuals having missing data

awk '$5 < 0.25 {print $1}' tmp.imiss | \
        bcftools view -S - wsct.snps.var.mac3.vcf.gz | \
        bcftools view -m2 -M2 -v snps -e 'MAC < 3 | F_MISSING > 0.25' -Oz -o wsct.snps.var.mac3.miss25.vcf.gz

bcftools query -l wsct.snps.var.mac3.miss25.vcf.gz > keeplist.miss25

rm tmp*

cut -f1 sample30.txt | bcftools view -S - wsct.snps.var.mac3.miss25.vcf.gz -Oz -o forfilter.vcf.gz
