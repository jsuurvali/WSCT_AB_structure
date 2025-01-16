#!/bin/bash

# remove Windows-specific characters if filtered list was created in windows, the add a tab to the end of each row
cat filtered.txt | tr -d '\r' | awk '{print $0"\t"}' - > tmp.varlist


# extract sites from the vcf based on the filter list, save as a new vcf
bcftools view -h wsct.snps.var.mac3.miss25.vcf.gz > tmp.header
bcftools view wsct.snps.var.mac3.miss25.vcf.gz | \
	grep -f tmp.varlist | \
	cat tmp.header - | \
	bcftools view -Oz -o wsct.filtered.all.vcf.gz
tabix -p vcf wsct.filtered.all.vcf.gz
rm tmp.varlist

# convert output to geno format
bcftools view -H wsct.filtered.all.vcf.gz | cut -f10- | sed -re 's|0/1|1|g' -e 's|1/0|1|g' -e 's|0/0|0|g' -e 's|1/1|2|g' -e 's|\./\.|9|g' | tr -d '\t' > wsct.filtered.all.geno
