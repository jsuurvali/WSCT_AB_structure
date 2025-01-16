#!/bin/bash
# usage: sh make_geno.sh <subset_name>

# subsets data based on an input keeplist, saves it in .vcf.gz and .geno formats.
# keeplist is a file containing a list of individuals to keep
# keeplist has to be formatted as keeplist.miss25.<subset_name>.sh

bcftools view -S keeplist.miss25.${1} wsct.filtered.all.vcf.gz -Oz -o wsct.filtered.${i}.vcf.gz
bcftools view -H wsct.filtered.${1}.vcf.gz | cut -f10- | sed -re 's|0/1|1|g' -e 's|1/0|1|g' -e 's|0/0|0|g' -e 's|1/1|2|g' -e 's|\./\.|9|g' | tr -d '\t' > wsct.filtered.${1}.geno
