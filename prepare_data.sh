#!/bin/bash

# ~~~~~ STEP 1 ~~~~~
# remove info and format fields that will not be used in the next steps anyway.
# the additional input file samples.list contains names of samples to keep, one per row

bcftools annotate -x 'INFO,^FORMAT/GT' AlbertaR1_2.vcf.gz | bcftools view -S samples.list --force-samples -Oz -o Alberta_R1.noannot.vcf.gz
bcftools annotate -x 'INFO,^FORMAT/GT' Alberta_2018_WCT_2ndRun_2_2.vcf.gz | bcftools view -S samples.list --force-samples -Oz -o Alberta_2018.noannot.vcf.gz

# ~~~~~ STEP 2 ~~~~~
# Extract chr, bp, ALT and REF from the files generated in previous step. Merge the outputs on chr and bp.
# Filter by removing any site that does not have the same REF and ALT alleles in both files. Also remove MNVs, multiallelics and indels.

bcftools view Alberta_R1.noannot.vcf.gz | awk 'BEGIN {FS=OFS="\t"} $0 !~/#/ {print $1":"$2,NR,$4,$5}' | sort -k1,1 > refalt1
bcftools view Alberta_2018.noannot.vcf.gz | awk 'BEGIN {FS=OFS="\t"} $0 !~/#/ {print $1":"$2,NR,$4,$5}' | sort -k1,1 > refalt2
join refalt1 refalt2 | \
	awk '$3 == $6 && $4 == $7 && ($3 == "A" || $3 == "C" ||$3 == "G" || $3 == "T") && ($4 == "A" || $4 == "C" || $4 == "G" || $4 == "T") {print $2"\t"$1}' | \
	sort -k1,1n | awk 'BEGIN {FS="[:\t]" ; print "#"} {print $2"\t"$3"\t"}' > okSNPs.forMerge
rm refalt*

# ~~~~~ STEP 3 ~~~~~
# Extract the snps in okSNPs.forMerge from the files, then write into new vcfs and index those.

bcftools view --no-version Alberta_2018.noannot.vcf.gz | grep -f okSNPs.forMerge | bcftools view -Oz -o Alberta_2018.okSNPs.vcf.gz
tabix -p vcf Alberta_2018.okSNPs.vcf.gz

bcftools view --no-version Alberta_R1.noannot.vcf.gz | grep -f okSNPs.forMerge | bcftools view -Oz -o Alberta_R1.okSNPs.vcf.gz
tabix -p vcf Alberta_R1.okSNPs.vcf.gz

# ~~~~~ STEP 4 ~~~~~
# Merge the files, order samples by the samples.list file that was also used above

bcftools merge Alberta_2018.okSNPs.vcf.gz Alberta_R1.okSNPs.vcf.gz | bcftools view -S samples.list -Oz -o WSCT.okSNPs.vcf.gz
tabix -p vcf WSCT.okSNPs.vcf.gz

# ~~~~~ STEP 5 ~~~~~
# Remove sites with 30% or more individuals having missing data
# Then remove individuals with 25% or more missing data
# Next, remove sites with 25% or more individuals having missing data, or with MAC less than 10
# Finally, minimize the impact of LD by thinning the vcf so that each snp is no closer than 1kb to any other SNP

bcftools view -i 'F_MISSING < 0.3' WSCT.okSNPs.vcf.gz -Oz -o tmp.vcf.gz
vcftools --missing-indv --gzvcf tmp.vcf.gz --out tmp

awk '$5 < 0.25 {print $1}' tmp.imiss > keeplist.samples.2
bcftools view -S keeplist.samples.2 WSCT.okSNPs.vcf.gz | bcftools view -i 'F_MISSING < 0.25 & AC <= AN-10 & AC >= 10' > tmp2.vcf
vcftools --gzvcf tmp2.vcf --thin 1000 --recode --out tmp.WSCT.okSNPs.miss25.mac10.thin1k
bcftools view tmp.WSCT.okSNPs.miss25.mac10.thin1k.recode.vcf -Oz -o WSCT.okSNPs.miss25.mac10.thin1k.vcf.gz

rm tmp*

# ~~~~~ STEP 6 ~~~~~
# Create a .geno file for snmf to use

bcftools view -H WSCT.okSNPs.miss25.mac10.thin1k.vcf.gz | cut -f10- | tr -d '/' | sed -re 's/00/0/g' -e 's/11/2/g' -e 's/(01|10)/1/g' -e 's/\.\./9/g' | tr -d '\t' > WSCT.okSNPs.miss25.mac10.thin1k.geno
