# Fine-scale population structure in Alberta's westslope cutthroat trout (_Oncorhynchus lewisi_)
**The full study was submitted to a peer-reviewed journal on Jan 15, 2025.**


## Analysis pipeline
- prepare_data_forfilter.sh
  - bash script to extract data with 25% missingness across both sites and individuals
  - uses a minor allele count threshold of at least 3
  - the output is ready for further filtering in R
  - requires [bcftools](https://samtools.github.io/bcftools/) and [vcftools](https://vcftools.github.io/index.html)
- filter.R
  - R script to filter on heterozygosity and prune based on linkage disequilibrium
  - outputs a file of loci to keep, "filtered.txt"
  - requires [dartR](https://green-striped-gecko.github.io/dartR/) and [tidyverse](https://www.tidyverse.org/)
- create_inputs.sh
  - bash script that filters data based on the file "filtered.txt" and outputs both .vcf.gz and .geno formatted files
  - requires [bcftools](https://samtools.github.io/bcftools/)
- make_geno.sh
  - bash script that extracts individuals from a vcf based on a keeplist of sample names
  - outputs both .vcf.gz and .geno formatted files
  - requires [bcftools](https://samtools.github.io/bcftools/)
- snmf.all.R
  - R script that uses .geno formatted files to calculate ancestry coefficients for differe K values, with 100 replicates
  - output plots are written into pdf files
  - requires [LEA](http://membres-timc.imag.fr/Olivier.Francois/LEA/) and [tidyverse](https://www.tidyverse.org/)
- seven scripts with names formatted as pca_<dataset>.R
  - R scripts that calculate principal components data and assign samples to user-defined groups to color the output plots
  - output plots are written into pdf files
  - requires [vcfR](https://github.com/knausb/vcfR), [adegenet](https://adegenet.r-forge.r-project.org/), [adegraphics](https://cran.r-project.org/web/packages/adegraphics/vignettes/adegraphics.html) and [tidyverse](https://www.tidyverse.org/)
- popstat.R
  - R script that calculated various population statistics based on an input .pdf and population assignment file
  - requires the file to assign sample to both sampling sites and populations (two different columns)
  - all statistics (heterozygosity, Fis, effective population size, Fst) are calculated both for sites and for populations
  - output is written into tables, in Rstudio generates plots as well
  - requires [dartR](https://green-striped-gecko.github.io/dartR/), [strataG](https://github.com/EricArcher/strataG) and [tidyverse](https://www.tidyverse.org/)



## Software versions used by the scripts
  - [adegenet](https://adegenet.r-forge.r-project.org/): 2.1.10
  - [adegraphics](https://cran.r-project.org/web/packages/adegraphics/vignettes/adegraphics.html): 1.0-21
  - [bcftools](https://mafft.cbrc.jp/): 1.13
  - [dartR](https://green-striped-gecko.github.io/dartR/): 2.9.8
  - [LEA](http://membres-timc.imag.fr/Olivier.Francois/LEA/): 3.16.0
  - [strataG](https://github.com/EricArcher/strataG): 2.5.01
  - [tidyverse](https://www.tidyverse.org/): 2.0.0
  - [vcfR](https://github.com/knausb/vcfR): 1.15.0
  - [vcftools](https://vcftools.github.io/index.html): 0.1.16
