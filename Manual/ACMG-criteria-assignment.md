# ACMG Criteria Assignment

## Pathogenic Criteria

### PVS1
Will be assigned to a variant if it is a stopgain or frameshift deletion/insertion located 50 bp further than the end of the final exon. (Based on the ExonicFunc column anf the REK_canon library) 
Will be assigned to a splicing variant with a dbscSNV score of more than 0.6 (ADA or RF) (Based on the Func column and the dbscSNV score annotation)

### PS1
Will be assigned if a variant have the same AA ref and AA alt as a known pathogenic variant.
Using all known pathogenic variants from clinvar

### PS2
Will be assigned if a variant is assumed de novo and parents are disease free. This requires trio data. See the [dedicated wiki page](https://github.com/a-xavier/tapes/wiki/Trio) for more informations on trio data.

### PS3
Will be assigned if clinvar classifies the variant as Pathogenic or drug reponse and the level of evidence is either ‘practice guideline’ or ‘reviewed by expert panel’

### PS4
Will be assigned if a variant is enriched in the samples provided. Requires either ‘output_with_samples.csv’ from the annotation to keep sample genotyping data or an annotated multi-sample vcf. PS4 will take the affected individuals with the mutations and the total number of individuals in the disease cohort and compare it to the data from gnomad_genome and gnomad_exome.
The number of individuals with and without variants in public data is extrapolated with the following formula:
Minor allele frequency in control population (MAF) = MAF_c=y×〖10〗^(-x)
Number of individuals with the variant in control population = n_c=⌈y⌉
Total number of individuals in control population = N_c=〖10〗^x/2-n_c
Then a fisher’s exact test is performed to calculate the odd ratios, the confidence interval and the p value.
PS4 will only be considered if at least 2 samples are affected by a variant. Otherwise, Intervar PS4 database, based on GWAS database will be used.
PS4 will be assigned if the Odd Ratio is superior to 20, the confidence interval does not cross one and the p value is under 0.01

### PM1
Will be assigned if the variant is a Missense variant (nonsynonymous SNV) and is located in a in a domain without benign variants (Using Intervar db) for benign domains

### PM2
Will be assigned if the variant is in a recessive gene and has a frequency under 0.005 or is in a dominant gene and has no frequency data available. Recessive and Dominant/Haploinsufficient genes were infered using Pli and Prec scores computed by Lek et al, 2016. A gene is considered dominant dominant with a pli >0.85 and recessive if prec >0.85

### PM4
Will be assigned if the variant is an in-frame deletion/insertion in a non-repeat region of the gene. Using the repeat_dict database.

### PM5
Will be assigned if a variant have the same AA ref and  a different AA alt as a known pathogenic variant.
Using all known pathogenic variants from clinvar

### PP2
Will be assigned if the variant is Missense (nonsynonymous SNV) in a gene where missense variants represents at least 80 percent of all known pathogenic variants (using PP2_BP1 database)

### PP3
Will be assigned if the variant is predicted to be pathogenic using various in-silico prediction tools (sift, lrt, mutationtaster, mutation assessor, fathmm, provean, meta svm, meta lr, mcap, mkl, genocanyon, gerp). Each prediction tool will add +1 or -1. A score over 3 will assign PP3 to a variant.

### PP5
Will be assigned the variant is classified as pathogenic or likely pathogenic by clinvar but the evidence is limited.

## Benign criteria

### BA1
Will be assigned to a variant if its frequency in gnomad_exome/exac or gnomad_genome is superior to 0.05

### BS1
Will be assigned to a variant if its frequency is superior to a cutoff (0.005) for a rare disease.

### BS2
Will be assigned if the variant was observed in a healthy individual as homozygous for a recessive disease and heterozygous for a dominant disease. (Using Intervar db BS2_hom_het)

### BS3
Will be assigned if clinvar classifies the variant as Benign or likely benign and the level of evidence is either ‘practice guideline’ or ‘reviewed by expert panel’

### BP1
Will be assigned if the variant is Missense (nonsynonymous SNV) in a gene where missense variants represents at most 10 percent of all known pathogenic variants (using PP2_BP1 database).

### BP3
Will be assigned if the variant is an in-frame deletion/insertion in a repeat region of the gene. (Using the repeat_dict database).

### BP4
Will be assigned if the variant is predicted to be benign using various in-silico prediction tools (sift, lrt, mutationtaster, mutation assessor, fathmm, provean, meta svm, meta lr, mcap, mkl, genocanyon, gerp). Each prediction tool will add +1 or -1. A score of 0 or under will assign BP4 to a variant.

### BP6
Will be assigned the variant is classified as Benign or likely benign by clinvar but the evidence is limited.

### BP7
Will be assigned if a variant if synonymous and no splicing impact is predicted by dbscSNV (score under 0.6)
 
