### ANNOVAR Necessary annotations

TAPES need a number of annotations to give the best possible estimation of a variant pathogenicity.
You can use the ```--acmg``` flag when annotating with TAPES to directly use them.

The necessary/recommended databases are the following:   

|Annotation Required	|ACMG Classification Criteria|
|---|---|
|Gene reference annotation (Refseq, ENSEMBL or UCSC)|	All|
|dbscSNV	|PVS1|
|gnomad_genome and either gnomad_exome or exac|	PS1|
|Clinvar (20151201 or above)	|PS3 / PP5 / BP6|
|DbNSFP (30a or above)|	PP3 /BP4|
|dbsnp or avsnp	|PS4|  

If one or more these annotations are missing. Do not use the ```--acmg``` flag while sorting.  
TAPES will prioritise the variant using as many annotation as possible but will obviously be less powerfull.  

### VEP Necessary annotations  

TAPES processing of VEP annotated vcf relies heavily on dbNSFP annotations. 
RefSeq reference annotations are recommended for VEP annotated vcf processing.

__Necessary dbNSFP annotations__ : 
- Interpro (for domain annotation)
- for in-silico prediction:
  - SIFT
  - LRT
  - MutationTaster 
  - MutationAssessor
  - FATHMM
  - fathmm-MKL
  - PROVEAN 
  - MetaSVM and MetaLR
  - M-CAP
  - GenoCanyon
  - GERP++
  - Polyphen-2
- clinvar
- gnomad_exome __or__ ExAC
- gnomad_genome
- dbscSNV


