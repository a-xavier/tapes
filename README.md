
# TAPES  [![Build Status](https://travis-ci.org/a-xavier/tapes.svg?branch=master)](https://travis-ci.org/a-xavier/tapes)
**TAPES** : a Tool for Assessment and Prioritisation in Exome Studies  

Please refer to the [WIKI](https://github.com/a-xavier/tapes/wiki) for more information on all TAPES functions.  
See [HERE for installation and dependencies resolutions](https://github.com/a-xavier/tapes/wiki/Installation-and-Dependencies).

## TAPES purpose
TAPES is a tools written in __python3__ designed to __predict pathogenicity__ of variants in __exome studies__.
It is especially designed for small scalle exome studies (from dozens to a few hundreds) sharing a common trait/phenotype.
TAPES can also handle __trio data__ with healthy parents and __calculate variant enrichment__ compared to the general population wwithout the use of controls. 

## TAPES WORKFLOW 

__Unnanotated VCF__ --```Annotate```--> __Annotated VCF__ --```Sort```--> __Reports__  
TAPES recommanded input is a multi-sample vcf or bcf file.

#### Toy dataset

[Follow the instructions here to try TAPES' workflow using a toy dataset.](https://github.com/a-xavier/tapes/wiki/Workflow)

![Workflow Diagram](https://raw.githubusercontent.com/a-xavier/tapes/testing/acmg_db/workflow.png)
## TAPES Functionalities
### ANNOVAR interface (Linux only)
TAPES can be used as an [ANNOVAR](annovar.openbioinformatics.org) wrapper for easy VCF/gzipped VCF/BCF annotation (this requires to download ANNOVAR). 

```python3 tapes.py db -s -A /path/to/annovar/```  to scan the annovar path for existing databases  
```python3 tapes.py db -b --acmg -a hg19```    to download all necessary databases  
```python3 tapes.py annotate -i /path/to/variant_file.bcf -o /path/to/annotated_file.vcf --amcg```  to annotate the bcf or vcf file of interest  

### Prioritise and Sort Annotated VCF
TAPES can predict the pathogenicity of variants (using [ANNOVAR](https://github.com/a-xavier/tapes/wiki/Necessary-Annotations#annovar-necessary-annotations) or [VEP](https://github.com/a-xavier/tapes/wiki/Necessary-Annotations#vep-necessary-annotations) annotated VCF):
  - by assigning the [ACMG/AMP Criteria](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)
  - provide a probability (from 0 to 1) of a variant to be pathogenic based on [Tavtigian et al model](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6336098/)  
  
  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --tab```  
  to sort and prioritise an annotated vcf.
  
  #### PS4 criteria/ Variant enrichment calculation
  TAPES can calculate the enrichment of variants in a population sharing the same phenotype without the need for controls, using public databases and extrapolating the number of individuals affected/non-affected by a variant using Minor Allele Frequencies from public databases.  
  [Check the details HERE](https://github.com/a-xavier/tapes/wiki/PS4-calculation).
  ### Reporting options
  TAPES provides an [array of useful reports]() on top of the main report (which will prioritise variants based on predicted pathogenicity)
  #### Filtering options   
  - ```--list```  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --by_gene```  
    to filter the main report with a user-supplied gene symbol list.
  - ```--disease```  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --by_gene```  
    to filter the main report based on a disease characteristics (eg. 'autosomal dominant' or 'colorectal')
  - ```--kegg```  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --by_gene```
     to filter the main report based on a kegg pathway (eg. "Pathways in cancer" or "base excision repair")
  
  #### By-gene report
  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --by_gene```    
  to create a by gene report with an associated gene mutational burden score.
  
  #### Pathway analysis
  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --enrichr GO_Biological_Process_2018```   
  to perfome pathway analysis with genes containing pathogenic variants.
  
  #### By-sample report
  ```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --by_sample```     
  to create a table of the most pathogenic variants per individual.
  
  ### Coming soon
  - ~~Support for VEP annotated vcf and BS2 criteria assignment optimisation~~ Done
- ~~Add ```--by_gene``` report to files with only one sample.~~ Done
- ~~Better enrichement calculation for PS4 criteria assignment~~ Done  
- ~~Polygenic Risk Score calculation without controls~~ Done 
- Output VCF files
- Analysing InterVar output files
- Using real control samples for OR calculation/PS4
- Docker image creation

### New Polygenic Risk Score (PRS) Module  

TAPES can now use a multi-sample vcf file to calculate a Polygenic Risk Score for a specific disease/trait.  
PRS can detect an elevated risk for a specific disease using common variants that add a very small but significant risk.

For now only cumulative PRS is calculated as follow

PRS=  sum { βi * SNPi }

Where βi is the beta-value associated with the SNPi and SNPi is the genotype at this locus (1 for heterozygous and 2 for homozygous)

PRS can help determine if the studied cohort is more at risk of a certain disease (vs n "healthy" controls from 1000genomephase3).

Samples are compared to public genotyped individuals from [ENSEMBL REST API](https://rest.ensembl.org/documentation/info/variation_id). Odds Ratios and Beta values for each disease/trait was obtained through [GWAS Catalog](https://www.ebi.ac.uk/gwas/). 

__CAVEATS__
- p-values from GWAS catalog are all under 9x10-6.  
- Linkage Desiquilibrium is __not__ adjusted so proceed with caution with your results (some SNPs might not be independants).  
- SNPs that are heavily linked to ethnicity are excluded (using a list from [Huang et al. 2015](https://www.ncbi.nlm.nih.gov/pubmed/26690364)  
- This function will work best with whole-exomes rather than targeted.
- If the is a very small number of SNP (n<20) overlapping between your cohort and the SNPs in GWAS catalog, the results will most likely be very inacurrate.
- As always, the more samples you use, the more precise this function will be.

__USAGE__ 

With a multi-sample annotated file, try:  
```python3 tapes.py sort -i ./annotated_multi_sample.vcf -o ./output_folder/ --prs "Crohn's disease" --prs_num 200```

__OUTPUT__

It will output:  

- A table with the genotype of your samples along with the control samples for the specified SNPs as well as informations on the SNPs (beta value, gene, frequency, etc). In addition, means from cases and controls are compared using a non-parametric test.

- A histogram showing ranked (from lowest to biggest) PRS. Blue bars are healthy public controls and red bars are the sample from the cohort studied. 


![PRS_picture](https://raw.githubusercontent.com/a-xavier/tapes/testing/acmg_db/PRS_CD.png)

  
  
