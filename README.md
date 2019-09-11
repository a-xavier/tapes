
# TAPES  [![Build Status](https://travis-ci.org/a-xavier/tapes.svg?branch=master)](https://travis-ci.org/a-xavier/tapes)
**TAPES** : a Tool for Assessment and Prioritisation in Exome Studies  

Please refer to the [WIKI](https://github.com/a-xavier/tapes/wiki) for more information on all TAPES functions.

## TAPES purpose
TAPES is a tools designed to predict pathogenicity of variants in exome studies.
It is especially designed for small scalle exome studies (from dozens to a few hundreds) sharing a common trait/phenotype.
TAPES can also handle trio data with healthy parents and calculate variant enrichment compared to the general population wwithout the use of controls.

## TAPES WORKFLOW 
![Workflow Diagram](https://raw.githubusercontent.com/a-xavier/tapes/testing/acmg_db/workflow.png)
## TAPES Functionalities
### ANNOVAR interface
TAPES can be used as an [ANNOVAR](annovar.openbioinformatics.org) wrapper for easy VCF/gzipped VCF/BCF annotation
### Prioritise and Sort Annotated VCF
TAPES can predict the pathogenicity of variants (using ANNOVAR(https://github.com/a-xavier/tapes/wiki/Necessary-Annotations#annovar-necessary-annotations) or [VEP](https://github.com/a-xavier/tapes/wiki/Necessary-Annotations#vep-necessary-annotations) annotated VCF):
  - by assigning the [ACMG/AMP Criteria](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)
  - provide a probability (from 0 to 1) of a variant to be pathogenic based on [Tavtigian et al model](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6336098/)
  #### PS4 criteria/ Variant enrichment calculation
  TAPES can calculate the enrichment of variants in a population sharing the same phenotype without the need for variants.
  Check the details [HERE].
  ### Reporting options
  
