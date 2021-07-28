
# TAPES ![BuildUpdates](https://www.travis-ci.com/a-xavier/tapes.svg?branch=master)
 
**TAPES** : a Tool for Assessment and Prioritisation in Exome Studies  

Please refer to the [WIKI](https://github.com/a-xavier/tapes/wiki) for more information on all TAPES functions.  
See [HERE for installation and dependencies resolutions](https://github.com/a-xavier/tapes/wiki/Installation-and-Dependencies).

## Warning  
The command ```perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avdblist ./``` does not seem to download the list of databases available for ANNOVAR anymore (this might be temporary).
If this happens, please download:  
[https://github.com/a-xavier/tapes/files/5507866/hg19_avdblist.txt](https://github.com/a-xavier/tapes/files/5507866/hg19_avdblist.txt)  
and  
[https://github.com/a-xavier/tapes/files/5507867/hg38_avdblist.txt](https://github.com/a-xavier/tapes/files/5507867/hg38_avdblist.txt)  
and place them in your ANNOVAR folder (where ```annotate_variation.pl``` is located).


## TAPES purpose
TAPES is an open-source tool written in __python3__ designed to __predict pathogenicity__ of variants in __exome studies__.
It is especially designed for small scalle exome studies (from dozens to a few hundreds) sharing a common trait/phenotype.
TAPES can also handle __trio data__ with healthy parents and __calculate variant enrichment__ compared to the general population wwithout the use of controls. 

## CONTRIBUTING
Have an awesome idea, correction, feedback for TAPES? Want to contribute, please add an entry in the [Issues](https://github.com/a-xavier/tapes/issues) or a [Pull Request](https://github.com/a-xavier/tapes/pulls).

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
- Output VCF files
- Analysing InterVar output files
- Using real control samples for OR calculation/PS4
- Docker image creation
- Polygenic Risk Score calculation without controls
   
### Reference
   
If you use TAPES, please cite the [reference paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007453):  
   
Xavier A, Scott RJ, Talseth-Palmer BA (2019) TAPES: A tool for assessment and prioritisation in exome studies. PLOS Computational Biology 15(10): e1007453. https://doi.org/10.1371/journal.pcbi.1007453
  
  
