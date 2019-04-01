# TAPES
**TAPES** : a Tool for Assessment and Prioritisation in Exome Studies

Please refer to the [Instruction Manual](https://github.com/a-xavier/tapes/blob/master/TAPES_Manual.pdf) for more detailed explanations on every option.   
  
TAPES can:  
- uses [ANNOVAR](annovar.openbioinformatics.org) or VEP annotated vcf/txt/csv files to prioritize variant, predict their pathogenicity (using both the ACMG Criteria and a direct pathogenicity probability ranked from 0 to 1).
- generate useful reports in multi-sample vcf/cohorty studies
- be used as an simplified interface for ANNOVAR  

TAPES is mostly designed for mendelian diseases and exome studies of cohorts sharing a common phenotype but can prioritise any ANNOVAR annotated vcf file.

TAPES is written for python3. You can use pip (pip3) for dependencies.
### Resolve dependencies 

run ```pip3 install -r requirements.txt --user ```

## Quick start : prioritisation

Using ```python``` or ```python3``` depending on your default installation:

```python tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output.txt --reporting_options```

##### Example of sorting and reporting options

* ```--acmg``` to perform a check of all necessary annotations | _flag_  
* ```--by_sample``` to generate a report of the 5 (or less) most pathogenic variants for each sample in the cohort | _flag_  
* ```--by_gene``` to generate a report of the gene mutation burden for each gene | _flag_  
* ```--enrichr [library]``` to generate a pathway analysis of the genes affected by pathogenic mutations | [default = GO_Biological_Process_2018]  

See [Instruction Manual](https://github.com/a-xavier/tapes/blob/master/TAPES_Manual.pdf) or the [Wiki](https://github.com/a-xavier/tapes/wiki) for the full range of options and functionalities.   

:warning: __Warning on VEP annotated vcf :__ VEP annotated vcf might give different results compared to ANNOVAR annotated vcf, due to the difference in interpretation of consequences. ANNOVAR is still the preferred annotation tool to use with TAPES.

___Coming Soon___...  
- ~~Support for VEP annotated vcf and BS2 criteria assignment optimisation~~ Done
- Add ```--by_gene``` report to files with only one sample. 

