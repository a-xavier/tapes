# TAPES
**TAPES** : a Tool for Assessment and Prioritisation in Exome Studies

Please refer to the [Instruction Manual](https://github.com/a-xavier/tapes/blob/master/TAPES_Manual.pdf) for more detailed explanations on every option.   
  
TAPES uses [ANNOVAR](annovar.openbioinformatics.org) annotated vcf/txt/csv files to prioritize variant, predict their pathogenicity and generates useful reports.  
TAPES is mostly designed for mendelian diseases and exome studies of cohorts sharing a common phenotype but can prioritise any vcf file.

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

___Coming Soon___...Support for VEP annotated vcf and BS2 criteria assignment optimisation
