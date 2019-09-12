# Quick Start (Choose a starting point)

## Unannotated VCF file

### I just want to use the necessary ANNOVAR annotations for TAPES sorting
 - Download [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
 - ```python3 tapes.py db -s -A /path/to/annovar/```  
 - ```python3 tapes.py db -b --acmg --assembly hg19```
 - ```python3 tapes.py annotate -i /path/to/file.vcf -o /path/to/output.vcf --acmg -a hg19```
### I want to customise the ANNOVAR annotations used
 - See detailed instructions on [Database Management](https://github.com/a-xavier/tapes/wiki/Database-Management#full-database-management) and [Custom annotations](https://github.com/a-xavier/tapes/wiki/Annotate-option#full-annotations-and-db_vcfjson)

## VEP or ANNOVAR annotated file
 ### I want to sort and prioritise my annotated VCF file
 - ```python3 tapes.py sort -i /to/annotated/file.vcf -o /to/output/folder/ --tab```
 ### I want to re-analyse TAPES main report
 - ```python3 tapes.py analyse -i /to/main/report.txt -o /to/output/report.txt --single_option```

# Using the toy dataset  

try:  
```python3 tapes.py sort -i ./toy_dataset/toy_annovar_multi.vcf -o ./Reports/ --tab --by_gene --by_sample --enrichr --list "MLH1 MSH6 MSH2" --disease "autosomal dominant" --kegg "pathways in cancer" ```
 