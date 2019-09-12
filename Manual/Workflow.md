# TAPES Workflow using the provided toy dataset
Here you can find an example for workflow starting from a __*multi-sample VCF*__ file and a fresh [ANNOVAR download](http://www.openbioinformatics.org/annovar/annovar_download_form.php)  

## Using ANNOVAR wrapping

First, ```cd``` to the main TAPES folder.  

First use:  
```python3 tapes.py db -s -A /path/to/annovar/```  
To determine the ANNOVAR folder.  

Then:  
```python3 tapes.py db -b --acmg```  
To download necessary databases (this might take some time).  

Then:  
```python3 tapes.py annotate -i ./toy_dataset/toy_multi.vcf -o ./annotated_multi.vcf --acmg```   

To annotate the original vcf.

Then either: 
* ```python3 tapes.py sort -i ./annotated_multi.vcf -o ./toy_report/ --tab```    

To prioritise the annotated vcf and generate the main output.  
  
**or**  
   
* ```python3 tapes.py sort -i ./annotated_multi.vcf -o ./toy_reports/ --tab --by_gene --by_sample --disease cancer --enrichr --list "TP53 BRCA1 BRCA2 PTEN KRAS" --trio ./toy_dataset/trio.txt   --kegg "pathways in cancer"```

To prioritise the annotated vcf, generate the main output and other few secondary reports

## Without ANNOVAR WRAPPING 

Try 

* ```python3 tapes.py sort -i ./toy_dataset/toy_annovar_multi.vcf -o ./toy_report/ --tab```    

To prioritise the annotated vcf and generate the main output.  
  
**or**  
   
* ```python3 tapes.py sort -i ./toy_dataset/toy_vep_multi.vcf -o ./toy_reports/ --tab --by_gene --by_sample --disease cancer --enrichr --list "TP53 BRCA1 BRCA2 PTEN KRAS" --trio ./toy_dataset/trio.txt   --kegg "pathways in cancer"```

To prioritise the annotated vcf, generate the main output and other few secondary reports

