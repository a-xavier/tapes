# ```sort```  

sort is the main function of TAPES. It will prioritise variants and generate reports.

## INPUT  
```-i``` or ```--input```
The input file will always be an ANNOVAR annotated file. It can be either:
 - A vcf file (multi or single sample)
 - A txt (tab separated) or csv (comma separated) file (ANNOVAR does not keep genotyping data in those format but you can use TAPES to annotate you original vcf files to the format of your liking)

## OUTPUT 
```-o``` or ```--output```
The output file can be:
- a txt (tab separated) or csv (comma separated) accompanied by an xlsx file
- a folder containing either multiple txt (tab separated) or csv (comma separated) files.

## EXAMPLES 
### Simple examples

```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output.txt ``` for txt+xlsx output  
```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output.csv ``` for csv+xlsx output  
```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ ``` for folder csv output  
```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --tab ``` for folder txt output  
### Example of full command
```python3 tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output/ --tab --acmg --by_sample --by_gene --enrichr --disease "autosomal dominant" --kegg "Pathways in cancer" --list "ATM APC MUTYH AKT1 CASP8 DMD MSH3 MSH6 MLH1 MSH2 POLE POLD POLQ KRAS" -a hg19 --BIG --trio ./trio.txt  ```

## Reports 

The different reports generated (if using the FOLDER mode) will have different suffixes appended to them.
For example:  

- report_**by_sample**.txt for the ```--by_sample``` option
- report_**GO_Biological_Process_2018**.txt for the ```--enrichr``` option using the GO_Biological_Process_2018 library

## OPTIONS  

* ```-a``` ```--assembly``` |_string_| The assembly used. Either hg19 or hg38 - default is hg19  
 * ```--acmg```|_flag_| Will check for all necessary annotations before proceeding. Will exit the program if alle the necessary annotations were not found. Do not use this unless your annotation file is fully "TAPES compliant"  
* ```--BIG```|_0 to 1 float_|If your file is particularly big, remove variants that are not likely to be pathogenic from the output.Take the percentage of chance to be pathogenic as an argument - default = 0.35  
*  ```--by_gene```|_flag_| Will generate a report of the mutational burden for each gene with associated samples for each variant.
* ```--by_sample```|_flag_|Will generate a report grouping all pathogenic variants per sample  
* ```--cutoff```|_0 to 1 float_| Select the cutoff for "rare" variant. Used in ACMG BS1 - default = 0.005  
* ```-d``` ```--disease``` |_string_| Will generate a report similar to the main report including only the variants involved in a particular disease. Based on the Disease column from ANNOVAR - default = "cancer"  
* ```-e``` ```--enrichr```|_string_| Implementation of the EnrichR API. Analyse top mutations (>0.85 in Pathology prediction) with EnrichR API to determine the pathways disrupted by pathogenic variants - default = GO_Biological_Process_2018 -
See [the EnrichR librairies page](https://github.com/a-xavier/tapes/wiki/EnrichR-Libraries) for the full list of available EnrichR Libraries  
* ```--kegg``` |_string_|Similar to list but if you do not know all the genes involved in a given pathway. Will generate a report similar to the main report including only variants found in the genes related to a given pathway. See [the KEGG pathways page](https://github.com/a-xavier/tapes/wiki/KEGG-pathways) for the full list of available pathways  


* ```--list``` |_string_ or _path_| Will generate a report similar to the main report including only variants found in genes list. Using gene symbols. Can use either a string containing gene symbols separated by a space and in quotes(eg. "MLH1 APC MLH2 BRCA2 PTEN") or a txt file with one gene symbol per line  
* ```-t``` ```--threads``` |int| Number of threads to use for sorting | Still not optimised and require a large amount of ram ( around 2Gb per thread + 2*size of the input file) 
* ```--tab``` |_flag_| If the output mode is folder, outputs txt (tab separated) files instead of csv (comma separated)

* ```--trio```|_path_| Path to a tab separated txt file. Use if you have trio data with unnaffected parents to detect de-novo mutations. Necessary for ACMG PS2 criteria. [See the trio wiki page for more information](https://github.com/a-xavier/tapes/wiki/Trio).

* ```--pp2_percent```|_0 to 100 int_|Threshold for PP2, considering PP2 positive if variant is missense in a gene where more than the threshold are pathogenic missense - default = 80  
* ```--pp2_min```|_int_| Number of minimum pathogenic missense variants in gene to consider PP2  
* ```--bp1_percent``` |_0 to 100 int_|Threshold for BP1, considering BP1 positive if variant is missense in a gene where less than the threshold are pathogenic missense - default = 15



