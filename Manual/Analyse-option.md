# ```analyse``` or ```analyze```  
analyse allow users to generate secondary reports without having to generate the main report again, saving time.

## Input
```-i``` ```--input```  
Any kind of main report from TAPES either txt (tab separated) or csv (comma separated) file.

## Output
```-o``` ```--output```
Will output a txt (tab separated) or csv (comma separated) file containing the desired report

## Examples
```python tapes.py analyse -i /path/to/main_report.csv -o /path/to/new_report.csv --by_gene```  
```python tapes.py analyze -i /path/to/main_report.csv -o /path/to/new_report.txt --by_sample```  
```python tapes.py analyse -i /path/to/main_report.txt -o /path/to/new_report.txt --enrichr GO_Molecular_Function_2018```  


## Warning
analyse will only work one report at a time. For example:  
```python tapes.py analyse -i /path/to/main_report.csv -o /path/to/new_report.csv --by_gene --by_sample``` will __not__ work.  
## Available Report Options
See the entire details of options at the [sort function page](https://github.com/a-xavier/tapes/wiki/Sort-option).  
* ```--by_gene```
* ```--by_sample```
* ```--list``` and ```--kegg```
* ```--enrichr```
* ```--disease```

### Analysis strategies 
Options such as  ```--disease```, ```--list``` and ```--kegg``` will output a file with a similar format to the main report. The output file will only keep genes passing a filter.

This means that those can be used with other options such as ```--by_gene``` or ```--by_sample```.

For example:  
```python tapes.py analyse -i ./main_report.txt -o ./only_cancer_genes.txt --disease cancer```  
then  
```python tapes.py analyse -i ./only_cancer_genes.txt -o ./cancer_sample.txt --by_sample```  
Will give, for each sample, the 5 most pathogenic variant only in genes that are involved in cancer pathways.