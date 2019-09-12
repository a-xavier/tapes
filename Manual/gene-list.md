# --list 

#### Used to filter the main report based on user-supplied genes.

usage:

```python3 tapes.py sort -i ./input.vcf -o ./output_folder/ --list "MLH1 PMS2 MSH6 MSH2 APC"```  
  
```python3 tapes.py sort -i ./input.vcf -o ./output_folder/ --list ./gene_list.txt```


By using the option ```--list``` with the ```sort``` or ```analyse``` option, a report containing only variants located in the user-supplied list of genes.  

```--list``` requires gene SYMBOLS, either __in quotes separated by a space__ or as a __.txt file with one gene symbol per lines__.

The format of this report will be identical to the main report.