# ```annotate```
annotate will use ANNOVAR to annotate a vcf file using a simplified command line.

If you have not already, you should read the entry about [database management](https://github.com/a-xavier/tapes/wiki/Database-Management) first.  
TAPES will accept **vcf**, **bcf**, **gzipped bcf** and **vcf** as well as **bgzipped vcf** as input for annotation.
All formats will be **converted to decomposed vcf** prior to annotation.

### Input file 
```-i``` ```--input```
TAPES will accept vcf, bcf, gzipped bcf and bgzipped vcf as input files for the ```annotate``` option
### Output file
```-o``` ```--output```
TAPES allow to output either annotated csv (comma separated), txt (tab separated) or vcf files. 

**WARNING** If your vcf is multi-sample and your output is csv or txt, two files will be created. One file without the sample genotyping and one file with the sample genotyping data, which will have the suffix: "_with_samples"

## Simplified annotation  
To annotate a vcf files using just the necessary annotation for TAPES, then run:  
```python3 tapes.py annotate -i /path/to/variant_file.bcf -o /path/to/annotated_file.vcf --amcg```

The ```--acmg``` tag will only use the [necessary annotations](https://github.com/a-xavier/tapes/wiki/Necessary-Annotations) for TAPES.

## Full annotations and db_vcf.json
If you want to annotate your files with other annotations, open the file **db_vcf.json** with any editor. It will show every databases detected in your annovar folder. Replace "**NO**" by "**YES**" for the databases you wish to use for annotation.   
Then run:  
```python3 tapes.py annotate -i /path/to/variant_file.vcf -o /path/to/annotated_file.vcf --ref_anno knownGene ```

## OPTIONS
* ```-a``` ```--assembly```|_string_| The assembly used. Either hg19 or hg38 - default = hg38
* ```--acmg```|_flag_| Bypassed the **db_vcf.json** file and only annotate the variant file with the necessary annotations for TAPES sorting.
* ```--ref_anno```|_string_| The reference annotation system used. Either 'refGene' for RefSeq, 'ensGene' for ENSEMBL or 'knownGene' for UCSC annotation - default = refGene

