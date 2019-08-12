# TAPES [![Build Status](https://travis-ci.org/a-xavier/tapes.svg?branch=testing)](https://travis-ci.org/a-xavier/tapes)
**TAPES** : a Tool for Assessment and Prioritisation in Exome Studies

Please refer to the [Instruction Manual](https://github.com/a-xavier/tapes/blob/master/TAPES_Manual.pdf) or the [Wiki](https://github.com/a-xavier/tapes/wiki) for more detailed explanations on every option.   
  
TAPES can:  
- use [ANNOVAR](annovar.openbioinformatics.org) (or [VEP](https://ensembl.org/info/docs/tools/vep/index.html)) annotated vcf/txt/csv files to prioritize variant, predict their pathogenicity (using both the ACMG Criteria and a direct pathogenicity probability ranked from 0 to 1).
- generate useful reports in multi-sample vcf/cohorty studies
- be used as an simplified interface for ANNOVAR.

Please note that TAPES does not require an ANNOVAR installation to sort already annotated files.  
TAPES is mostly designed for mendelian diseases and exome studies of cohorts sharing a common phenotype but can prioritise any ANNOVAR annotated vcf file.

TAPES is written for python3. You can use pip (pip3) for dependencies.

### Download  

Either use the download button [from the github website](https://github.com/a-xavier/tapes/archive/master.zip) in a browser or use the command:  
```git clone https://github.com/a-xavier/tapes```

### Resolve dependencies 

run ```pip3 install -r requirements.txt --user ``` on Linux/macOS  

run ```pip3 install -r winrequirements.txt --user ``` on Windows
  
:warning:Depending on your distribution, installing **tkinter** might be necessary to generate graphs.  
It can be called ```python3-tk``` for debian/ubuntu/fedora or just ```tk``` on arch-based linux.

## Quick start : prioritisation  

See the [QUICK START](https://github.com/a-xavier/tapes/wiki/Quick-Start) section on the wiki.

Using ```python``` or ```python3``` depending on your default installation:

```python tapes.py sort -i /path/to/annotated_file.vcf -o /path/to/output.txt --reporting_options```

##### Example of sorting and reporting options

* ```--acmg``` to perform a check of all necessary annotations | _flag_  
* ```--by_sample``` to generate a report of the 5 (or less) most pathogenic variants for each sample in the cohort | _flag_  
* ```--by_gene``` to generate a report of the gene mutation burden for each gene | _flag_  
* ```--enrichr [library]``` to generate a pathway analysis of the genes affected by pathogenic mutations | [default = GO_Biological_Process_2018]  

See [Instruction Manual](https://github.com/a-xavier/tapes/blob/master/TAPES_Manual.pdf) or the [Wiki](https://github.com/a-xavier/tapes/wiki) for the full range of options, functionalities and [required annotations](https://github.com/a-xavier/tapes/wiki/Necessary-Annotations).   

#### Toy Dataset  
You can test TAPES using the simultated datased located in the folder _toy_dataset_.  
Open a terminal and ```cd``` to the TAPES folder.   

Then try pasting this command in the terminal and press enter:    

```python tapes.py sort -i ./toy_dataset/toy_vep_multi.vcf -o ./report/ --tab --by_gene --by_sample --enrichr --disease "autosomal dominant" --kegg "Pathways in cancer"```  

If this does not work, try: 

```python3 tapes.py sort -i ./toy_dataset/toy_vep_multi.vcf -o ./report/ --tab --by_gene --by_sample --enrichr --disease "autosomal dominant" --kegg "Pathways in cancer"```     

This should create a folder called Toy_dataset containing various reports.


:warning: __Warning on VEP annotated vcf :__ VEP annotated vcf might give different results compared to ANNOVAR annotated vcf, due to the difference in interpretation of consequences. ANNOVAR is still the preferred annotation tool to use with TAPES.

___Coming Soon___...  
- ~~Support for VEP annotated vcf and BS2 criteria assignment optimisation~~ Done
- ~~Add ```--by_gene``` report to files with only one sample.~~ Done
- ~~Better enrichement calculation for PS4 criteria assignment~~ Done  
- Output VCF files
- Analysing InterVar output files
- Using real control samples for OR calculation/PS4

__Odds Ratio (OR) calculation (for PS4 criteria) was recently changed__ to be closer to the reality. In short the old extrapolation from MAF is first calculated and then another OR calculation is made using a smaller frequency in the control population. Then a mean between the 2 results is calculated. This represents (around) 30% less difference between the extrapolated OR and the reality. Keep in mind that TAPES OR calculation will always be more stringent than the normal calculation to avoid excessive false positive, meaning OR will (nearly) always be lower (only enrichment is tested).
(see graph and table below)

| Affected in general population | Unaffected in general population | Frequency | Normal OR calculation | Old TAPES extrapolation | New TAPES extrapolation |
|---|---|---|---|---|---|
|1	|9999|	1.00E-04|	256.38|	128.18|	128.165|
|2	|9998|	2.00E-04|	128.18|	64.08|	96.105|
|3	|9997|	3.00E-04|	85.44|	42.71|	85.405|
|4	|9996|	4.00E-04|	64.08|	32.03|	48.03|
|5	|9995|	5.00E-04|	51.26|	25.62|	44.815|
|6	|9994|	6.00E-04|	42.71|	21.34|	42.67|
|7	|9993|	7.00E-04|	36.6|	18.29|	30.47|
|8	|9992|	8.01E-04|	32.03|	16|	29.32|
|9	|9991|	9.01E-04|	28.46|	14.22|	28.425|
|99	|9901|	1.00E-02|	2.56|	1.26|	2.5|
|100|	9900|	1.01E-02|	2.54|	1.26| 2.465|
|101|	9899|	1.02E-02|	2.51|	1.26|	2.465|
|102|	9898|	1.03E-02|	2.49|	1.26|	2.465|
|103|	9897|	1.04E-02|	2.46|	1.26|	2.465|
|104|	9896|	1.05E-02|	2.44|	1.26|	2.465|
|105|	9895|	1.06E-02|	2.42|	0.62|	1.215|
|106|	9894|	1.07E-02|	2.39|	0.62|	1.215|
|107|	9893|	1.08E-02|	2.37|	0.62|	1.215|

**Assuming a frequency of 0.025 in the case cohort**

![alt text](https://raw.githubusercontent.com/a-xavier/tapes/master/Example_Output/New%20PS4%20calc.png "OR calculation")


