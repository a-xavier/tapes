# --kegg  

#### Used to filter the main report based on functional kegg ([Kyoto Encyclopedia of Genes and Genome](https://www.genome.jp/kegg/pathway.html)) pathways.

usage:

```python3 tapes.py sort -i ./input.vcf -o ./output_folder/ --kegg "fc epsilon ri signaling pathway"```

By using the option ```--kegg``` with the ```sort``` or ```analyse``` option, a report containing only variants located in the genes contained in the designated pathway will be created. 

A list of [KEGG patwhays can be found here](https://github.com/a-xavier/tapes/wiki/KEGG-pathways).  
The format of this report will be identical to the main report.