# This toy dataset is composed of made-up vcf files with the data from: 

Zhang J, Walsh MF, Wu G, Edmonson MN, Gruber TA, Easton J, et al. Germline Mutations in Predisposition Genes in Pediatric Cancer. N Engl J Med. 2015;373(24):2336-46.  

This dataset was used for benchmarking TAPES.

 - toy.vcf is a simple unannotated vcf file
 - toy_annovar.vcf is the same toy.vcf annotated with ANNOVAR
 - toy_vep.vcf is the same toy.vcf annotated with VEP
 
 - toy_multi.vcf is a simple unannotated multi-sample vcf file
 - toy_annovar_multi.vcf is the same multi-sampletoy.vcf annotated with ANNOVAR
 - toy_vep_multi.vcf is the same multi-sample toy.vcf annotated with VEP  

# TRY
-- Unannotated VCF file
  - I just want to use the necessary ANNOVAR annotations for TAPES sorting
    - Download [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
    - ```python tapes.py db -s -A /path/to/annovar/```
    - ```python tapes.py db -b --acmg --assembly hg19```
    - ```python tapes.py annotate -i ./toy_dataset/toy.vcf -o ./toy_dataset/anno.vcf --acmg –a hg19```
  
-- VEP or ANNOVAR annotated file
  - I want to sort and prioritise my VCF file
    - ```python tapes.py sort -i ./toy_dataset/toy_annovar_multi.vcf  –o ./toy_report/ --tab --by_gene --by_sample --enrichr --disease      "autosomal dominant" --kegg "Pathways in cancer" --list "MLH1 MSH2 MSH6" --BIG```

 
