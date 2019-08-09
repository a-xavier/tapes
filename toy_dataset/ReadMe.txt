This toy dataset is composed of made-up vcf files with the data from: 

Zhang J, Walsh MF, Wu G, Edmonson MN, Gruber TA, Easton J, et al. Germline Mutations in Predisposition Genes in Pediatric Cancer. N Engl J Med. 2015;373(24):2336-46.

 - toy.vcf is a simple unannotated vcf file
 - toy_annovar.vcf is the same toy.vcf annotated with ANNOVAR
 - toy_vep.vcf is the same toy.vcf annotated with VEP
 
 - toy_multi.vcf is a simple unannotated multi-sample vcf file
 - toy_annovar_multi.vcf is the same multi-sampletoy.vcf annotated with ANNOVAR
 - toy_vep_multi.vcf is the same multi-sample toy.vcf annotated with VEP  
 
 Try running from TAPES main folder:  
 
 python tapes.py sort -i ./toy_vep_multi.vcf -o ./Report/ --tab --by_gene --by_sample --enrichr --disease "autosomal dominant" --kegg "Pathways in cancer"
 
 or  
 
 python tapes.py sort -i ./toy_annovar_multi.csv -o ./Report.txt