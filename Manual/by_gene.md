# --by_gene report

By using the tag ```--by_gene``` with the ```analyse``` or ```sort``` option, you will generate a report that will group variants by gene.

Each gene will be assigned a "mutational burden" score that is useful to assess how much this gene is affected by variants in a cohort (if the input file is a multi-sample vcf).   
  
The gene burden score is the sum of the probability of pathogenicity (__0.8 and above__ otherwise excluded) or variant __i__ multiplied by the number of individuals affected by variant __i__.

A warning will be associated if the genes is a FLAG genes (frequently mutated in exome studies), if the number of sample affected seems to be to high (half of the variants in a gene with more than half of the individuals affected) or if the gene is especially long (more than 250,000bp, which are expected to harbour more variants). 

This report is useful to detect genes frequently mutated in a cohort sharing a same phenotype (but where several different variants are present).

### Example

| ERCC2 | 50.726   |          |     |     |        |                    |                  |                       |                                                                                                                                                                                |
|-------|----------|----------|-----|-----|--------|--------------------|------------------|-----------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Chr   | Start    | End      | Ref | Alt | Allele | ExonicFunc.refGene | Probability_Path | Prediction_ACMG_tapes | Samples                                                                                                                                                                        |
| 19    | 45867550 | 45867550 | C   | T   | T      | missense_variant   | 0.9993           | Pathogenic            | Sample_12, Sample_15, Sample_18, Sample_21, Sample_24, Sample_27, Sample_28, Sample_34, Sample_35, Sample_37, Sample_38, Sample_40, Sample_41, Sample_42, Sample_45, Sample_47 |
| 19    | 45856371 | 45856371 | G   | A   | A      | missense_variant   | 0.9993           | Pathogenic            | Sample_12, Sample_15, Sample_18, Sample_21, Sample_24, Sample_27, Sample_28, Sample_34, Sample_35, Sample_37, Sample_38, Sample_40, Sample_41, Sample_42, Sample_45, Sample_47 |
| 19    | 45860626 | 45860626 | G   | C   | C      | missense_variant   | 0.9941           | Pathogenic            | Sample_12                                                                                                                                                                      |
| 19    | 45867532 | 45867532 | C   | T   | T      | missense_variant   | 0.9941           | Pathogenic            | Sample_16, Sample_36                                                                                                                                                           |
| 19    | 45855805 | 45855807 | CTG | CG  | G      | frameshift_variant | 0.9941           | Likely Pathogenic     | Sample_27                                                                                                                                                                      |
| 19    | 45855574 | 45855574 | G   | A   | A      | missense_variant   | 0.9878           | Likely Pathogenic     | Sample_16, Sample_36                                                                                                                                                           |
| 19    | 45860760 | 45860760 | C   | T   | T      | missense_variant   | 0.9878           | Likely Pathogenic     | Sample_3, Sample_23, Sample_31                                                                                                                                                 |
| 19    | 45856059 | 45856059 | C   | G   | G      | missense_variant   | 0.8999           | Likely Pathogenic     |                                                                                                                                                                                |
| 19    | 45860760 | 45860760 | C   | T   | T      | missense_variant   | 0.8121           | VUS                   | Sample_21                                                                                                                                                                      |
| 19    | 45858047 | 45858047 | C   | T   | T      | missense_variant   | 0.8121           | VUS                   | Sample_26                                                                                                                                                                      |
| 19    | 45868168 | 45868168 | G   | C   | C      | missense_variant   | 0.8121           | VUS                   | Sample_3, Sample_21, Sample_23, Sample_27, Sample_31, Sample_32, Sample_45, Sample_46                                                                                          |
| 19    | 45867139 | 45867139 | T   | C   | C      | missense_variant   | 0.8121           | VUS                   | Sample_21                                                                                                                                                                      |