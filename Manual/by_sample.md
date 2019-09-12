# -- by_sample report

By using the tag ```--by_sample``` with the ```analyse``` or ```sort``` option, you will generate a report that will show for each sample the 5 most probably pathogenic variants.

This obviously requires multi-sample vcf.

## Example 
| Sample_9  |          |          |     |     |        |                    |              |                  |                       |
|-----------|----------|----------|-----|-----|--------|--------------------|--------------|------------------|-----------------------|
| Chr       | Start    | End      | Ref | Alt | Allele | ExonicFunc.refGene | Gene.refGene | Probability_Path | Prediction_ACMG_tapes |
| 17        | 7574017  | 7574017  | C   | T   | T      | missense_variant   | TP53         | 0.9999           | Pathogenic            |
| 17        | 7574017  | 7574017  | C   | T   | T      | missense_variant   | TP53         | 0.9999           | Pathogenic            |
| 13        | 32953932 | 32953932 | T   | A   | A      | stop_gained        | BRCA2        | 0.9999           | Pathogenic            |
| 17        | 7574017  | 7574017  | C   | T   | T      | missense_variant   | TP53         | 0.9999           | Pathogenic            |
| 14        | 45606305 | 45606305 | G   | A   | A      | stop_gained        | FANCM        | 0.9998           | Pathogenic            |
|           |          |          |     |     |        |                    |              |                  |                       |
| Sample_12 |          |          |     |     |        |                    |              |                  |                       |
| Chr       | Start    | End      | Ref | Alt | Allele | ExonicFunc.refGene | Gene.refGene | Probability_Path | Prediction_ACMG_tapes |
| 13        | 32968863 | 32968863 | C   | G   | G      | stop_gained        | BRCA2        | 1                | Pathogenic            |
| 17        | 59793412 | 59793412 | G   | A   | A      | stop_gained        | BRIP1        | 1                | Pathogenic            |
| 17        | 7574017  | 7574017  | C   | T   | T      | missense_variant   | TP53         | 0.9999           | Pathogenic            |
| 17        | 7574017  | 7574017  | C   | T   | T      | missense_variant   | TP53         | 0.9999           | Pathogenic            |
| 17        | 7574017  | 7574017  | C   | T   | T      | missense_variant   | TP53         | 0.9999           | Pathogenic            |