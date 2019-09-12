## Warnings  

Here are a few things to keep in mind when using TAPES

### ACMG classification and Pathogenicity prediction  
The ACMG criteria and classification (from Benign to Pathogenic or class 1 to 5) reflect the **PROBABILITY** of a variant to be pathogenic. A variant of class 4 (Likely Pathogenic) is **NOT MORE** pathogenic than a variant of class 5 (Pathogenic). Instead, it is just **MORE LIKELY** to be pathogenic.  

The same goes with the ACMG criteria modeling from [Tavtigian et al](https://www.ncbi.nlm.nih.gov/pubmed/29300386) that TAPES uses.  

This has implications in the ```--by_gene``` score used in TAPES. This score is the sum of all probabilities multiplied by the number of sample affected. It reflect the probability of the genes to be affected by a pathogenic variant in the cohort. 

