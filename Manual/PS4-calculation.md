# PS4/Enrichment calculation

TAPES can calculates the PS4 ACMG criteria (enrichment of a variant in a diseased population vs general population) without the need of sequenced controls.

Calculating Odds Ratio using Fisher's exact test require integers while most of the databases only gives minor allele frequencies. 
## OLD METHOD
By assuming that most very rare variants are heterozygous, we extrapolate the number of individuals with the variant and the number of individuals without the variant in public databases (ExAC or gNomad).

If MAF in the controls population  
__MAF  = y^-x__   
then the number of individuals with the variant is   
__n = ceil(y)__  
and the number of individuals without the variant is 
__N = (x/2) - n__

__If MAF = 3.23^-5 then n = 4 and N = 49996.__

Then Odds ratios are calculated using Fisher's exact test (one sided to test for enrichment only).

## NEW METHOD

To compensate for the difference between extrapolating from MAF and using exact number of individuals from public databases, TAPES now first calculate OR the old method. 

Then another OR is calculated using __n = ceil(y/3)*10 and N = (((x/2) - n)*10) - n__

The final OR is a mean between the OLD and NEW OR calculation.

## Graphs
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

![PS4 calculation](https://raw.githubusercontent.com/a-xavier/tapes/master/toy_dataset/New%20PS4%20calc.png "OR calculation")