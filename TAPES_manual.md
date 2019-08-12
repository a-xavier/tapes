**TAPES - INSTRUCTION MANUAL**

TAPES: a Tool for Assessment and Prioritisation in Exome Studies, is a script written in python 3.7 which serves three purposes:

1. Be a simplified interface to ANNOVAR ([http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/)) with easy database management and easy commands for annotation.

1. Prioritize variants using the ACMG 2015 (DOI: 10.1038/gim.2015.30) criteria and probability of pathogenicity to classify variants from pathogenic to benign.

1. Create appropriate reports for researcher based on relevant criteria.

TAPES focuses on multi-sample VCFs files and disease cohorts but any file annotated with ANNOVAR can be used.

COMPATIBILITY

TAPES used Travis CI, testing on Linux, Windows and macOS.

TAPES main function: sort will work on both UNIX and Windows.

The ANNOVAR interface will only work on UNIX due to ANNOVAR compatibility.

TAPES was written and tested on python3.7 and will work on any python3 version.

Setting up a virtual environment using python3 on Linux or windows: 

 - ```cd``` to the TAPES directory
 - ```python3 -m venv tapes_env``` creating a virtual environment called 'tapes_env'
 - ```source env/bin/activate``` Activate the virtual environment (```.\env\Scripts\activate``` on windows)
 - ```pip3 install -r requirements.txt --user``` to install dependencies (winrequirements.txt on windows)
 - [Use TAPES](https://github.com/a-xavier/tapes/wiki/Quick-Start)
 - ```deactivate``` to leave the virtual environment

1. 1)QUICK START

Choose a starting point:

- --Unannotated VCF file
  - I just want to use the necessary ANNOVAR annotations for TAPES sorting
    - python tapes.py db -s -A /path/to/annovar/
    - python tapes.py db -b --acmg --assembly hg19
    - python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.vcf --acmg –a hg19
  - I want to customise the ANNOVAR annotations used
    - See detailed instructions on section **ANNOVAR INTERFACE**
- --VEP or ANNOVAR annotated file
  - I want to sort and prioritise my VCF file
    - python tapes.py sort -i /to/annotated/file.vcf –o /to/output/folder/ --tab
  - I want to re-analyse TAPES main report
    - python tapes.py analyse -i /to/main/report.txt –o /to/output/report.txt –-single\_option





Workflow Example

Here you can find an example for workflow starting from a **multi-sample VCF** file and a fresh [ANNOVAR download](http://www.openbioinformatics.org/annovar/annovar_download_form.php)

First use:
python tapes.py db -s -A /path/to/annovar/
To determine the ANNOVAR folder.

Then:
python tapes.py db -b --acmg
To download necessary databases for TAPES (this might take some time).

Then:
python tapes.py annotate -i /path/to/multi-sample.vcf -o /path/to/annotated.vcf --acmg
To annotate the original vcf.

Then either:

- python tapes.py sort -i /path/to/annotated.vcf -o /path/to/output/ --tab --acmg
To prioritise the annotated vcf and generate the main output.

**or**

- python tapes.py sort -i /path/to/annotated.vcf -o /path/to/output/ --tab --acmg --by\_gene --by\_sample --disease cancer --enrichr --list &quot;TP53 BRCA1 BRCA2 PTEN KRAS&quot;

To prioritise the annotated vcf, generate the main output and other few secondary reports

1) INSTALLATION

TAPES does not require installation, just download the repository at https://github.com/a-xavier/tapes and extract it to any convenient location.

If pip is not installed on your system you can install it easily:

Install PIP On Debian/Ubuntu

apt install python3-pip

Install PIP on Fedora

dnf install python3-pip

Install PIP on Arch Linux

pacman -S python-pip

Install PIP on Windows

First install python3 from [https://www.python.org/downloads/](https://www.python.org/downloads/) and add python and pip to your path in the environment variable menu.

(On windows 7 :  Control Panel -\&gt; System - \&gt; Advanced System Settings -\&gt; Environment variables  then  under System Variables double click on path and add the python installation path separated by a semicolon &quot; ; &quot; )

Then use either cmd.exe or Windows Powershell to use TAPES.

Using pip you can install all the requirement with:

cd path/to/TAPES

pip install --upgrade -r requirements.txt



This will install all the required python modules.

On windows, use:

pip install --upgrade -r winrequirements.txt

Note that since TAPES is written in python 3, you might need to run it using python3 instead of python depending on your system.

If you plan on using TAPES as an ANNOVAR wrapper, please download ANNOVAR first here :

[http://annovar.openbioinformatics.org/en/latest/](http://annovar.openbioinformatics.org/en/latest/)

2) RUNNING A PRIORITISATION JOB

1) Prioritise the annotated file

To prioritise your variant use the sort option.

There is two main output mode: FOLDER and CSV/TXT+XLSX

When writing the output, just specify a folder or a csv file to choose the mode (see examples below).

In both mode, the flag --acmg can be added.

- Using the --acmg tag will ensure all the main annotations for ACMG classification are present before the sorting process. If you are not sure your annotated file in fully compliant with TAPES, you can remove the –acmg tag.
- If the --acmg tag is not present, TAPES will annotate as much as it can based on the present annotation. This ensures that even older files annotated with ANNOVAR can be prioritised to a certain extent.

        a) FOLDER MODE

This mode will output a folder with different csv files and figures based on the options:

python tapes.py sort -i /path/to/annotated/file.csv –o/path/to/output/folder/

Will output csv files

python tapes.py sort -i /to/annotated/file.csv –o /to/output/folder/ --tab

Will output tab-separated files

The output must be either an empty folder or a non-existent folder.

        b) CSV/TXT+XLSX MODE

This mode will output a csv file and an xlsx report containing different spreadsheets based on the options:

python tapes.py sort -i /path/to/annotated/file.csv -o /path/to/output.csv

will output csv + xlsx files

python tapes.py sort -i /path/to/annotated/file.csv -o /path/to/output.txt

will output a tab-separated txt + xlsx files (.tsv also works)



2) Sorting and reporting Options

| **Option        ** | **Type** | **Description** | **Default** |
| --- | --- | --- | --- |
| --acmg | flag | Perform check for main annotations before sorting |   |
| --trio | Path to txt file | A trio text file (see specification) |   |
| --by\_sample | flag | Create output with the 5 most pathogenic variants per sample |   |
| --enrichr | str | Use enrichr to analayse the pathways impacted by pathogenic variants | GO\_Biological\_Process\_2018 |
| --disease | str | Check in the &#39;disease&#39; column the presence of a term | cancer |
| --list | str or path to txt file | A list of gene of interest (in quotes separated by a space) or a text file with one gene symbol per line |   |
| --kegg | str | Similar to list but when you do not know all genes of interest. Select a pathways and a report will be created with only genes involved in that pathway (see Appendix for the full list of available pathways) |   |
| --by\_gene | flag | Create output ranking each gene based on a simple Gene-burden metrics |   |



Notes on --trio :

The trio file must be a tab delimited file with the following info:

Role in family: **m f o** in no particular order for mother, father and offspring

Trio id: any string without space

Sample name: as they appear on the original vcf file

Only use UNIQUE trio IDs; if there are several trios in one family, use different IDs or the result will be incorrect.

Note on --by\_gene flag:

        The --by\_gene flag will create a report grouping all variants **that are predicted to be pathogenic** contained in a single gene. The metrics used to measure gene burden is quite simple:   where  is the probability for a variant to be pathogenic and  the number of Sample affected by this variant.

Since this score does not account for several other parameters, a number of warnings are also present:

- --Number of sample warning: If more than half of the variant of each gene are present in more than half of the samples. It means that the number of sample affected is suspiciously high. This can happened in misaligned reads in X and Y homologous regions for example.
- --Long gene warning: If the gene is long (more than 250,000 bp), more variants are expected.
- --FLAGS Gene: FLAGS genes are the most frequently mutated genes in Exome sequencing. See [https://doi.org/10.1186/s12920-014-0064-y](https://doi.org/10.1186/s12920-014-0064-y) for more details.



3) Output Explained

Please note that apart from the main report, all other reports will have suffixes appended to the main report name such as &quot;Report\_by\_sample.txt&quot; or &quot;Report\_GO\_Biological\_Process\_2018.txt&quot;



1. a)Main Output

The ouptut files will always be sorted csv/txt/tsv or xlsx files. The variants are sorted from most pathogenic to most benign. Apart from the classical ACMG classification (see original paper for infos, S Richards et al - ‎2015), TAPES will also provide an estimated probability of pathogenicity calculated based on S.V. Tavtigian et al 2018. To be simple it outputs the probability that this particular variant is pathogenic based on the ACMG criteria.

The default of Prior\_P = 0.1, exponent X = 2 and OPVST=350 are used.



1.
b)

By-Sample report


By-Sample report

This report will contain the 5 most pathogenic variants per sample.


Eg.

.

1.
c)

By-Gene report


By-Gene report

Every table has, above the header, the name of the gene, the gene burden score and, in certain cases, a warning.



1.
d)

EnrichR report


EnrichR report

The 11 most relevant pathway will be in the EnrichR report. Only pathways with significant **adjusted** p-values should be considered





1. e)Kegg, List and Disease reports

Kegg, list and Disease report will look very similar to the main output. Kegg and list will only show variant that belong to either a determined keg pathway (see list in appendix) or a list of user-provided genes.

The disease report will only show variant that have a certain term in the Disease column of the annotation. Eg. &quot;Autosomal dominant&quot;, &quot;cancer&quot;, &quot;Colorectal&quot;

3) ANNOVAR INTERFACE

 Note that TAPES accepts for annotation: vcf files, bcf files, bgzipped vcf files and gzipped bcf files. They will automatically converted to vcf files prior to annotation.

 Users should also have downloaded ANNOVAR first (free for non-commercial use) : [http://www.openbioinformatics.org/annovar/annovar\_download\_form.php](http://www.openbioinformatics.org/annovar/annovar_download_form.php)

1) First Use

When using TAPES for the first time, you need to indicate the location of your local ANNOVAR folder:

python tapes.py db -s -A /path/to/annovar/


The -s stands for --see-db, a tag used to see all databases present on your system. The output should look like this:

2) Simplified database management and annotation: Using the --acmg tag



        a) DOWNLOADING DATABASES

Use db -b --acmg or db --build\_db --acmg to start downloading the necessary databases for the ACMG criteria assignment. You can specify the assembly to use (either hg19 or hg38) with the --assembly option (default is hg19)

The necessary databases for all possible criteria assignment are:

- gnomad\_genome
- gnomad\_exome or exac03 (gnomad\_exome is the default)
- avsnp150
- clinvar\_20180603
- dbnsfp35c
- ne of the genome annotation : refGene, ensGene, knownGene

python tapes.py db -b --acmg --assembly hg19

This command will download the databases in the /humandb directory located in the ANNOVAR folder.

You can then check that all the databases have been downloaded using:

python tapes.py db –s

        b) ANNOTATING VCF FILE

To annotate a VCF file, use the annotate option with --acmg tag to easily annotate your vcf with all the relevant databases for ACMG classification. One again use --assembly to specify the assembly version

python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.csv --acmg –assembly hg19

This will produce the annotated file **output.csv** and if the vcf is multi-sample, the file **output\_with\_samples.csv** will also be created.

python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.txt --acmg –assembly hg19

This will produce the annotated file **output.txt** and if the vcf is multi-sample, the file **output\_with\_samples.txt** will also be created.



python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.vcf --acmg –assembly hg19

This will produce the annotated file **output.vcf.**



3)  Advanced database management and annotation

        a) DATABASE MANAGMENT

TAPES provides two files to easily manage databases and ANNOVAR annotations.

db\_config.json is an easily readable json file which shows all (most of the) available ANNOVAR databases.

Those files are generated after the first use.

Missing databases are flagged &quot;MISSING&quot;, downloaded databases are flagged OK.

To flag a database for download, replace &quot;MISSING&quot; by &quot;DOWNLOAD&quot; or &quot;DOWN&quot;.

Then run:

python tapes.py db -b

This will download all databases flagged for download.

        b) ANNOTATION

db\_vcf.json is an easily readable json file which shows all downloaded databases and which databases are used to annotate vcf\_files.

Databases flagged &quot;YES&quot; will be used for annotation and databases flagged &quot;NO&quot; will be ignored.

Flag &quot;YES&quot; for all databases you want to use for annotation then run:

 python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.csv

This will output two file: a standard annotated output.csv file and an output\_with\_samples.csv containing sample genotyping data.

python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.txt

This will output two file: a standard annotated output.txt file and an output\_with\_samples.txt containing sample genotyping data.

python tapes.py annotate -i /path/to/file.vcf -o /path/to/output.vcf

This will produce the annotated file **output.vcf**



        c) ANNOTATION OPTIONS

| **Option** | **Type** | **Description** | **Default** |
| --- | --- | --- | --- |
| --assembly | str | Assembly version : either hg19 or hg38 | hg19 |
| --ref\_anno | str | Genome annotation : either refGene for RefSeq,ensGene for ENSEMBL qnd knownGene for UCSC | refGene |

4) DECOMPOSING VCF

TAPES will automatically decompose VCFs files before annotation. But TAPES can decompose a VCF file without annotating it using:

python tapes.py decompose –i /original.vcf –o /decomposed.vcf

5) RE-ANALYSING TAPES OUTPUTS

        If you want to generate a report from previously sorted file. You can use the analyse (or analyze) option.

        For example:

python tapes.py analyse -i /path/to/sorted\_output.txt -o /path/to/output\_report.txt  --by\_sample

Will output a by-sample report

python tapes.py analyse -i /path/to/sorted\_output.txt -o /path/to/output\_report.txt  --by\_gene

Will output a by-gene report

Please note that you can only output one report at a time. For example

python tapes.py analyse -i /path/to/sorted\_output.txt -o /path/to/output\_report.txt –--by\_gene –by-sample –enrichr –list &quot;MLH1 MSH2 APC&quot;

will not work.













APPENDIX

KEGG Pathways keys

- 2-oxocarboxylic acid metabolism
- abc transporters
- acute myeloid leukemia
- adherens junction
- adipocytokine signaling pathway
- adrenergic signaling in cardiomyocytes
- african trypanosomiasis
- age-rage signaling pathway in diabetic complications
- alanine, aspartate and glutamate metabolism
- alcoholism
- aldosterone synthesis and secretion
- aldosterone-regulated sodium reabsorption
- allograft rejection
- alpha-linolenic acid metabolism
- alzheimer disease
- amino sugar and nucleotide sugar metabolism
- aminoacyl-trna biosynthesis
- amoebiasis
- amphetamine addiction
- ampk signaling pathway
- amyotrophic lateral sclerosis
- antifolate resistance
- antigen processing and presentation
- apelin signaling pathway
- apoptosis
- apoptosis - multiple species
- arachidonic acid metabolism
- arginine and proline metabolism
- arginine biosynthesis
- arrhythmogenic right ventricular cardiomyopathy
- ascorbate and aldarate metabolism
- asthma
- autoimmune thyroid disease
- autophagy - animal
- autophagy - other
- axon guidance
- b cell receptor signaling pathway
- bacterial invasion of epithelial cells
- basal cell carcinoma
- basal transcription factors
- base excision repair
- beta-alanine metabolism
- bile secretion
- biosynthesis of amino acids
- biosynthesis of unsaturated fatty acids
- biotin metabolism
- bladder cancer
- breast cancer
- butanoate metabolism
- c-type lectin receptor signaling pathway
- caffeine metabolism
- calcium signaling pathway
- camp signaling pathway
- carbohydrate digestion and absorption
- carbon metabolism
- cardiac muscle contraction
- cell adhesion molecules
- cell cycle
- cellular senescence
- central carbon metabolism in cancer
- cgmp-pkg signaling pathway
- chagas disease
- chemical carcinogenesis
- chemokine signaling pathway
- cholesterol metabolism
- choline metabolism in cancer
- cholinergic synapse
- chronic myeloid leukemia
- circadian entrainment
- circadian rhythm
- citrate cycle
- cocaine addiction
- collecting duct acid secretion
- colorectal cancer
- complement and coagulation cascades
- cortisol synthesis and secretion
- cushing syndrome
- cysteine and methionine metabolism
- cytokine-cytokine receptor interaction
- cytosolic dna-sensing pathway
- d-arginine and d-ornithine metabolism
- d-glutamine and d-glutamate metabolism
- dilated cardiomyopathy
- dna replication
- dopaminergic synapse
- drug metabolism - cytochrome p450
- drug metabolism - other enzymes
- ecm-receptor interaction
- egfr tyrosine kinase inhibitor resistance
- endocrine and other factor-regulated calcium reabsorption
- endocrine resistance
- endocytosis
- endometrial cancer
- epithelial cell signaling in helicobacter pylori infection
- epstein-barr virus infection
- erbb signaling pathway
- estrogen signaling pathway
- ether lipid metabolism
- fanconi anemia pathway
- fat digestion and absorption
- fatty acid biosynthesis
- fatty acid degradation
- fatty acid elongation
- fatty acid metabolism
- fc epsilon ri signaling pathway
- fc gamma r-mediated phagocytosis
- ferroptosis
- fluid shear stress and atherosclerosis
- focal adhesion
- folate biosynthesis
- foxo signaling pathway
- fructose and mannose metabolism
- gabaergic synapse
- galactose metabolism
- gap junction
- gastric acid secretion
- gastric cancer
- glioma
- glucagon signaling pathway
- glutamatergic synapse
- glutathione metabolism
- glycerolipid metabolism
- glycerophospholipid metabolism
- glycine, serine and threonine metabolism
- glycolysis / gluconeogenesis
- glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate
- glycosaminoglycan biosynthesis - heparan sulfate / heparin
- glycosaminoglycan biosynthesis - keratan sulfate
- glycosaminoglycan degradation
- glycosphingolipid biosynthesis - ganglio series
- glycosphingolipid biosynthesis - globo and isoglobo series
- glycosphingolipid biosynthesis - lacto and neolacto series
- glycosylphosphatidylinositol
- glyoxylate and dicarboxylate metabolism
- gnrh signaling pathway
- graft-versus-host disease
- hedgehog signaling pathway
- hematopoietic cell lineage
- hepatitis b
- hepatitis c
- hepatocellular carcinoma
- herpes simplex infection
- hif-1 signaling pathway
- hippo signaling pathway
- hippo signaling pathway - multiple species
- histidine metabolism
- homologous recombination
- human cytomegalovirus infection
- human immunodeficiency virus 1 infection
- human papillomavirus infection
- human t-cell leukemia virus 1 infection
- huntington disease
- hypertrophic cardiomyopathy
- il-17 signaling pathway
- inflammatory bowel disease
- inflammatory mediator regulation of trp channels
- influenza a
- inositol phosphate metabolism
- insulin resistance
- insulin secretion
- insulin signaling pathway
- intestinal immune network for iga production
- jak-stat signaling pathway
- kaposi sarcoma-associated herpesvirus infection
- legionellosis
- leishmaniasis
- leukocyte transendothelial migration
- linoleic acid metabolism
- lipoic acid metabolism
- long-term depression
- long-term potentiation
- longevity regulating pathway
- longevity regulating pathway - multiple species
- lysine degradation
- lysosome
- malaria
- mannose type o-glycan biosynthesis
- mapk signaling pathway
- maturity onset diabetes of the young
- measles
- melanogenesis
- melanoma
- metabolic pathways
- metabolism of xenobiotics by cytochrome p450
- micrornas in cancer
- mineral absorption
- mismatch repair
- mitophagy - animal
- morphine addiction
- mrna surveillance pathway
- mtor signaling pathway
- mucin type o-glycan biosynthesis
- n-glycan biosynthesis
- natural killer cell mediated cytotoxicity
- necroptosis
- neomycin, kanamycin and gentamicin biosynthesis
- neuroactive ligand-receptor interaction
- neurotrophin signaling pathway
- nf-kappa b signaling pathway
- nicotinate and nicotinamide metabolism
- nicotine addiction
- nitrogen metabolism
- nod-like receptor signaling pathway
- non-alcoholic fatty liver disease
- non-homologous end-joining
- non-small cell lung cancer
- notch signaling pathway
- nucleotide excision repair
- lfactory transduction
- ne carbon pool by folate
- cyte meiosis
- steoclast differentiation
- ther glycan degradation
- ther types of o-glycan biosynthesis
- varian steroidogenesis
- xidative phosphorylation
- xytocin signaling pathway
- p53 signaling pathway
- pancreatic cancer
- pancreatic secretion
- pantothenate and coa biosynthesis
- parathyroid hormone synthesis, secretion and action
- parkinson disease
- pathogenic escherichia coli infection
- pathways in cancer
- pentose and glucuronate interconversions
- pentose phosphate pathway
- peroxisome
- pertussis
- phagosome
- phenylalanine metabolism
- phenylalanine, tyrosine and tryptophan biosynthesis
- phosphatidylinositol signaling system
- phospholipase d signaling pathway
- phosphonate and phosphinate metabolism
- phototransduction
- pi3k-akt signaling pathway
- platelet activation
- platinum drug resistance
- porphyrin and chlorophyll metabolism
- ppar signaling pathway
- primary bile acid biosynthesis
- primary immunodeficiency
- prion diseases
- progesterone-mediated oocyte maturation
- prolactin signaling pathway
- propanoate metabolism
- prostate cancer
- proteasome
- protein digestion and absorption
- protein export
- protein processing in endoplasmic reticulum
- proteoglycans in cancer
- proximal tubule bicarbonate reclamation
- purine metabolism
- pyrimidine metabolism
- pyruvate metabolism
- rap1 signaling pathway
- ras signaling pathway
- regulation of actin cytoskeleton
- regulation of lipolysis in adipocytes
- relaxin signaling pathway
- renal cell carcinoma
- renin secretion
- renin-angiotensin system
- retinol metabolism
- retrograde endocannabinoid signaling
- rheumatoid arthritis
- riboflavin metabolism
- ribosome
- ribosome biogenesis in eukaryotes
- rig-i-like receptor signaling pathway
- rna degradation
- rna polymerase
- rna transport
- salivary secretion
- salmonella infection
- selenocompound metabolism
- serotonergic synapse
- shigellosis
- signaling pathways regulating pluripotency of stem cells
- small cell lung cancer
- snare interactions in vesicular transport
- sphingolipid metabolism
- sphingolipid signaling pathway
- spliceosome
- staphylococcus aureus infection
- starch and sucrose metabolism
- steroid biosynthesis
- steroid hormone biosynthesis
- sulfur metabolism
- sulfur relay system
- synaptic vesicle cycle
- synthesis and degradation of ketone bodies
- systemic lupus erythematosus
- t cell receptor signaling pathway
- taste transduction
- taurine and hypotaurine metabolism
- terpenoid backbone biosynthesis
- tgf-beta signaling pathway
- th1 and th2 cell differentiation
- th17 cell differentiation
- thermogenesis
- thiamine metabolism
- thyroid cancer
- thyroid hormone signaling pathway
- thyroid hormone synthesis
- tight junction
- tnf signaling pathway
- toll-like receptor signaling pathway
- toxoplasmosis
- transcriptional misregulation in cancer
- tryptophan metabolism
- tuberculosis
- type i diabetes mellitus
- type ii diabetes mellitus
- tyrosine metabolism
- ubiquinone and other terpenoid-quinone biosynthesis
- ubiquitin mediated proteolysis
- valine, leucine and isoleucine biosynthesis
- valine, leucine and isoleucine degradation
- vascular smooth muscle contraction
- vasopressin-regulated water reabsorption
- vegf signaling pathway
- vibrio cholerae infection
- viral carcinogenesis
- viral myocarditis
- vitamin b6 metabolism
- vitamin digestion and absorption
- wnt signaling pathway



EnrichR Libraries

- Genes\_Associated\_with\_NIH\_Grants
- Cancer\_Cell\_Line\_Encyclopedia
- Achilles\_fitness\_decrease
- Achilles\_fitness\_increase
- Aging\_Perturbations\_from\_GEO\_down
- Aging\_Perturbations\_from\_GEO\_up
- Allen\_Brain\_Atlas\_down
- Allen\_Brain\_Atlas\_up
- ARCHS4\_Cell-lines
- ARCHS4\_IDG\_Coexp
- ARCHS4\_Kinases\_Coexp
- ARCHS4\_TFs\_Coexp
- ARCHS4\_Tissues
- BioCarta\_2013
- BioCarta\_2015
- BioCarta\_2016
- BioPlex\_2017
- ChEA\_2013
- ChEA\_2015
- ChEA\_2016
- Chromosome\_Location
- Chromosome\_Location\_hg19
- CORUM
- Data\_Acquisition\_Method\_Most\_Popular\_Genes
- dbGaP
- Disease\_Perturbations\_from\_GEO\_down
- Disease\_Perturbations\_from\_GEO\_up
- Disease\_Signatures\_from\_GEO\_down\_2014
- Disease\_Signatures\_from\_GEO\_up\_2014
- Drug\_Perturbations\_from\_GEO\_2014
- Drug\_Perturbations\_from\_GEO\_down
- Drug\_Perturbations\_from\_GEO\_up
- DrugMatrix
- DSigDB
- ENCODE\_and\_ChEA\_Consensus\_TFs\_from\_ChIP-X
- ENCODE\_Histone\_Modifications\_2013
- ENCODE\_Histone\_Modifications\_2015
- ENCODE\_TF\_ChIP-seq\_2014
- ENCODE\_TF\_ChIP-seq\_2015
- Enrichr\_Libraries\_Most\_Popular\_Genes
- Enrichr\_Submissions\_TF-Gene\_Coocurrence
- Epigenomics\_Roadmap\_HM\_ChIP-seq
- ESCAPE
- GeneSigDB
- Genome\_Browser\_PWMs
- GO\_Biological\_Process\_2013
- GO\_Biological\_Process\_2015
- GO\_Biological\_Process\_2017
- GO\_Biological\_Process\_2017b
- GO\_Biological\_Process\_2018
- GO\_Cellular\_Component\_2013
- GO\_Cellular\_Component\_2015
- GO\_Cellular\_Component\_2017
- GO\_Cellular\_Component\_2017b
- GO\_Cellular\_Component\_2018
- GO\_Molecular\_Function\_2013
- GO\_Molecular\_Function\_2015
- GO\_Molecular\_Function\_2017
- GO\_Molecular\_Function\_2017b
- GO\_Molecular\_Function\_2018
- GTEx\_Tissue\_Sample\_Gene\_Expression\_Profiles\_down
- GTEx\_Tissue\_Sample\_Gene\_Expression\_Profiles\_up
- HMDB\_Metabolites
- HomoloGene
- Human\_Gene\_Atlas
- Human\_Phenotype\_Ontology
- HumanCyc\_2015
- HumanCyc\_2016
- huMAP
- Jensen\_COMPARTMENTS
- Jensen\_DISEASES
- Jensen\_TISSUES
- KEA\_2013
- KEA\_2015
- KEGG\_2013
- KEGG\_2015
- KEGG\_2016
- Kinase\_Perturbations\_from\_GEO\_down
- Kinase\_Perturbations\_from\_GEO\_up
- Ligand\_Perturbations\_from\_GEO\_down
- Ligand\_Perturbations\_from\_GEO\_up
- LINCS\_L1000\_Chem\_Pert\_down
- LINCS\_L1000\_Chem\_Pert\_up
- LINCS\_L1000\_Kinase\_Perturbations\_down
- LINCS\_L1000\_Kinase\_Perturbations\_up
- LINCS\_L1000\_Ligand\_Perturbations\_down
- LINCS\_L1000\_Ligand\_Perturbations\_up
- MCF7\_Perturbations\_from\_GEO\_down
- MCF7\_Perturbations\_from\_GEO\_up
- MGI\_Mammalian\_Phenotype\_2013
- MGI\_Mammalian\_Phenotype\_2017
- MGI\_Mammalian\_Phenotype\_Level\_3
- MGI\_Mammalian\_Phenotype\_Level\_4
- Microbe\_Perturbations\_from\_GEO\_down
- Microbe\_Perturbations\_from\_GEO\_up
- miRTarBase\_2017
- Mouse\_Gene\_Atlas
- MSigDB\_Computational
- MSigDB\_Oncogenic\_Signatures
- NCI-60\_Cancer\_Cell\_Lines
- NCI-Nature\_2015
- NCI-Nature\_2016
- NURSA\_Human\_Endogenous\_Complexome
- Old\_CMAP\_down
- Old\_CMAP\_up
- OMIM\_Disease
- OMIM\_Expanded
- Panther\_2015
- Panther\_2016
- Pfam\_InterPro\_Domains
- Phosphatase\_Substrates\_from\_DEPOD
- PPI\_Hub\_Proteins
- Reactome\_2013
- Reactome\_2015
- Reactome\_2016
- RNA-Seq\_Disease\_Gene\_and\_Drug\_Signatures\_from\_GEO
- SILAC\_Phosphoproteomics
- Single\_Gene\_Perturbations\_from\_GEO\_down
- Single\_Gene\_Perturbations\_from\_GEO\_up
- SysMyo\_Muscle\_Gene\_Sets
- TargetScan\_microRNA
- TargetScan\_microRNA\_2017
- TF-LOF\_Expression\_from\_GEO
- TF\_Perturbations\_Followed\_by\_Expression
- Tissue\_Protein\_Expression\_from\_Human\_Proteome\_Map
- Tissue\_Protein\_Expression\_from\_ProteomicsDB
- Transcription\_Factor\_PPIs
- TRANSFAC\_and\_JASPAR\_PWMs
- Virus\_Perturbations\_from\_GEO\_down
- Virus\_Perturbations\_from\_GEO\_up
- VirusMINT
- WikiPathways\_2013
- WikiPathways\_2015
- WikiPathways\_2016

ACMG Criteria assignment (refer to S Richards et al - ‎2015 for a description of the criteria)

Pathogenic Criteria

PVS1

Will be assigned to a variant if it is a stopgain or frameshift deletion/insertion located 50 bp further than the end of the final exon. (Based on the ExonicFunc column anf the REK\_canon library)

Will be assigned to a splicing variant with a dbscSNV score of more than 0.6 (ADA or RF) (Based on the Func column and the dbscSNV score annotation)

PS1

Will be assigned if a variant have the same AA ref and AA alt as a known pathogenic variant.

Using all known pathogenic variants from clinvar

PS2

Will be assigned if a variant is assumed de novo and parents are disease free. This requires trio data.

PS3

Will be assigned if clinvar classifies the variant as Pathogenic or drug reponse and the level of evidence is either &#39;practice guideline&#39; or &#39;reviewed by expert panel&#39;

PS4

Will be assigned if a variant is enriched in the samples provided. Requires either &#39;output\_with\_samples.csv&#39; from the annotation to keep sample genotyping data or an annotated multi-sample vcf. PS4 will take the affected individuals with the mutations and the total number of individuals in the disease cohort and compare it to the data from gnomad\_genome and gnomad\_exome.

The number of individuals with and without variants in public data is extrapolated with the following formula:

Minor allele frequency in control population (MAF) =

Number of individuals with the variant in control population =

Total number of individuals in control population =

Then a fisher&#39;s exact test is performed to calculate the odd ratios, the confidence interval and the p value.

PS4 will only be considered if at least 2 samples are affected by a variant. Otherwise, Intervar PS4 database, based on GWAS database will be used.

PS4 will be assigned if the Odd Ratio is superior to 20, the confidence interval does not cross one and the p value is under 0.01

PM1

Will be assigned if the variant is a Missense variant (nonsynonymous SNV) and is located in a in a domain without benign variants (Using Intervar db) for benign domains

PM2

Will be assigned if the variant is in a recessive gene and has a frequency under 0.005 or is in a dominant gene and has no frequency data available. Recessive and Dominant/Haploinsufficient genes were infered using Pli and Prec scores computed by Lek et al, 2016. A gene is considered dominant dominant with a pli \&gt;0.85 and recessive if prec \&gt;0.85

PM4

Will be assigned if the variant is an in-frame deletion/insertion in a non-repeat region of the gene. Using the repeat\_dict database.

PM5

Will be assigned if a variant have the same AA ref and  a different AA alt as a known pathogenic variant.

Using all known pathogenic variants from clinvar

PP2

Will be assigned if the variant is Missense (nonsynonymous SNV) in a gene where missense variants represents at least 80 percent of all known pathogenic variants (using PP2\_BP1 database)

PP3

Will be assigned if the variant is predicted to be pathogenic using various in-silico prediction tools (sift, lrt, mutationtaster, mutation assessor, fathmm, provean, meta svm, meta lr, mcap, mkl, genocanyon, gerp)

PP5

Will be assigned the variant is classified as pathogenic or likely pathogenic by clinvar but the evidence is limited.

Benign criteria

BA1

Will be assigned to a variant if its frequency in gnomad\_exome/exac or gnomad\_genome is superior to 0.05

BS1

Will be assigned to a variant if its frequency is superior to a cutoff (0.005) for a rare disease.

BS2

Will be assigned if the variant was observed in a healthy individual as homozygous for a recessive disease and heterozygous for a dominant disease. (Using Intervar db BS2\_hom\_het)

BS3

Will be assigned if clinvar classifies the variant as Benign or likely benign and the level of evidence is either &#39;practice guideline&#39; or &#39;reviewed by expert panel&#39;

BP1

Will be assigned if the variant is Missense (nonsynonymous SNV) in a gene where missense variants represents at most 10 percent of all known pathogenic variants (using PP2\_BP1 database).

BP3

Will be assigned if the variant is an in-frame deletion/insertion in a repeat region of the gene. (Using the repeat\_dict database).

BP4

Will be assigned if the variant is predicted to be benign using various in-silico prediction tools (sift, lrt, mutationtaster, mutation assessor, fathmm, provean, meta svm, meta lr, mcap, mkl, genocanyon, gerp)

BP6

Will be assigned the variant is classified as Benign or likely benign by clinvar but the evidence is limited.

BP7

Will be assigned if a variant if synonymous and no splicing impact is predicted by dbscSNV (score under 0.6)
