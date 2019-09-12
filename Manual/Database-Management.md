# ```db``` and ANNOVAR Database Management

db is the option to manage your ANNOVAR databases.
First you need to [download ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/download/).

## First use 
If you use TAPES to manage ANNOVAR databases for the first time you need to provide the location of your ANNOVAR folder using:  
```python3 tapes.py db -s -A /path/to/annovar/```

This will create 2 useful files in the src folder : **db_config.json** and **db_vcf.json** in the src folder.

## Simplified database download
You can download only the databases needed by TAPES using:  
```python3 tapes.py db -b --acmg -a hg19```

## Full database management

Using the file db_config.json, you can see which databases are downloaded and which are missing.
To flag a database for download, open **db_config.json** using any editor and replace "**MISSING**" by "**DOWNLOAD**" or "**DOWN**".
then run:  
```python3 tapes.py db -b -a hg38```  

The databases flagged for download will all be downloaded

## Options 

* ```-s``` ```--see_db```|_flag_|See which databases are present in the ANNOVAR folder  
* ```-b``` ```--build_db```|_flag_| Download flagged databases 
* ```--acmg``` |_flag_| Bypasses the db_config.json file and download only necessary databases for TAPES
* ```-a```  ```--assembly```|_string_| Either hg19 or hg38 - default is hg19

