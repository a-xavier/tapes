# --disease  
#### Used to filter the main report based on the associated disease description.

usage:  
```python3 tapes.py sort -i ./input.vcf -o ./output_folder/ --disease "autosomal recessive"```  
```python3 tapes.py sort -i ./input.vcf -o ./output_folder/ -d "folate metabolism" ```

By using the option ```--disease``` or ```-d``` with the ```sort``` or ```analyse``` option, a report will be created containing only variant with the search term in their "disease" column. 

User can use a variety of terms such as "autosomal dominant", "multiple sclerosis", "cancer", etc.

Quotes are expected if the search term contains multiple words. 

The format of this report will be identical to the main report.