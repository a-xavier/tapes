# ```decompose```

decompose allow users to manually decompose their vcf files prior to annotation.  
Decomposing vcf files allow remove multi-alleles (eg. REF-A ALT-T,C,GT) loci from vcf files to have a 'one variant per line' file.

# Usage

```python tapes.py decompose -i ./input.vcf -o ./output_decomposed.vcf```

# Credits

The ```decompose``` function is using the [vcf_parser](https://github.com/moonso/vcf_parser) module.