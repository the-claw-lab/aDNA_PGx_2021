# Pharmacogenetic variation in Neanderthals and Denisovans and the implications for evolutionary medicine

## Genotyping pipeline

**Samples were genotyped using using [pypgx - version 0.1.37 *bam2vcf*](https://github.com/sbslee/pypgx)**

`CODE BLOCK`

**Genotype and read depth files were inputted into Stargazer with the following command**

`python3 stargazer.py "genotype"
	-o "output_preffix"
	-d wgs
	-t "supported_gene"
	--vcf "input.vcf"
	-c "control_gene #can be vdr, ryr1, or egfr
	--gdf "input_rd.gdf"`

