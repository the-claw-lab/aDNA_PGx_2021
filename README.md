# Pharmacogenetic variation in Neanderthals and Denisovans and the implications for evolutionary medicine

## Genotyping pipeline

**Samples were genotyped and read depth was extracted using using [pypgx - version 0.1.37](https://github.com/sbslee/pypgx)**

`pypgx bam2vcf gatk "reference_genome" "gene" "output_name.vcf" hg19 --bam_dir "${d}/bam_Files.txt"`

`pypgx bam2gdf hg19 "gene" "control_gene" "output_name.gdf" --bam_list "${d}/bamFiles.txt"`

## Analysis pipeline

**Genotype and read depth files were inputted into [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/) [(Lee et al. 2019)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.1552) with the following command**

`python3 stargazer.py "genotype" -o "output_preffix" -d wgs -t "supported_gene" --vcf "input.vcf" -c "control_gene --gdf "input_rd.gdf"`

**Stargazer**

Output Processing: [2021_Archaic_PGx_SGE_Output_Analysis.R]()

**SNV Analysis**

Select gene positions (start-stop gene loci + 2000bp upstream for the promoter) and filter the QUAL >= 40

`bcftools view -t "chr:strt-end" "input.vcf" > "output.vcf"`

`bcftools view -i '%QUAL>=40' "output.vcf" > "output_filtered40.vcf"`

SNVs are annotated with [Annovar](), [CADD](), [SIFT](), and [PolyPhen3]()

SNV Analysis Script: [2021_Archaic_PGx_variantSNV_Analysis.R]()

**Potentially Damaging Variants**
