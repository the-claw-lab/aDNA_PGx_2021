This pipeline was developed to identify pharmacogenetic variation in publicly avalailable archaic genomes. Contact Teddy Wroblewski (tadeusz.wroblewski@cuanschutz.edu) for more information.

## Required Software

- RStudio
- Python3
- pypgx
- Stargazer
- bcftools

## Genotyping pipeline

**Samples were genotyped and read depth was extracted using using [pypgx - version 0.1.37](https://github.com/sbslee/pypgx)**

`pypgx bam2vcf gatk "reference_genome" "gene" "output_name.vcf" hg19 --bam_dir "${d}/bam_Files.txt"`

`pypgx bam2gdf hg19 "gene" "control_gene" "output_name.gdf" --bam_list "${d}/bamFiles.txt"`

Examples for *CYP1A2* and *CYP2A6* can be found in `input/genotype/` and `input/depth_of_coverage/`

## Analysis pipeline

**Genotype and read depth files were inputted into [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/) [(Lee et al. 2019)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.1552) with the following command**

`python3 stargazer.py "genotype" -o "output_preffix" -d wgs -t "supported_gene" --vcf "input.vcf" -c "control_gene --gdf "input_rd.gdf"`

Examples for *CYP1A2* and *CYP2A6* can be found in `input/sge_output/`

**Stargazer**

Output Processing: [2021_Archaic_PGx_SGE_Output_Analysis.R](https://github.com/the-claw-lab/aDNA_PGx_2021/blob/main/2021_Archaic_PGx_SGE_Output_Analysis.R)

Example output found [here](https://github.com/the-claw-lab/aDNA_PGx_2021/tree/main/output/2021_Archaic_PGx_SGE_Output_Analysis_EXAMPLE_OUTPUT).

**SNV Analysis**

Select gene positions (start-stop gene loci + 2000bp upstream for the promoter) and filter the QUAL >= 40

`bcftools view -t "chr:strt-end" "input.vcf" > "output.vcf"`

`bcftools view -i '%QUAL>=40' "output.vcf" > "output_filtered40.vcf"`

SNVs are annotated with [Annovar](https://annovar.openbioinformatics.org/en/latest/), [CADD](https://cadd.gs.washington.edu), [SIFT](https://sift.bii.a-star.edu.sg), and [PolyPhen2](http://genetics.bwh.harvard.edu/pph2/)

SNV Analysis Script: [2021_Archaic_PGx_variantSNV_Analysis.R](https://github.com/the-claw-lab/aDNA_PGx_2021/blob/main/2021_Archaic_PGx_variantSNV_Analysis.R)

Example output found [here](https://github.com/the-claw-lab/aDNA_PGx_2021/tree/main/output/2021_Archaic_PGx_variantSNV_Analysis_EXAMPLE_OUTPUT).

**Potentially Damaging Variants**

Potentially Damaging Variant Script: [2021_Archaic_PGx_potentiallyDamagingSNVs_Analysis.R](https://github.com/the-claw-lab/aDNA_PGx_2021/blob/main/2021_Archaic_PGx_potentiallyDamagingSNVs_Analysis.R)

Example output found [here](https://github.com/the-claw-lab/aDNA_PGx_2021/tree/main/output/2021_Archaic_PGx_potentiallyDamagingSNVs_Analysis_EXAMPLE_OUTPUT).
