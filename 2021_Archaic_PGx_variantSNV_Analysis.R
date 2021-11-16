####Requires####
#Stargazer installed on computer
#Stargazer output for genes of interest
#Read depth for genes of interest
#Annovar table
#Genotype table

####MUST RUN STARGAZER ANALYSIS SCRIPT FIRST####

####Libraries####
suppressMessages(library(dplyr)) #Data Wrangling
suppressMessages(library(reshape2)) #Data Wrangling
suppressMessages(library(stringr)) #String Management
suppressMessages(library(ggplot2)) #Plotting
suppressMessages(library(RColorBrewer)) #Plotting 

#Gene List
geneList <- c("cyp1a2", "cyp2a6") #CAN BE "cyp1a2", "cyp2a6","cyp2b6","cyp2c8","cyp2c9","cyp2c19", "cyp2d6","cyp2e1", "cyp2j2", "cyp3a4", "cyp3a5"

mapQual <- 40
preffixPath <- "~/Desktop/2021_Archaic_PGx_Github/" #"~/path/to/working/dir/"
SGEPath = "~/Software/Stargazer/" #"~/path/to/Stargazer/"

#Create directory
if(!dir.exists(paste(preffixPath, "output/", sep = ""))) dir.create(paste(preffixPath, "output/", sep = ""))
path.output = paste(preffixPath, "output/2021_Archaic_PGx_variantSNV_Analysis_EXAMPLE_OUTPUT/", sep = "")
if(!dir.exists(path.output)) dir.create(path.output)

######VARIANT TABLE SET-UP#######
#SGE Calls
if(file.exists(paste(preffixPath, "2021_Archaic_PGX_Stargazer_modified_startable.csv", sep = ""))){
  sge_calls <- read.table(paste(preffixPath, "2021_Archaic_PGX_Stargazer_modified_startable.csv", sep = ""), header = T, sep = ",", fill = T)
}else{
  stop("Must download modified Stargazer Star Allele Table from Github and place in same WD: https://github.com/the-claw-lab/aDNA_PGx_2021/blob/main/2021_Archaic_PGX_Stargazer_modified_startable.csv")
}

#Clean up sge_calls table
sge_calls <- filter(sge_calls, sge_calls$type != "ref")
for(POS in unique(sge_calls$pos)){
  temp <- filter(sge_calls, sge_calls$pos == POS)
  gene <- temp$gene[1]
  pos <- POS
  type <- paste(temp$type, collapse = ",")
  name <- paste(temp$name, collapse = ",")
  mutation <- paste(temp$mutation, collapse = ",")
  
  out <- data.frame(gene, name, type, pos, mutation)
  
  if(!exists("sge_calls_out")){
    sge_calls_out <- out
  }else{
    sge_calls_out <- rbind(sge_calls_out, out)
  }
}

#Remove variants
rm(temp)
rm(out)

#Read in Starzazer Analysis
path.SGE = paste(preffixPath, "output/2021_Archaic_PGx_SGE_Output_Analysis_EXAMPLE_OUTPUT/out_table_variants_full.csv", sep = "")
if(file.exists(path.SGE)){
  sge_output <- read.table(path.SGE, sep = ",", header = TRUE)
}else{
  stop("PLEASE RUN STARGAZER ANALYSIS SCRIPT FIRST: 2021_Archaic_PGx_SGE_Output_Analysis.R")
}

#Read in annotated vcf, cadd scores, and read depth
for(gene in geneList){
  anno_path_preffix = paste(preffixPath, "input/genotype/", sep = "")
  #Read in annotated vcf
  anno_path <- paste(anno_path_preffix, "vcf_annovar/neanderthal_full_", gene, "_out_filtered", mapQual, "_region_annotated.hg19_multianno.vcf", sep = "")
  
  #Read in the columns first then data
  tmp_vcf<-readLines(anno_path)
  tmp_vcf_data<-read.table(anno_path, stringsAsFactors = FALSE)
  
  # filter for the columns names
  tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
  tmp_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
  names(tmp_vcf_data)<-tmp_names
  vcfIn <- tmp_vcf_data
  
  #Parse Annovar Annotations
  annotation_table <- subset(vcfIn, select = c(INFO))
  for(i in 1:nrow(annotation_table)){
    #avsnp150
    rsID <- str_match(annotation_table$INFO[i], "avsnp150=*(.*?)\\s*\\;")[,2]
    #RefGene
    Function.refGene <- str_match(annotation_table$INFO[i], "Func.refGeneWithVer=*(.*?)\\s*\\;")[,2]
    Gene.refGene <- str_match(annotation_table$INFO[i], "Gene.refGeneWithVer=*(.*?)\\s*\\;")[,2]
    Exonic_Function.refGene <- str_match(annotation_table$INFO[i], "ExonicFunc.refGeneWithVer=*(.*?)\\s*\\;")[,2]
    AminoAcid_Change.refGene <- str_match(annotation_table$INFO[i], "AAChange.refGeneWithVer=*(.*?)\\s*\\;")[,2]
    #KnownGene
    Function.UCSC <- str_match(annotation_table$INFO[i], "Func.knownGene=*(.*?)\\s*\\;")[,2]
    Gene.UCSC <- str_match(annotation_table$INFO[i], "Gene.knownGene=*(.*?)\\s*\\;")[,2]
    Exonic_Function.UCSC <- str_match(annotation_table$INFO[i], "ExonicFunc.knownGene=*(.*?)\\s*\\;")[,2]
    AminoAcid_Change.UCSC <- str_match(annotation_table$INFO[i], "AAChange.knownGene=*(.*?)\\s*\\;")[,2]
    #PolyPhen Score
    Polyphen_HVAR_Score <- str_match(annotation_table$INFO[i], "Polyphen2_HVAR_score=*(.*?)\\s*\\;")[,2]
    Polyphen2_HVAR_Prediction <- str_match(annotation_table$INFO[i], "Polyphen2_HVAR_pred=*(.*?)\\s*\\;")[,2]
    #Alt Allele Frequency
    Minor_Allele_Frequency <- str_match(annotation_table$INFO[i], "AF=*(.*?)\\s*\\;")[,2]
    #Annovar Date
    ANNOVAR.Date <- str_match(annotation_table$INFO[i], "ANNOVAR_DATE=*(.*?)\\s*\\;")[,2]
    
    out <- data.frame(rsID, 
                      Function.refGene, Gene.refGene, Exonic_Function.refGene, AminoAcid_Change.refGene,
                      Function.UCSC, Gene.UCSC, Exonic_Function.UCSC, AminoAcid_Change.UCSC,
                      Polyphen_HVAR_Score, Polyphen2_HVAR_Prediction,
                      Minor_Allele_Frequency,
                      ANNOVAR.Date)
    
    if(!exists("annot_out")){
      annot_out <- out
    }else{
      annot_out <- rbind(annot_out, out)
    }
  }
  
  vcfIn <- subset(vcfIn, select = -c(INFO, ID, FORMAT))
  colnames(vcfIn) <- c("CHR", "POS", "REF", "ALT", "QUAL", "FILTER", 
                       "Altai_Neanderthal.gt", "Chagyrskaya_8_Neanderthal.gt", "Denisovan.gt", "Vindija_33.19_Neanderthal.gt")
  vcf_In_AD <- subset(vcfIn, select = c(Altai_Neanderthal.gt, Chagyrskaya_8_Neanderthal.gt, Denisovan.gt, Vindija_33.19_Neanderthal.gt))
  colnames(vcf_In_AD) <- c("Altai_Neanderthal.AD", "Chagyrskaya_8_Neanderthal.AD", "Denisovan.AD", "Vindija_33.19_Neanderthal.AD")
  
  vcfIn <- cbind(vcfIn, vcf_In_AD, annot_out)
  rm(annot_out)
  
  #Add CADD Scores
  cadd_path <- paste(anno_path_preffix, "CADD_annotation/GRCh37-v1.6_", gene, ".tsv", sep = "")
  cadd_in <- read.table(cadd_path, header = FALSE, sep = "\t")
  colnames(cadd_in) <- c("CHR", "POS", "REF", "ALT", "CADD.Raw_Score", "CADD.Phred_Normalized")
  cadd_in <- subset(cadd_in, select = c(POS, CADD.Raw_Score, CADD.Phred_Normalized))
  
  #Add SIFT Scores
  sift_path <- paste(anno_path_preffix, "SIFT_annotation/neanderthal_full_", gene, "_out_filtered", mapQual, "_region_annotated.hg19_multianno_SIFTannotations.xls", sep = "")
  sift_in <- read.table(sift_path, header = TRUE, sep = "\t")
  colnames(sift_in) <- c("CHROM", "POS", "REF", "ALT", "TRANSCRIPT_ID",
                         "GENE_ID",  "GENE_NAME", "REGION", "VARIANT_TYPE",
                         "REF_AMINO", "ALT_AMINO", "AMINO_POS", "SIFT_Score",
                         "SIFT_Median", "NUM_SEQS", "dbSNP", "SIFT_prediction")
  sift_in <- subset(sift_in, select = c(POS, SIFT_Score, SIFT_Median, SIFT_prediction))
  sift_in$SIFT_prediction <- gsub("TOLERATED", "T", sift_in$SIFT_prediction)
  sift_in$SIFT_prediction <- gsub("DELETERIOUS", "D", sift_in$SIFT_prediction)
  
  #Add read depth
  CHR <- as.character(vcfIn$CHR[1])
  rd_path <- paste(preffixPath, "input/depth_of_coverage/neanderthal_vdr_", gene, "_out.gdf", sep = "")
  rd_in <- read.table(rd_path, header = TRUE, sep = "\t")
  temp <- data.frame(do.call('rbind', strsplit(as.character(rd_in$Locus),':',fixed=TRUE)))
  rd_in <- cbind(temp, rd_in)
  rd_in <- filter(rd_in, rd_in$X1 == CHR)
  rd_in <- subset(rd_in, select = c(X2, Total_Depth, Average_Depth_sample, Depth_for_AltaiNea, Depth_for_ChagyrskayaPhalanx, 
                                    Depth_for_DNK02, Depth_for_Vi33.19))
  colnames(rd_in) <- c("POS", "Total_Depth", "Average_Depth", "Altai_Neanderthal.dp", "Chagyrskaya_8_Neanderthal.dp", "Denisovan.dp", "Vindija_33.19_Neanderthal.dp")
  
  #Add in Stargazer star alleles
  #SGE Star Table
  sge_call_byGene <- filter(sge_calls_out,sge_calls_out$gene == tolower(gene))
  sge_call_byGene <- subset(sge_call_byGene, select = c(pos, name, mutation, type))
  sge_call_byGene <- filter(sge_call_byGene, sge_call_byGene$type != "ref")
  colnames(sge_call_byGene) <- c("POS", "SGE.name", "SGE.mutation", "SGE.variant")
  
  #Merge all files
  out <- merge(x = vcfIn, y = cadd_in, by = "POS", all.x = TRUE)
  out <- merge(x = out, y = sift_in, by = "POS", all.x = TRUE)
  out <- merge(x = out, y = sge_call_byGene, by = "POS", all.x = TRUE)
  out <- merge(x = out, y = rd_in, by = "POS", all.x = TRUE)
  
  out$GENE_CALL <- toupper(gene)
  
  if(!exists("out_table")){
    out_table <- out
  }else{
    out_table <- rbind(out_table, out)
  }
}

#Clean up genotype columns and add in allelic depths
for(i in 1:nrow(out_table)){
  ref <- out_table$REF[i]
  alt <- unlist(strsplit(out_table$ALT[i], ","))
  for(sample in c("Altai_Neanderthal","Chagyrskaya_8_Neanderthal","Denisovan","Vindija_33.19_Neanderthal")){
    #First gt
    colname.gt <- paste(sample, ".gt", sep = "")
    genotypes <- unlist(strsplit(out_table[i,colname.gt], ":"))[1]
    genotypes <- gsub("\\|", "/", genotypes)
    genotypes <- gsub("0", ref, genotypes)
    genotypes <- gsub("1", alt[1], genotypes)
    genotypes <- gsub("2", alt[2], genotypes)
    genotypes <- gsub("3", alt[3], genotypes)
    genotypes <- gsub("4", alt[4], genotypes)
    genotypes <- gsub("5", alt[5], genotypes)
    out_table[i,colname.gt] <- genotypes
    
    #Allelic Depth
    colname.AD <- paste(sample, ".AD", sep = "")
    AD <- unlist(strsplit(out_table[i,colname.AD], ":"))[2]
    out_table[i,colname.AD] <- AD
  }
}

#Check for dupliate positions due to two different CADD scores and choose the higher cadd score (manually checked and none of these are over 10)
dup_pos <- out_table[duplicated(out_table[,1]),]$POS
dup_table <- filter(out_table, out_table$POS %in% dup_pos)
no_dup_table <- out_table %>% dplyr::filter(!out_table$POS %in% dup_pos)

splitByPOS <- split(dup_table, dup_table$POS)
for(i in splitByPOS){
  CADD.Raw_Score <- max(i$CADD.Raw_Score)
  CADD.Phred_Normalized <- max(i$CADD.Phred_Normalized)
  out <- i[1,]
  out$CADD.Raw_Score[1] <- CADD.Raw_Score
  out$CADD.Phred_Normalized[1] <- CADD.Phred_Normalized
  
  if(!exists("dedup_table")){
    dedup_table <- out
  }else{
    dedup_table <- rbind(dedup_table, out)
  }
}

out_table <- rbind(no_dup_table, dedup_table)

#Fix rsid names
for(i in 1:nrow(out_table)){
  if(out_table$rsID[i] == "."){
    out_table$rsID[i] <- paste("pos", out_table$POS[i], sep = "")
  }
}

#Remove all unnecssary tables
rm(out)
rm(rd_in)
rm(sge_call_byGene)
rm(sge_calls_out)
rm(sge_calls)
rm(temp)
rm(dup_table)
rm(dedup_table)
rm(no_dup_table)

out_table[is.na(out_table)]="."

out_table$Gene.refGene <- gsub("P\\\\x3b", ",", out_table$Gene.refGene)
out_table$Gene.refGene <- gsub("\\\\x3b", ",", out_table$Gene.refGene)
out_table$Gene.UCSC <- gsub("P\\\\x3b", ",", out_table$Gene.UCSC)
out_table$Gene.UCSC <- gsub("\\\\x3b", ",", out_table$Gene.UCSC)

#Make ANNOVAR ordering table
annovar_order <- data.frame(c("exonic", "splicing", "ncRNA_exonic", "ncRNA_intronic",
                              "UTR5", "UTR3", "intronic", "upstream", "downstream", "intergenic"),
                            c(1, 1, 2, 2, 3, 3, 4, 5, 5, 6))
colnames(annovar_order) <- c("location", "score")

#Check if genes and annotations agree
out_table$Function.finalCall <- NA
out_table$check.gene <- NA
out_table$check.annot <- NA
for(i in 1:nrow(out_table)){
  #determine if genes are the same
  gene <- out_table$GENE_CALL[i]
  gene_UCSC <- grepl(gene, out_table$Gene.UCSC[i])
  gene_refGene <- grepl(gene, out_table$Gene.refGene[i])
  
  if(gene_refGene & gene_UCSC){
    out_table$check.gene[i] <- TRUE
  }else{
    out_table$check.gene[i] <- paste("Gene:", gene, ",",
                                     "refGene:",out_table$Gene.refGene[i], ",",
                                     "UCSC:", out_table$Gene.UCSC[i], sep = "")
  }
  
  #Determine if the calls are the same
  if(out_table$Function.refGene[i] == out_table$Function.UCSC[i]){
    out_table$check.annot[i] <- TRUE
    out_table$Function.finalCall[i] <- out_table$Function.refGene[i]
  }else{
    out_table$check.annot[i] <- paste("refGene:",out_table$Function.refGene[i], ",",
                                      "UCSC:", out_table$Function.UCSC[i], sep = "")
    #Check to see if genes agree or not
    if(!isTRUE(out_table$check.gene[i])){#First find the easy ones that just dont have the same gene
      if(gene_refGene){
        out_table$Function.finalCall[i] <- out_table$Function.refGene[i]
      }else if(gene_UCSC){
        out_table$Function.finalCall[i] <- out_table$Function.UCSC[i]
      }
    }
    #Now look for ones that are still blank
    if(is.na(out_table$Function.finalCall[i])){#Now the ordering from ANNOVAR
      var_locations <- append(out_table$Function.refGene[i], out_table$Function.UCSC[i])
      var_compare <- filter(annovar_order, annovar_order$location %in% var_locations)
      min_score <- min(var_compare$score)
      var_compare <- filter(var_compare, var_compare$score == min_score)
      out_table$Function.finalCall[i] <- var_compare$location[1]
    }
  }
}

#Split by gene
splitByGene <- split(out_table, out_table$GENE_CALL)

#Now make a list of variants for each sample for each gene
for(i in splitByGene){
  GENE <- i$GENE_CALL[1]
  
  #Make List of positions that are the alternative for each sample
  for(sample_name in c("Altai_Neanderthal","Chagyrskaya_8_Neanderthal","Denisovan","Vindija_33.19_Neanderthal")){
    #Make Sample Table
    sample_table <- as.data.frame(i[,paste(sample_name, ".gt", sep = "")])
    colnames(sample_table) <- "Genotype"
    
    ad_table <- as.data.frame(i[,paste(sample_name, ".AD", sep = "")])
    colnames(ad_table) <- "Genotype_AllelicDepth"
    
    ALT_table <- subset(i, select = c(ALT))
    sample_table <- cbind(sample_table, ALT_table)
    sample_table$Genotype_Numeric <- NA
    
    #Remove uncalled variants
    sample_table$Genotype <- gsub("\\.\\/\\.", NA, sample_table$Genotype)
    
    for(k in 1:nrow(sample_table)){
      geno <- sample_table$Genotype[k]
      if(!is.na(geno)){
        geno <- unlist(str_split(geno, "/"))
        alt <- unlist(str_split(sample_table$ALT[k], ","))
        temp <- 0
        for(l in alt){
          temp <- temp + as.numeric(sum(geno==l))
        }
      }else{
        temp <- NA
      }
      sample_table$Genotype_Numeric[k] <- temp
    }
    
    sample_table <- subset(sample_table, select = -c(ALT))
    
    rd_table <- as.data.frame(i[,paste(sample_name, ".dp", sep = "")])
    colnames(rd_table) <- "Read_Depth"
    
    info_table <- subset(i, select = c(CHR, POS, REF, ALT, rsID, QUAL, 
                                       Function.refGene, Exonic_Function.refGene, AminoAcid_Change.refGene,
                                       Function.UCSC, Exonic_Function.UCSC, AminoAcid_Change.UCSC,
                                       Function.finalCall, Minor_Allele_Frequency,
                                       SGE.name, SGE.mutation, SGE.variant,
                                       CADD.Raw_Score, CADD.Phred_Normalized,
                                       SIFT_Score, SIFT_prediction, SIFT_Median, Polyphen_HVAR_Score, Polyphen2_HVAR_Prediction))
    
    info_table$SGE.alt_var <- gsub(".>", "", info_table$SGE.mutation)
    
    temp_table <- cbind(sample_table, ad_table, rd_table,info_table)
    
    temp_table$Sample <- sample_name
    temp_table$Gene <- GENE
    
    if(!exists("variants_oi_table")){
      variants_oi_table <- temp_table
    }else{
      variants_oi_table <- rbind(variants_oi_table, temp_table)
    }
  }
}

#Add in stargazer final call
variants_oi_table$Diplotype <- NA
variants_oi_table$Phenotype <- NA
for(i in 1:nrow(variants_oi_table)){
  sample <- variants_oi_table$Sample[i]
  gene <- variants_oi_table$Gene[i]
  diplo_in <- filter(sge_output, sge_output$Gene == gene)
  diplo_in <- filter(diplo_in, diplo_in$Sample == sample)
  diplo <- diplo_in$Diplotype[1]
  pheno <- diplo_in$Phenotype[1]
  variants_oi_table$Diplotype[i] <- diplo
  variants_oi_table$Phenotype[i] <- pheno
}

#Reorder data frame for output
variants_oi_table <- variants_oi_table[c("Gene", "Sample", "Diplotype", "Phenotype", "CHR", "POS", "REF", "ALT", "Genotype", "Genotype_Numeric","Genotype_AllelicDepth",
                                         #"Genotype_MaxPlank", 
                                         "rsID", "Read_Depth", "QUAL", "Function.finalCall", "Minor_Allele_Frequency",
                                         "Function.refGene", "Exonic_Function.refGene","AminoAcid_Change.refGene", 
                                         "Function.UCSC", "Exonic_Function.UCSC","AminoAcid_Change.UCSC", 
                                         "SGE.name", "SGE.mutation", "SGE.variant", "SGE.alt_var",
                                         "CADD.Raw_Score", "CADD.Phred_Normalized", "SIFT_Score", "SIFT_prediction", "SIFT_Median", "Polyphen_HVAR_Score",
                                         "Polyphen2_HVAR_Prediction")]



variants_oi_table$Sample <- gsub("_", " ", variants_oi_table$Sample)
variants_oi_table$Function.refGene <- toupper(variants_oi_table$Function.refGene)
variants_oi_table$Function.UCSC <- toupper(variants_oi_table$Function.UCSC)
variants_oi_table$Function.finalCall <- toupper(variants_oi_table$Function.finalCall)
variants_oi_table <- unique(variants_oi_table)

#####Allele balance#######
allelicBiasRatio <- 0.2
variants_oi_table$AlleleBalance <- NA

for(i in 1:nrow(variants_oi_table)){
  row_rem <- FALSE
  if(!is.na(variants_oi_table$Genotype[i])){
    alleles <- append(variants_oi_table$REF[i], unlist(strsplit(variants_oi_table$ALT[i], ","))) #Reference is always the first position
    genotypes <- unlist(strsplit(variants_oi_table$Genotype[i], "/"))
    #Set allelic depths
    allelic_depth <- as.numeric(unlist(strsplit(variants_oi_table$Genotype_AllelicDepth[i], ",")))
    
    #Determine zygosity composition
    heterozygous = ifelse(genotypes[1] == genotypes[2], FALSE, TRUE)
    
    #Homozygous case
    if(!heterozygous){
      index_1 = match(genotypes[1], alleles)
      ad_1 = allelic_depth[index_1]
      ad_2 = min(as.numeric(allelic_depth))
    }else{
      index_1 = match(genotypes[1], alleles)
      ad_1 = allelic_depth[index_1]
      index_2 = match(genotypes[2], alleles)
      ad_2 = allelic_depth[index_2]
    }
    
    total_depth = sum(as.numeric(allelic_depth))
    ad_min = min(ad_1, ad_2)
    
    #Find the ratio
    ratio = ifelse(ad_min == 0, 1, (ad_min/total_depth))
    variants_oi_table$AlleleBalance[i] <- ratio
    
    #Now look for heterozyous sites which have an allele depth value bias
    if(genotypes[1] != genotypes[2]){
      #First, remove heterozygous sites with a single read for one allele
      if(ad_1 == 1 | ad_2 == 1){
        removed_row <- variants_oi_table[i,]
        variants_oi_table$Genotype[i] <- NA
        variants_oi_table$Genotype_Numeric[i] <- NA
        row_rem <- TRUE
      }else if(ratio < allelicBiasRatio | ratio > 0.5){
        removed_row <- variants_oi_table[i,]
        variants_oi_table$Genotype[i] <- NA
        variants_oi_table$Genotype_Numeric[i] <- NA
        row_rem <- TRUE
      }
      
      if(row_rem){
        #Make table of removed variants
        if(!exists("removed_vars_table")){
          removed_vars_table <- removed_row
        }else{
          removed_vars_table <- rbind(removed_vars_table, removed_row)
        }
      }
    }
  }
}

#Make Histogram of allele balance for heterozygous sites
DF.ab.full = rbind(variants_oi_table, removed_vars_table) %>%
  filter(Genotype_Numeric == 1) %>%
  dplyr::select(Gene, Genotype, Genotype_Numeric, AlleleBalance, Genotype_AllelicDepth, REF, ALT)
DF.ab.filtered = variants_oi_table %>%
  filter(Genotype_Numeric == 1) %>%
  dplyr::select(Gene, Genotype, Genotype_Numeric, AlleleBalance, Genotype_AllelicDepth, REF, ALT)

#Plot
GG.ab.full = ggplot(DF.ab.full, aes(x=AlleleBalance)) + 
  geom_histogram(binwidth = 0.01, fill = "blue") +
  geom_vline(xintercept = allelicBiasRatio, color = "grey70") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,250) +
  labs(x = "Allele Balance", y = "Count") +
  theme(
    axis.text = element_text(color = "black", face= "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank()
  )
GG.ab.filtered = ggplot(DF.ab.filtered, aes(x=AlleleBalance)) + 
  geom_histogram(binwidth = 0.01, fill = "blue") +
  geom_vline(xintercept = allelicBiasRatio, color = "grey70") +
  theme_bw() +
  xlim(0,1) +
  ylim(0,250) +
  labs(x = "Allele Balance", y = "Count") +
  theme(
    axis.text = element_text(color = "black", face= "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.major.x = element_blank()
  )

#Save plots
ab.full.Path <- paste(path.output, "nea_whole_snp_abDistribution_full.png", sep = "")
ab.filtered.Path <- paste(path.output, "nea_whole_snp_abDistribution_filtered.png", sep = "")

suppressWarnings(ggsave(ab.full.Path, GG.ab.full, width = 8, height = 5, units = "in", dpi = 400, bg = "white"))
suppressWarnings(ggsave(ab.filtered.Path, GG.ab.filtered, width = 8, height = 5, units = "in", dpi = 400, bg = "white"))

#Rm unnecissary points
rm(ALT_table)
rm(diplo_in)
rm(info_table)
rm(rd_table)
rm(sample_table)
rm(temp_table)

#######MAKE BASIC STATS TABLE AND CLEAN Variant oi Table###########
basic_stats_in <- variants_oi_table
basic_stats_in$Genotype_Numeric <- gsub(0, NA, basic_stats_in$Genotype_Numeric)

#Change names to make less categories
basic_stats_in$Function.finalCall <- gsub("UTR5", "UTR", basic_stats_in$Function.finalCall)
basic_stats_in$Function.finalCall <- gsub("UTR3", "UTR", basic_stats_in$Function.finalCall)
basic_stats_in$Function.finalCall <- gsub("NCRNA_EXONIC", "NCRNA", basic_stats_in$Function.finalCall)
basic_stats_in$Function.finalCall <- gsub("NCRNA_INTRONIC", "NCRNA", basic_stats_in$Function.finalCall)
basic_stats_in$Function.finalCall <- gsub("UPSTREAM", "PROMOTER", basic_stats_in$Function.finalCall)
basic_stats_in$Function.finalCall <- gsub("INTERGENIC", "PROMOTER", basic_stats_in$Function.finalCall)

splitByGene <- split(basic_stats_in, basic_stats_in$Gene)
for(i in splitByGene){
  Gene <- i$Gene[1]
  
  info_table <- unique(subset(i, select = c(POS, Function.finalCall)))
  gt <- subset(i, select = c(POS, Sample, Genotype_Numeric))
  gt <- gt %>% dcast(POS ~ Sample, value.var = "Genotype_Numeric")
  
  for(j in 1:nrow(gt)){
    if(is.na(gt$`Altai Neanderthal`[j]) & is.na(gt$`Chagyrskaya 8 Neanderthal`[j]) & 
       is.na(gt$Denisovan[j]) & is.na(gt$`Vindija 33.19 Neanderthal`[j])){
    }else{
      #Add this row to the cleaned table
      if(!exists("gt_out")){
        gt_out <- gt[j,]
      }else{
        gt_out <- rbind(gt_out, gt[j,])
      }
    }
  }
  
  #Now do stats on this cleaned table
  gt_stat <- merge(x = gt_out, y = info_table, by = "POS", all.x = TRUE)
  rm(gt_out)
  
  #Filter variant oi table by list of positions
  variant_list <- as.data.frame(gt_stat$POS)
  variant_list$Gene <- Gene
  colnames(variant_list) <- c("POS", "Gene")
  
  if(!exists("variant_list_out")){
    variant_list_out <- variant_list
  }else{
    variant_list_out <- rbind(variant_list_out, variant_list)
  }
  
  #Split by Population
  for(Sample in c("Denisovan", "Altai Neanderthal", "Chagyrskaya 8 Neanderthal", "Vindija 33.19 Neanderthal")){
    population_in <- as.data.frame(gt_stat[,Sample])
    colnames(population_in) <- "Sample"
    info_in <- subset(gt_stat, select = c(POS, Function.finalCall))
    in_table <- cbind(population_in, info_in)
    in_table <- na.omit(in_table)
    
    Function <- "TOTAL"
    Count <- as.numeric(nrow(in_table))
    out <- data.frame(Sample, Gene, Function, Count)
    if(!exists("gene_wide_summary_table_bySample")){
      gene_wide_summary_table_bySample <- out
    }else{
      gene_wide_summary_table_bySample <- rbind(gene_wide_summary_table_bySample, out)
    }
    
    #By Function full
    splitByFunction <- split(in_table, in_table$Function.finalCall)
    for(j in splitByFunction){
      Function <- j$Function.finalCall[1]
      Count <- as.numeric(nrow(j))
      out <- data.frame(Sample, Gene, Function, Count)
      
      if(!exists("gene_wide_summary_table_bySample")){
        gene_wide_summary_table_bySample <- out
      }else{
        gene_wide_summary_table_bySample <- rbind(gene_wide_summary_table_bySample, out)
      }
    }
  }
  
  #Now make stats table full
  Function <- "TOTAL"
  Count <- as.numeric(nrow(gt_stat))
  out <- data.frame(Gene, Function, Count)
  if(!exists("gene_wide_summary_table")){
    gene_wide_summary_table <- out
  }else{
    gene_wide_summary_table <- rbind(gene_wide_summary_table, out)
  }
  
  #By Function full
  splitByFunction <- split(gt_stat, gt_stat$Function.finalCall)
  for(j in splitByFunction){
    Function <- j$Function.finalCall[1]
    Count <- as.numeric(nrow(j))
    out <- data.frame(Gene, Function, Count)
    
    if(!exists("gene_wide_summary_table")){
      gene_wide_summary_table <- out
    }else{
      gene_wide_summary_table <- rbind(gene_wide_summary_table, out)
    }
  }
}

#Clean up gene wide summary table
gene_wide_summary_table$Function <- str_to_title(gene_wide_summary_table$Function)
gene_wide_summary_table_cast <- gene_wide_summary_table %>% dcast(Function ~ Gene, value.var = "Count")
gene_wide_summary_table_cast[is.na(gene_wide_summary_table_cast)] = 0
gene_wide_summary_table_cast$Function <- str_to_title(gene_wide_summary_table_cast$Function)

#Make a quick summary table
#First find total 
total <- filter(gene_wide_summary_table, gene_wide_summary_table$Function == "Total")
total <- sum(total$Count)
out_summary <- data.frame(Desc=c("Total"), Count=c(total), Prcent=c(100))

#Find counts for each subset
sub_func <- filter(gene_wide_summary_table, gene_wide_summary_table$Function != "Total")
splitByFunc <- split(sub_func, sub_func$Function)
for(i in splitByFunc){
  Desc <- paste(i$Function[1], "_Total", sep = "")
  Count <- sum(i$Count)
  Prcent <- (Count/total)*100
  
  out <- data.frame(Desc, Count, Prcent)
  out_summary <- rbind(out_summary, out)
}

splitByGene <- split(gene_wide_summary_table, gene_wide_summary_table$Gene)
for(i in splitByGene){
  gene <- i$Gene[1]
  #First find total
  t_sub <- filter(i, i$Function == "Total")
  t_sub <- t_sub$Count[1]
  Count <- t_sub
  Desc <- paste(gene, "Total", sep = "_")
  Prcent <- (Count/total)*100
  out <- data.frame(Desc, Count, Prcent)
  out_summary <- rbind(out_summary, out)
  
  #Now look at each location
  temp <- filter(i, i$Function != "Total")
  for(n in 1:nrow(temp)){
    Desc <- paste(gene, temp$Function[n], sep = "_")
    Count <- temp$Count[n]
    Prcent <- (Count/t_sub)*100
    
    out <- data.frame(Desc, Count, Prcent)
    out_summary <- rbind(out_summary, out)
  }
}

##Add in which sample has the most variants
#Find counts for each sample
sub_sample <- filter(gene_wide_summary_table_bySample, gene_wide_summary_table_bySample$Function == "TOTAL")
splitBySample <- split(sub_sample, sub_sample$Sample)
for(i in splitBySample){
  Desc <- paste(i$Sample[1], "_Total", sep = "")
  Count <- sum(i$Count)
  Prcent <- NA
  
  out <- data.frame(Desc, Count, Prcent)
  out_summary <- rbind(out_summary, out)
}

#Now look at which sample has most variants for each sample
splitByGene <- split(sub_sample, sub_sample$Gene)
for(i in splitByGene){
  Desc <- paste(i$Gene[1], "_MostVars", sep = "")
  max_var <- max(i$Count)
  Count <- filter(i, i$Count == max_var)$Count[1]
  Prcent <- filter(i, i$Count == max_var)$Sample[1]
  
  out <- data.frame(Desc, Count, Prcent)
  out_summary <- rbind(out_summary, out)
  
  #Now add all stats
  temp <- i
  temp$Desc <- paste(temp$Gene, "_", temp$Sample, sep = "")
  temp <- subset(temp, select = c(Desc, Count))
  temp$Prcent <- NA
  
  out_summary <- rbind(out_summary, temp)
}

######Clean up variants oi table#####
splitByGene <- split(variants_oi_table, variants_oi_table$Gene)
for(i in splitByGene){
  gene <- i$Gene[1]
  
  filter_list <- filter(variant_list_out, variant_list_out$Gene == gene)
  filter_list <- filter_list$POS
  
  out <- filter(i, i$POS %in% filter_list)
  
  if(!exists("variants_oi_table_out")){
    variants_oi_table_out <- out
  }else{
    variants_oi_table_out <- rbind(variants_oi_table_out, out)
  }
}

#Save Spreadsheets to csv
path_full <- paste(path.output, "nea_whole_snp_overview.csv", sep = "")
path_cast <- paste(path.output, "nea_whole_snp_cast.csv", sep = "")
path_sample <- paste(path.output, "nea_whole_snp_bySample.csv", sep = "")
path_fvt <- paste(path.output, "full_variant_table_alleleBalanceRatio_" , allelicBiasRatio, ".csv", sep = "")
path_removed <- paste(path.output, "full_variant_table_alleleBalanceRatio_" , allelicBiasRatio, "_REMOVED_VARS.csv", sep = "")
path_summary <- paste(path.output, "full_variant_summary.csv", sep = "")

write.table(variants_oi_table_out, path_fvt, sep = ",", row.names = FALSE, col.names = TRUE)
write.table(removed_vars_table, path_removed, sep = ",", row.names = FALSE, col.names = TRUE)
write.table(gene_wide_summary_table, path_full, sep = ",", row.names = FALSE, col.names = TRUE)
write.table(gene_wide_summary_table_cast, path_cast, sep = ",", row.names = FALSE, col.names = TRUE)
write.table(gene_wide_summary_table_bySample, path_sample, sep = ",", row.names = FALSE, col.names = TRUE)
write.table(out_summary, path_summary, sep = ",", row.names = FALSE, col.names = TRUE)

#######Now make bar plot for frequencies######
#First clean table
variant_summary_toPlot <- filter(gene_wide_summary_table_bySample, gene_wide_summary_table_bySample$Function != "TOTAL")
variant_summary_toPlot <- filter(variant_summary_toPlot, variant_summary_toPlot$Function != "NCRNA")
variant_summary_toPlot$Function <- str_to_title(variant_summary_toPlot$Function)
variant_summary_toPlot$Function <- gsub("Utr", "UTR", variant_summary_toPlot$Function)
variant_summary_toPlot$Sample <- gsub(" Neanderthal", "", variant_summary_toPlot$Sample)
variant_summary_toPlot$Sample <- gsub(" 8", "", variant_summary_toPlot$Sample)
variant_summary_toPlot$Sample <- gsub(" 33.19", "", variant_summary_toPlot$Sample)

#Now factor
variant_summary_toPlot <- variant_summary_toPlot %>%
  mutate(Sample = factor(Sample, levels=c("Denisovan", "Altai", "Chagyrskaya", "Vindija")))

variant_summary_toPlot$Gene_f = factor(variant_summary_toPlot$Gene, levels=c('CYP1A2','CYP2A6', 'CYP2B6', 'CYP2C8', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP2E1', 'CYP2J2', 'CYP3A4', 'CYP3A5'))
variant_summary_toPlot$Function = factor(variant_summary_toPlot$Function, levels=c('Exonic','Splicing', 'UTR', 'Intronic', 'Promoter'))

#Make a nice color palette
colPalette <- brewer.pal(n = 5, name = "Set1")

p <- ggplot(variant_summary_toPlot, aes(fill=Function, y=Count, x=as.factor(Sample))) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Gene_f) +
  scale_fill_manual(values = colPalette) +
  theme_bw()+
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "") +
  theme(plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 20, face = "bold", hjust = 1),
        axis.text.x = element_text(color = "black", size = 15, face = "bold", hjust = 1, vjust = 0.5, angle = 90),
        axis.title = element_text(size = 12, face = "bold"),
        legend.key.size = unit(2, "line"), #adjusts legend size
        legend.text = element_text(size = 15, face = "bold"),
        legend.title.align = 0.2,
        legend.title = element_text(size = 12, face = "bold", hjust = 0, vjust = 0.5),
        rect = element_rect(fill = "transparent"), # all rectangles
        panel.background = element_rect(fill = "transparent"), 
        panel.border = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 15, face = "bold.italic", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1.2, "line"))

plotPath <- paste(path.output, "nea_whole_snp_overview.png", sep = "")
suppressWarnings(ggsave(plotPath, p, width = 15, height = 7, units = "in", dpi = 400, bg = "transparent"))


