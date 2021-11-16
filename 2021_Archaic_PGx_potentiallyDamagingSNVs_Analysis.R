####Requires####
#Stargazer installed on computer
#Stargazer output for genes of interest
#Read depth for genes of interest
#Annovar table
#Genotype table

####MUST RUN SNV ANALYSIS SCRIPT FIRST####

####Libraries####
suppressMessages(library(dplyr)) #Data Wrangling
suppressMessages(library(reshape2)) #Data Wrangling
suppressMessages(library(stringr)) #String Management
suppressMessages(library(ggplot2)) #Plotting

#####Parameters######
ABR <- "0.2" #Allele Balance Ratio
CADDCutoff <- 20
BOOL.save = T

####File Paths####
preffix <- "~/Desktop/2021_Archaic_PGx_Github/" #"~/path/to/working/dir/"

#Save path
if(BOOL.save){
  path.save.preffix = paste(preffix, "output/2021_Archaic_PGx_potentiallyDamagingSNVs_Analysis_EXAMPLE_OUTPUT/", sep = "")
  if(!dir.exists(path.save.preffix)) dir.create(path.save.preffix)
}

#Import files
path_fvt <- paste(preffix, "output/2021_Archaic_PGx_variantSNV_Analysis_EXAMPLE_OUTPUT/full_variant_table_alleleBalanceRatio_" , ABR, ".csv", sep = "")
path_removed <- paste(preffix, "output/2021_Archaic_PGx_variantSNV_Analysis_EXAMPLE_OUTPUT/full_variant_table_alleleBalanceRatio_" , ABR, "_REMOVED_VARS.csv", sep = "")
path_poten_novel <- paste(preffix, "2021_Archaic_PGX_Gnomad-dbSNP_manualAnnotation.csv", sep = "")

#Import files from SNV script
if(file.exists(path_fvt)){
  variants_oi_table <- read.table(path_fvt, header = TRUE, sep = ",")
  removed_variants_table <- read.table(path_removed, header = TRUE, sep = ",")
}else{
  stop("PLEASE RUN SNV ANALYSIS SCRIPT FIRST: 2021_Archaic_PGx_variantSNV_Analysis.R")
}

#Import Manual Annotation File
if(file.exists(path_poten_novel)){
  manual_annotation <- read.table(path_poten_novel, header = T, sep = ",", fill = T)
}else{
  stop("Must download Gnomad / dbSNP Manual Annotation file from Github and place in same WD: https://github.com/the-claw-lab/aDNA_PGx_2021/blob/main/2021_Archaic_PGX_Gnomad-dbSNP_manualAnnotation.csv")
}

manual_annotation <- unique(subset(manual_annotation, select = c(CHR, POS, RSID))) %>% filter(RSID != ".")

#Add in other variants from Gnomad v2.1.1
variants_oi_table_manual <- filter(variants_oi_table, variants_oi_table$POS %in% manual_annotation$POS)
variants_oi_table_no_man <- variants_oi_table %>% dplyr::filter(!POS %in% manual_annotation$POS)

for(i in 1:nrow(variants_oi_table_manual)){
  pos <- as.character(variants_oi_table_manual$POS[i])
  chr <- as.character(variants_oi_table_manual$CHR[i])

  newID <- filter(manual_annotation, manual_annotation$CHR == chr & manual_annotation$POS == pos)$RSID[1]

  if(newID != "."){
    variants_oi_table_manual$rsID[i] <- newID
  }
}

variants_oi_table <- rbind(variants_oi_table_manual, variants_oi_table_no_man)

#######Look for variants of Interest visually########
#Find variants associated with star alleles
variants_oi_table$is.SGE <- "."
for(i in 1:nrow(variants_oi_table)){
  if(variants_oi_table$SGE.alt_var[i] != "."){
    variants_oi_table$is.SGE[i] <- "yes"
  }
}
#Convert CADD Score to numeric
variants_oi_table$CADD.Phred_Normalized <- as.numeric(variants_oi_table$CADD.Phred_Normalized)

#Add in conditions based on CADD, SIFT, and polyphen
variants_oi_table$Conditions <- NA
for(i in 1:nrow(variants_oi_table)){
  temp <- NA
  # if(variants_oi_table$is.SGE[i] == "yes"){
  #   variants_oi_table$rsID[i] <-  paste(variants_oi_table$rsID[i], "*", sep = "")
  # }
  if(variants_oi_table$CADD.Phred_Normalized[i] >= CADDCutoff){
    if(!is.na(temp)){
      temp <- paste(temp, "C", sep = "|")
    }else{
      temp <- "C"
    }
  }
  if(variants_oi_table$SIFT_prediction[i] == "D"){
    if(!is.na(temp)){
      temp <- paste(temp, "S", sep = "|")
    }else{
      temp <- "S"
    }
  }
  if(variants_oi_table$Polyphen2_HVAR_Prediction[i] == "D" | variants_oi_table$Polyphen2_HVAR_Prediction[i] == "P"){
    if(!is.na(temp)){
      temp <- paste(temp, "P", sep = "|")
    }else{
      temp <- "P"
    }
  }
  
  variants_oi_table$Conditions[i] <- temp
}

#Make variant annotations for plotting
for(i in 1:nrow(variants_oi_table)){
  if(!is.na(variants_oi_table$Conditions[i])){
    variants_oi_table$rsID_annot[i] <- paste(variants_oi_table$rsID[i], " (", variants_oi_table$Conditions[i], ")", sep = "")
  }else{
    variants_oi_table$rsID_annot[i] <- variants_oi_table$rsID[i]
  }
}

variant_oi_table_plot <- variants_oi_table

######Look at variants of interest########
#Clean up variant table before looking for variants of interest
variants_oi_table[is.na(variants_oi_table)] <- "."
variants_oi_table$Genotype <- gsub("^.$", "./.", variants_oi_table$Genotype)
variants_oi_table$rsID <- gsub("\\*", "", variants_oi_table$rsID)
variants_oi_table$rsID[grep("pos", variants_oi_table$rsID)] <- "."
variants_oi_table <- subset(variants_oi_table, select = -c(is.SGE, Conditions, rsID_annot))
variants_oi_table$Exonic_Function.refGene <- str_to_title(variants_oi_table$Exonic_Function.refGene)
variants_oi_table$Exonic_Function.refGene <- gsub("_", " ", variants_oi_table$Exonic_Function.refGene)
variants_oi_table$Exonic_Function.refGene <- gsub(" snv", "", variants_oi_table$Exonic_Function.refGene)
variants_oi_table$SGE.name <- gsub("\\:tag", "", variants_oi_table$SGE.name)
variants_oi_table$SGE.name <- gsub("\\:core", "", variants_oi_table$SGE.name)


variant_zoom_table <- filter(variants_oi_table, variants_oi_table$CADD.Phred_Normalized >= CADDCutoff |
                               variants_oi_table$SIFT_prediction == "D" |
                               variants_oi_table$Polyphen2_HVAR_Prediction == "D" |
                               variants_oi_table$Polyphen2_HVAR_Prediction == "P" |
                               variants_oi_table$SGE.name != ".")

removedVar_zoom_table <- filter(removed_variants_table, removed_variants_table$CADD.Phred_Normalized >= CADDCutoff |
                                  removed_variants_table$SIFT_prediction == "D" |
                                  removed_variants_table$Polyphen2_HVAR_Prediction == "D" |
                                  removed_variants_table$Polyphen2_HVAR_Prediction == "P" |
                                  removed_variants_table$SGE.name != ".")

#Remove variants that have allele dose of 0 or . across all variants
splitByGene <- split(variant_zoom_table, variant_zoom_table$Gene)
for(j in splitByGene){
  splitByPOS <- split(j, j$POS)
  for(i in splitByPOS){
    if(1 %in% i$Genotype_Numeric | 2 %in% i$Genotype_Numeric){
      if(!exists("variant_zoom_table_out")){
        variant_zoom_table_out <- i
      }else{
        variant_zoom_table_out <- rbind(variant_zoom_table_out, i)
      }
    }
  }
}

#Make summary table
splitByGene <- split(variant_zoom_table_out, variant_zoom_table_out$Gene)
for(i in splitByGene){
  gt_table <- subset(i, select = c(Sample, CHR, POS, Genotype, Genotype_AllelicDepth, Read_Depth))
  gt_table$SITE <- paste(gt_table$CHR, gt_table$POS, sep = ":")
  gt_table$Genotype_out <- paste(gt_table$Genotype, gt_table$Genotype_AllelicDepth, gt_table$Read_Depth, sep = ":")
  gt_table <- subset(gt_table, select = c(Sample, SITE, Genotype_out))
  gt_table <- gt_table %>% dcast(SITE ~ Sample, value.var = "Genotype_out")
  
  info_table <- subset(i, select = c(Gene, CHR, POS, REF, ALT, rsID, Function.refGene, Exonic_Function.refGene, SGE.name,
                                     CADD.Phred_Normalized, SIFT_Score, SIFT_prediction, Polyphen_HVAR_Score, Polyphen2_HVAR_Prediction))
  info_table <- unique(info_table)
  info_table$SITE <- paste(info_table$CHR, info_table$POS, sep = ":")
  info_table$SIFT_pred <- paste(info_table$SIFT_prediction, " (", info_table$SIFT_Score, ")", sep = "")
  info_table$Polyphen_pred <- paste(info_table$Polyphen2_HVAR_Prediction, " (", info_table$Polyphen_HVAR_Score, ")", sep = "")
  info_table <- subset(info_table, select = -c(CHR, POS, SIFT_prediction, SIFT_Score, Polyphen2_HVAR_Prediction, Polyphen_HVAR_Score))
  
  out <- merge(x = info_table, y = gt_table, by = "SITE", all = TRUE)
  
  if(!exists("variants_oi_summary")){
    variants_oi_summary <- out
  }else{
    variants_oi_summary <- rbind(variants_oi_summary, out)
  }
}

variants_oi_summary$SIFT_pred <- gsub("\\.\\ \\(\\.\\)", ".", variants_oi_summary$SIFT_pred)
variants_oi_summary$Polyphen_pred <- gsub("\\.\\ \\(\\.\\)", ".", variants_oi_summary$Polyphen_pred)

colnames(variants_oi_summary) <- c("Position", "Gene", "Ref", "Alt", "Variant ID", "Location",
                                   "Consequence", "Star Allele", "CADD Score (Phred Normalized)",
                                   "SIFT Prediction", "Polyphen2 Prediction","Altai Neanderthal", 
                                   "Chagyrskaya 8 Neanderthal", "Vindija 33.19 Neanderthal", "Denisovan")

#Reorder
variants_oi_summary <- variants_oi_summary[c("Gene", "Position", "Variant ID", "Location", "Consequence", "Ref", "Alt", 
                                             "Altai Neanderthal", "Chagyrskaya 8 Neanderthal", "Vindija 33.19 Neanderthal", "Denisovan",
                                             "CADD Score (Phred Normalized)","SIFT Prediction", "Polyphen2 Prediction","Star Allele")]

#Select only deleterious variants
variants_oi_summary_deleterious <- filter(variants_oi_summary, variants_oi_summary$`CADD Score (Phred Normalized)` >= CADDCutoff |
                                            grepl("D", variants_oi_summary$`SIFT Prediction`) |
                                            grepl("D", variants_oi_summary$`Polyphen2 Prediction`) |
                                            grepl("P", variants_oi_summary$`Polyphen2 Prediction`))


#Save full variant table, variant zoom table, and removed varian zoom table
if(BOOL.save){
  #Save paths
  path_summary <- paste(path.save.preffix, "nea_variantsOfInterest_summary_snpPlot_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".csv", sep = "")
  path_summary_deleterious <- paste(path.save.preffix, "nea_variantsOfInterest_summary_DELETERIOUS_snpPlot_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".csv", sep = "")
  path_full <- paste(path.save.preffix, "nea_variantsOfInterest_Full_snpPlot_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".csv", sep = "")
  path_zoom <- paste(path.save.preffix, "nea_variantsOfInterest_zoom_snpPlot_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".csv", sep = "")
  path_removed <- paste(path.save.preffix, "nea_variantsOfInterest_removed_snpPlot_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".csv", sep = "")
  
  write.table(variants_oi_summary, path_summary, sep = ",", row.names = FALSE, col.names = TRUE)
  write.table(variants_oi_summary_deleterious, path_summary_deleterious, sep = ",", row.names = FALSE, col.names = TRUE)
  write.table(variant_zoom_table_out, path_zoom, sep = ",", row.names = FALSE, col.names = TRUE)
  write.table(variants_oi_table, path_full, sep = ",", row.names = FALSE, col.names = TRUE)
  write.table(removedVar_zoom_table, path_removed, sep = ",", row.names = FALSE, col.names = TRUE)
}

#####PLOT VARIANTS#######
#Select only deleterious variants
variant_oi_table_plot <- filter(variant_oi_table_plot, !is.na(variant_oi_table_plot$Conditions))

#Make Denisovan/Neanderthal Facet
species_table <- data.frame(c("Altai Neanderthal", "Chagyrskaya 8 Neanderthal", "Denisovan", "Vindija 33.19 Neanderthal"),
                            c("Neanderthal", "Neanderthal", "Denisovan", "Neanderthal"))
colnames(species_table) <- c("Sample", "Species")
variant_oi_table_plot$Species <- species_table$Species[match(variant_oi_table_plot$Sample, species_table$Sample)]

#Clean Names
variant_oi_table_plot$Sample <- gsub(" Neanderthal", "", variant_oi_table_plot$Sample)
variant_oi_table_plot$Sample <- gsub(" 8", "", variant_oi_table_plot$Sample)
variant_oi_table_plot$Sample <- gsub(" 33.19", "", variant_oi_table_plot$Sample)

#Factor genes
variant_oi_table_plot <- variant_oi_table_plot %>%
  mutate(Gene = factor(Gene, levels=c("CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8",
                                      "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2",
                                      "CYP3A4", "CYP3A5")))

#Make custom color palette
colPalette <- c("#d3d3d3", #Light Grey (0)
                "#ff0000", #Red (1)
                "#0000FF", #Blue (2)
                "#FFFFFF") #White

#Add in Labelling
variant_oi_table_plot$Label <- variant_oi_table_plot$rsID
for(i in 1:nrow(variant_oi_table_plot)){
  input <- variant_oi_table_plot$Label[i]
  if(grepl("pos", input)){
    CHR <- variant_oi_table_plot$CHR[i]
    POS <- variant_oi_table_plot$POS[i]
    REF <- variant_oi_table_plot$REF[i]
    ALT <- variant_oi_table_plot$ALT[i]
    out <- paste(POS, ":", REF, ">", ALT, sep = "")
  }else{
    out <- input
  }
  variant_oi_table_plot$Label[i] <- out
}

#Change NA's to 'no data'
variant_oi_table_plot$Genotype_Numeric[is.na(variant_oi_table_plot$Genotype_Numeric)] <- "No Data"

#Make full plot (circle map)
alpha <- ifelse((!is.na(variant_oi_table_plot$Conditions) & variant_oi_table_plot$Genotype_Numeric > 0), "YES", "NO")
fullPlot_circleMap <- suppressWarnings(ggplot(variant_oi_table_plot, aes(Label, Sample)) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = colPalette)+
  scale_alpha_discrete(range = c(0.8, 1.0), guide=FALSE)+
  geom_point(aes(colour = as.character(Genotype_Numeric),
                 size=8,
                 alpha = alpha)) +
  geom_point(shape = 1, size = 6, colour = ifelse(variant_oi_table_plot$Genotype_Numeric == "No Data", "black", NA)) +
  scale_size(range = c(4,8))+
  facet_grid(Species~Gene, scales = "free", space = "free")+
  xlab(label = "") +
  ylab(label = "") +
  theme_bw() +
  labs(colour = "Allele Dose\n") +
  theme(plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 14, face = "bold", hjust = 1),
        axis.text.x = element_text(color = "black", size = 12, face = "bold", hjust = 1, vjust = 0.5, angle = 90),
        axis.title = element_text(size = 12, face = "bold"),
        axis.ticks  = element_blank(),
        legend.position = "right",
        panel.background = element_rect(colour = "white", fill = "white"),
        panel.border = element_rect(colour = "white"),
        strip.background.x = element_rect(fill = "white", color = "white"),
        strip.background.y = element_rect(fill = "white", color = "white"),
        strip.text.y = element_blank(),
        strip.text.x = element_text(color = "black", size = 12, face = "bold.italic", hjust = 0.5, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(0.5, "line"),
        panel.spacing.x = unit(1.1, "line")))

#Save the two plots
if(BOOL.save){
  circle_path <- paste(path.save.preffix, "variant_circlemap_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".png", sep = "")
  suppressWarnings(ggplot2::ggsave(circle_path, plot = fullPlot_circleMap, width = 11, height = 4, units = "in", dpi = 300, device = "png", limitsize = FALSE))
}

#Now count number of variants per sample/gene
variants_per_sample <- filter(variant_oi_table_plot, variant_oi_table_plot$Genotype_Numeric != 0)
splitByGene <- split(variants_per_sample, variants_per_sample$Gene)
for(i in splitByGene){
  gene <- as.character(i$Gene[1])
  #Check each of the four samples
  for(sample in c("Altai", "Chagyrskaya", "Vindija", "Denisovan")){
    vars_sample <- filter(i, i$Sample == sample)
    variant_count <- as.numeric(nrow(vars_sample))
    variant_pos <- paste(vars_sample$POS, collapse = ",")
    
    out <- data.frame(gene, sample, variant_count, variant_pos)
    
    if(!exists("out_table_vars_per_sample")){
      out_table_vars_per_sample <- out
    }else{
      out_table_vars_per_sample <- rbind(out_table_vars_per_sample, out)
    }
  }
}

out_table_vars_per_sample$gene[is.na(out_table_vars_per_sample$gene)] <- "CYP2A6"
colnames(out_table_vars_per_sample) <- c("Gene", "Sample", "Count", "Position")

#Differentiate by species
species_table <- data.frame(c("Altai", "Chagyrskaya", "Denisovan", "Vindija"),
                            c("Neanderthal", "Neanderthal", "Denisovan", "Neanderthal"))
colnames(species_table) <- c("Sample", "Species")
out_table_vars_per_sample$Species <- species_table$Species[match(out_table_vars_per_sample$Sample, species_table$Sample)]

#Factor genes
out_table_vars_per_sample <- out_table_vars_per_sample %>%
  mutate(Gene = factor(Gene, levels=c("CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8",
                                      "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2",
                                      "CYP3A4", "CYP3A5")))

fullPlot_snv <- ggplot(out_table_vars_per_sample, aes(Sample, Gene, fill= Count)) +
  geom_tile() +
  # geom_text(aes(label=Diplotype, fontface="bold"), size = 5) +
  theme_bw() +
  facet_grid(.~Species, scales = "free_x", space = "free") +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "Variant\nCount") +
  scale_y_discrete(limits = rev(levels(out_table_vars_per_sample$Gene))) +
  theme(plot.title = element_text(face = "bold", size = 25, hjust = 0.5),
        axis.text.y = element_text(color = "black", size = 15, face = "bold.italic", hjust = 1),
        axis.text.x = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.ticks  = element_blank(),
        legend.key.size = unit(2, "line"), #adjusts legend size
        legend.text = element_text(size = 12, face = "bold"),
        legend.title.align = 0.2,
        legend.title = element_text(size = 14, face = "bold", hjust = 0, vjust = 0.5),
        panel.border = element_rect(colour = "white", size = 2),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "line"))

if(BOOL.save){
  fullPlot_path <- paste(path.save.preffix, "variant_summaryMap_CADD_", CADDCutoff, "_alleleBalance_", ABR, ".png", sep = "")
  ggsave(fullPlot_path, fullPlot_snv, width = 10, height = 8, units = "in", dpi = 250)
}








