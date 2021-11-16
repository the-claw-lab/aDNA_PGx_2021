####Requires####
#Stargazer installed on computer
#Stargazer output for genes of interest
#Read depth for genes of interest

#######Libraries#######
suppressMessages(library(dplyr)) #Data Wrangling
suppressMessages(library(reshape2)) #Data Wrangling
suppressMessages(library(stringr)) #String Management
suppressMessages(library(ggplot2)) #Plotting

#Paths
preffixPath <- "~/Desktop/2021_Archaic_PGx_Github/" #"~/path/to/working/dir/"
SGEPath = "~/Software/Stargazer/" #"~/path/to/Stargazer/"

#########Import Data########
geneList <- c("cyp1a2", "cyp2a6") #Choose from "cyp1a2", "cyp2a6", "cyp2b6", "cyp2c8", "cyp2c9", "cyp2c19", "cyp2d6", "cyp2e1", "cyp2j2", "cyp3a4", "cyp3a5"
controlGenes <- c("vdr", "egfr") #can be vdr, egfr, or ryr1, must be more than 1

#GeneTable
gene_table <- read.table(paste(SGEPath, "gene_table.txt", sep = ""), sep = "\t", header = TRUE)
gene_table$chr <- gsub("chr","", gene_table$chr)

#Stargazer output
for(gene in geneList){
  for(control in controlGenes){
    path <- paste(preffixPath, "input/sge_output/sge_output_full_", control, "/", gene, "/full_", gene, "_", control, ".stargazer-genotype.txt", sep = "")
    sge_in <- read.table(path, sep = "\t", header = TRUE)
    sge_in$gene <- gene
    sge_in$control <- control
    
    if(!exists("out_table_full")){
      out_table_full <- sge_in
    }else{
      out_table_full <- rbind(out_table_full, sge_in)
    }
  }
}

#Read Depth Data
for(gene in geneList){
  path <- paste(preffixPath, "input/depth_of_coverage/neanderthal_vdr_", gene, "_out.gdf", sep = "")
  rd_in <- read.table(path, sep = "\t", header = TRUE)
  
  if(!exists("out_table_readDepth")){
    out_table_readDepth <- rd_in
  }else{
    out_table_readDepth <- rbind(out_table_readDepth, rd_in)
  }
}

#Naming Conventions
colnames(out_table_readDepth) <- c("Locus", "Total_Depth", "Average_Depth", "Altai_Neanderthal",
                                   "Chagyrskaya_8_Neanderthal", "Denisovan", "Vindija_33.19_Neanderthal")

#Clean full data table
colnames(out_table_full) <- c("Sample", "Genotype_Status", "Main_Star_Allele_1", "Main_Star_Allele_2", "Candidate_Star_Allele_1", 
                              "Candidate_Star_Allele_2", "Activity_Score_1", "Activity_Score_2", "Combined_Activity_Score",
                              "Predicted_Phenotype", "Combined_Unphased_SV", "SV_Call_1", "SV_Call_2", "Sum_Squared_Residuals",
                              "All_Possible_Candidate_Star_Alleles", "Core_SNPs_1",
                              "Core_SNPs_2", "Tag_SNPs_1", "Tag_SNPs_2", "Mean_Allele_Fraction_Gene_1", "Mean_Allele_Fraction_Gene_2",
                              "Mean_Allele_Fraction_Allele_1", "Mean_Allele_Fraction_Allele_2", "Gene", "Control")

out_table_full$Sample <- gsub("AltaiNea", "Altai_Neanderthal", out_table_full$Sample)
out_table_full$Sample <- gsub("ChagyrskayaPhalanx", "Chagyrskaya_8_Neanderthal", out_table_full$Sample)
out_table_full$Sample <- gsub("DNK02", "Denisovan", out_table_full$Sample)
out_table_full$Sample <- gsub("Vi33.19", "Vindija_33.19_Neanderthal", out_table_full$Sample)

#Make sum of squared residuals supplement table
in_table_ssr <- subset(out_table_full, select = c(Sample, Gene, Sum_Squared_Residuals, Control))
splitByControl <- split(in_table_ssr, in_table_ssr$Control)
for(i in splitByControl){
  out <- subset(i, select = -c(Control))
  out$Gene <- toupper(out$Gene)
  out <- out %>% dcast(Sample ~ Gene, value.var = "Sum_Squared_Residuals")
  out$Control_Gene <- i$Control[1]
  
  if(!exists("out_table_ssr")){
    out_table_ssr <- out
  }else{
    out_table_ssr <- rbind(out_table_ssr, out)
  }
}

########Combine identical rows from different control genes#######
#Remove ssr column
out_table_working <- subset(out_table_full, select = -c(Sum_Squared_Residuals))
splitByGene <- split(out_table_working, out_table_working$Gene)
for(i in splitByGene){
  gene <- i$Gene[1]
  #Determine if different calls are made based on control gene
  splitBySample <- split(i, i$Sample)
  for(k in splitBySample){
    #For each sample find main and candidate call based on the number of calls
    temp <- unique(k[c("Main_Star_Allele_1", "Main_Star_Allele_2", "Candidate_Star_Allele_1", "Candidate_Star_Allele_2")])
    for(l in 1:nrow(temp)){
      out <- filter(k, k$Main_Star_Allele_1 == temp$Main_Star_Allele_1[l] & 
                      k$Main_Star_Allele_2 == temp$Main_Star_Allele_2[l] &
                      k$Candidate_Star_Allele_1 == temp$Candidate_Star_Allele_1[l] &
                      k$Candidate_Star_Allele_2 == temp$Candidate_Star_Allele_2[l])
      
      controlGenes <- paste(out$Control, collapse = ",")
      if(nrow(out) > 1){
        call <- "main"
      }else{
        call <- "candidate"
      }
      out <- as.data.frame(out[1,])
      out$Control[1] <- controlGenes
      out$Final_Call <- call
      
      if(!exists("out_table_final")){
        out_table_final <- out
      }else{
        out_table_final <- rbind(out_table_final, out)
      }
    }
  }
}


#Clean the final Table
out_table_final$Activity_Score_1 <- gsub("unknown", "Unknown", out_table_final$Activity_Score_1)
out_table_final$Activity_Score_2 <- gsub("unknown", "Unknown", out_table_final$Activity_Score_2)
out_table_final$Combined_Activity_Score <- gsub("unknown", "Unknown", out_table_final$Combined_Activity_Score)
out_table_final$Predicted_Phenotype <- gsub("unknown", "Unknown", out_table_final$Predicted_Phenotype)
out_table_final$Predicted_Phenotype <- gsub("_", " ", out_table_final$Predicted_Phenotype)
out_table_final$Predicted_Phenotype <- str_to_title(out_table_final$Predicted_Phenotype)
out_table_final$Combined_Unphased_SV <- gsub(",", ", ", out_table_final$Combined_Unphased_SV)
out_table_final$Combined_Unphased_SV <- gsub("_", " ", out_table_final$Combined_Unphased_SV)
out_table_final$SV_Call_1 <- gsub("_", " ", out_table_final$SV_Call_1)
out_table_final$SV_Call_2 <- gsub("_", " ", out_table_final$SV_Call_2)
out_table_final$Gene <- toupper(out_table_final$Gene)
out_table_final$Control <- toupper(out_table_final$Control)
out_table_final$Final_Call <- str_to_title(out_table_final$Final_Call)

#Manually change vinidija cyp2c9 call to *1/*1 (this is for CYP2C9)
if("cyp2c9" %in% geneList){
  out_table_final$Main_Star_Allele_1[27] <- "*1"
  out_table_final$Main_Star_Allele_2[27] <- "*1"
  out_table_final$Candidate_Star_Allele_1[27] <- "*1"
  out_table_final$Candidate_Star_Allele_2[27] <- "*1"
  out_table_final$Activity_Score_1[27] <- "1"
  out_table_final$Activity_Score_2[27] <- "1"
  out_table_final$Combined_Activity_Score[27] <- "2"
  out_table_final$Predicted_Phenotype[27] <- "Normal Metabolizer"
  out_table_final$All_Possible_Candidate_Star_Alleles[27] <- "*1"
  out_table_final$Core_SNPs_1[27] <- "."
  out_table_final$Core_SNPs_2[27] <- "."
}

#Save overview table with every call and SSR table
if(!dir.exists(paste(preffixPath, "output/", sep = ""))) dir.create(paste(preffixPath, "output/", sep = ""))
path.output = paste(preffixPath, "output/2021_Archaic_PGx_SGE_Output_Analysis_EXAMPLE_OUTPUT/", sep = "")
if(!dir.exists(path.output)) dir.create(path.output)
path_sge <- paste(path.output,"full_sge_table.csv", sep = "")
path_ssr <- paste(path.output, "ssr_sge_table.csv", sep = "")

write.table(out_table_final, path_sge, sep = ",", col.names = TRUE, row.names = FALSE)
write.table(out_table_ssr, path_ssr, sep = ",", col.names = TRUE, row.names = FALSE)

#######Create Overview Table#######
out_table_simple <- subset(out_table_final, select = c(Sample, Main_Star_Allele_1, Main_Star_Allele_2, Activity_Score_1,
                                                       Activity_Score_2, Combined_Activity_Score, Predicted_Phenotype, SV_Call_1,
                                                       SV_Call_2, All_Possible_Candidate_Star_Alleles, Core_SNPs_1, Core_SNPs_2,
                                                       Tag_SNPs_1, Tag_SNPs_2, Gene, Control, Final_Call))
out_table_simple$Diplotype <- paste(out_table_simple$Main_Star_Allele_1, out_table_simple$Main_Star_Allele_2, sep = "/")

########Create final table for each variants with read depth and call with only main calls#######
out_table_zoom <- filter(out_table_simple, out_table_simple$Final_Call == "Main")
for(i in 1:nrow(out_table_zoom)){
  #Retrieve Applicable Items
  gene <- out_table_zoom$Gene[i]
  control <- out_table_zoom$Control[i]
  sample <- out_table_zoom$Sample[i]
  diplotype <- out_table_zoom$Diplotype[i]
  main.1 <- out_table_zoom$Main_Star_Allele_1[i]
  main.2 <- out_table_zoom$Main_Star_Allele_2[i]
  as.1 <- out_table_zoom$Activity_Score_1[i]
  as.2 <- out_table_zoom$Activity_Score_2[i]
  as.com <- out_table_zoom$Combined_Activity_Score[i]
  phen <- out_table_zoom$Predicted_Phenotype[i]
  sv.1 <- out_table_zoom$SV_Call_1[i]
  sv.2 <- out_table_zoom$SV_Call_2[i]
  dip <- out_table_zoom$All_Possible_Candidate_Star_Alleles[i]
  core.1 <- out_table_zoom$Core_SNPs_1[i]
  core.2 <- out_table_zoom$Core_SNPs_2[i]
  tag.1 <- out_table_zoom$Tag_SNPs_1[i]
  tag.2 <- out_table_zoom$Tag_SNPs_2[i]
  
  #Find Chromosome
  chr <- filter(gene_table, gene_table$name == tolower(gene))$chr[1]
  
  #Make row for each unique variant
  if(main.1 == main.2){
    allele <- main.1
    as <- as.1
    sv <- sv.1
    core <- core.1
    tag <- tag.1
    
    #Check to see if there are variants
    if(core != "."){
      #First determine core variants
      cor.vars <- unlist(str_split(core, ","))
      for(var in cor.vars){
        var.type <- "Core"
        var.temp <- unlist(str_split(var, ":"))
        pos <- as.numeric(gsub("<", "", var.temp[1]))
        mut <- var.temp[2]
        location <- var.temp[5]
        impact <- var.temp[6]
        aa <- gsub(">", "", var.temp[7])
        
        #Get read depth
        loc <- paste(chr, pos, sep = ":")
        rd <- filter(out_table_readDepth, out_table_readDepth$Locus == loc)
        rd <- as.numeric(rd[sample])
        
        #out df
        out <- data.frame(gene, control, sample, diplotype, as.com, phen, dip, allele, as, sv, pos, mut, location, impact, aa, rd, var.type)
        
        if(!exists("out_table_variants")){
          out_table_variants <- out
        }else{
          out_table_variants <- rbind(out_table_variants, out)
        }
      }
      #Now determine tag variants
      if(tag != "."){
        #First determine core variants
        tag.vars <- unlist(str_split(tag, ","))
        for(var in tag.vars){
          var.type <- "Tag"
          var.temp <- unlist(str_split(var, ":"))
          pos <- as.numeric(gsub("<", "", var.temp[1]))
          mut <- var.temp[2]
          location <- var.temp[5]
          impact <- var.temp[6]
          aa <- gsub(">", "", var.temp[7])
          
          #Get read depth
          loc <- paste(chr, pos, sep = ":")
          rd <- filter(out_table_readDepth, out_table_readDepth$Locus == loc)
          rd <- as.numeric(rd[sample])
          
          #out df
          out <- data.frame(gene, control, sample, diplotype, as.com, phen, dip, allele, as, sv, pos, mut, location, impact, aa, rd, var.type)
          
          if(!exists("out_table_variants")){
            out_table_variants <- out
          }else{
            out_table_variants <- rbind(out_table_variants, out)
          }
        }
      }
    }else{
      pos <- "."
      mut <- "."
      location <- "."
      impact <- "."
      aa <- "."
      rd <- "."
      var.type <- "."
      
      #out df
      out <- data.frame(gene, control, sample, diplotype, as.com, phen, dip, allele, as, sv, pos, mut, location, impact, aa, rd, var.type)
      
      if(!exists("out_table_variants")){
        out_table_variants <- out
      }else{
        out_table_variants <- rbind(out_table_variants, out)
      }
    }
  }else{
    for(num in c(1,2)){
      allele <- eval(parse(text = paste("main.", num, sep = "")))
      as <- eval(parse(text = paste("as.", num, sep = "")))
      sv <- eval(parse(text = paste("sv.", num, sep = "")))
      core <- eval(parse(text = paste("core.", num, sep = "")))
      tag <- eval(parse(text = paste("tag.", num, sep = "")))
      
      #Check to see if there are variants
      if(core != "."){
        #First determine core variants
        cor.vars <- unlist(str_split(core, ","))
        for(var in cor.vars){
          var.type <- "Core"
          var.temp <- unlist(str_split(var, ":"))
          pos <- as.numeric(gsub("<", "", var.temp[1]))
          mut <- var.temp[2]
          location <- var.temp[5]
          impact <- var.temp[6]
          aa <- gsub(">", "", var.temp[7])
          
          #Get read depth
          loc <- paste(chr, pos, sep = ":")
          rd <- filter(out_table_readDepth, out_table_readDepth$Locus == loc)
          rd <- as.numeric(rd[sample])
          
          #out df
          out <- data.frame(gene, control, sample, diplotype, as.com, phen, dip, allele, as, sv, pos, mut, location, impact, aa, rd, var.type)
          
          if(!exists("out_table_variants")){
            out_table_variants <- out
          }else{
            out_table_variants <- rbind(out_table_variants, out)
          }
        }
        #Now determine tag variants
        if(tag != "."){
          #First determine core variants
          tag.vars <- unlist(str_split(tag, ","))
          for(var in tag.vars){
            var.type <- "Tag"
            var.temp <- unlist(str_split(var, ":"))
            pos <- as.numeric(gsub("<", "", var.temp[1]))
            mut <- var.temp[2]
            location <- var.temp[5]
            impact <- var.temp[6]
            aa <- gsub(">", "", var.temp[7])
            
            #Get read depth
            loc <- paste(chr, pos, sep = ":")
            rd <- filter(out_table_readDepth, out_table_readDepth$Locus == loc)
            rd <- as.numeric(rd[sample])
            
            #out df
            out <- data.frame(gene, control, sample, diplotype, as.com, phen, dip, allele, as, sv, pos, mut, location, impact, aa, rd, var.type)
            
            if(!exists("out_table_variants")){
              out_table_variants <- out
            }else{
              out_table_variants <- rbind(out_table_variants, out)
            }
          }
        }
      }else{
        pos <- "."
        mut <- "."
        location <- "."
        impact <- "."
        aa <- "."
        rd <- "."
        var.type <- "."
        
        #out df
        out <- data.frame(gene, control, sample, diplotype, as.com, phen, dip, allele, as, sv, pos, mut, location, impact, aa, rd, var.type)
        
        if(!exists("out_table_variants")){
          out_table_variants <- out
        }else{
          out_table_variants <- rbind(out_table_variants, out)
        }
      }
    }
  }
}


#Clean out table
out_table_variants$location <- gsub("_", " ", out_table_variants$location)
out_table_variants$location <- str_to_title(out_table_variants$location)
out_table_variants$impact <- gsub("_", " ", out_table_variants$impact)
out_table_variants$impact <- str_to_title(out_table_variants$impact)
out_table_variants$aa <- gsub("_", " ", out_table_variants$aa)
out_table_variants$aa <- gsub("no effect", "No Effect", out_table_variants$aa)
out_table_variants$sample <- gsub("_", " ", out_table_variants$sample)

colnames(out_table_variants) <- c("Gene","Consensus Control Genes", "Sample", "Diplotype", "Combined Activity Score",
                                  "Phenotype", "All Possible Candidate Star Alleles", "Star Allele", "Activity Score",
                                  "Structural Variation", "Position", "Mutation", "Variant Location",  "Impact",
                                  "Amino Acid Change", "Position Read Depth", "Variant Type")

#Work on making final table
out_table_variants_summary <- out_table_variants

#Clean up table before formatting
out_table_variants_summary$`Structural Variation` <- toupper(out_table_variants_summary$`Structural Variation`)
out_table_variants_summary$`Structural Variation` <- gsub("NO", "No", out_table_variants_summary$`Structural Variation`)
out_table_variants_summary$`Structural Variation` <- gsub("CNV.", "CNV", out_table_variants_summary$`Structural Variation`)
out_table_variants_summary$`Structural Variation` <- gsub("GC E1E2", "Gene-Hybrid", out_table_variants_summary$`Structural Variation`)
out_table_variants_summary$`Consensus Control Genes` <- gsub("EGFR", "E", out_table_variants_summary$`Consensus Control Genes`)
out_table_variants_summary$`Consensus Control Genes` <- gsub("VDR", "V", out_table_variants_summary$`Consensus Control Genes`)
out_table_variants_summary$`Consensus Control Genes` <- gsub("RYR1", "R", out_table_variants_summary$`Consensus Control Genes`)

col_order <- c("Gene", "Sample", "Diplotype", "Phenotype", "Combined Activity Score", "All Possible Candidate Star Alleles",  "Consensus Control Genes",          
               "Star Allele", "Activity Score", "Structural Variation", "Position", "Mutation", "Variant Location", "Impact", "Amino Acid Change",                  
               "Position Read Depth", "Variant Type")
out_table_variants_summary <- out_table_variants_summary[, col_order]

splitByGene <- split(out_table_variants_summary, out_table_variants_summary$Gene)
for(i in splitByGene){
  #Now make rows for each sample
  splitBySample <- split(i, i$Sample)
  for(j in splitBySample){
    #Make sample row
    sample_row <- j[1,]
    sample_row[1, c(8:17)] <- ""
    #Make variant specific rows
    for(k in 1:nrow(j)){
      var_row <- j[k,]
      var_row[1, c(1:7)] <- ""
      
      if(!exists("var_rows")){
        var_rows <- var_row
      }else{
        var_rows <- rbind(var_rows, var_row)
      }
    }
    
    #Now merge sample header and var rows
    sample_rows <- rbind(sample_row, var_rows)
    rm(var_rows)
    
    if(!exists("gene_rows")){
      gene_rows <- sample_rows
    }else{
      gene_rows <- rbind(gene_rows, sample_rows)
    }
  }
  #Now make final set of rows per gene
  rows_per_gene <- gene_rows
  rm(gene_rows)
  #save everything in out table
  if(!exists("out_table_variants_clean")){
    out_table_variants_clean <- rows_per_gene
  }else{
    out_table_variants_clean <- rbind(out_table_variants_clean, rows_per_gene)
  }
}

#Save summary / publication table
path_sum <- paste(path.output, "summary_sge_table.csv", sep = "")
path_full <- paste(path.output,"out_table_variants_full.csv", sep = "")

write.table(out_table_variants_clean, path_sum, sep = ",", col.names = TRUE, row.names = FALSE)
write.table(out_table_variants, path_full, sep = ",", col.names = TRUE, row.names = FALSE)

#########Make Final Plot########
#Set Factor
phenotypes <- c("Poor Metabolizer","Poor Function", "Slow Metabolizer",
                "Intermediate Metabolizer",
                "Unfavorable Response","Decreased Function",
                "Normal Metabolizer", "Normal Function",
                "Increased Function", "Favorable Response",
                "Rapid Metabolizer",
                "Ultrarapid Metabolizer",
                "Unknown")
out_table_variants$Phenotype <- factor(out_table_variants$Phenotype, levels = phenotypes)

#Make my color palette
colPalette <- c("#D55E00", "#D55E00", #Red
                "#E69F00", #Orange
                "#F0E442", "#F0E442", "#F0E442", #Yellow
                "#009E73", "#009E73", #Green
                "#56B4E9", "#56B4E9", "#56B4E9","#56B4E9", #Blue
                "#999999") #Grey
colPalette = setNames(colPalette, levels(out_table_variants$Phenotype))

#Make Denisovan/Neanderthal Facet
species_table <- data.frame(c("Altai Neanderthal", "Chagyrskaya 8 Neanderthal", "Denisovan", "Vindija 33.19 Neanderthal"),
                            c("Neanderthal", "Neanderthal", "Denisovan", "Neanderthal"))
colnames(species_table) <- c("Sample", "Species")
out_table_variants$Species <- species_table$Species[match(out_table_variants$Sample, species_table$Sample)]

#Change Sample names
out_table_variants$Sample <- gsub(" Neanderthal", "", out_table_variants$Sample)
out_table_variants$Sample <- gsub(" 8", "", out_table_variants$Sample)
out_table_variants$Sample <- gsub(" 33.19", "", out_table_variants$Sample)

out_table_plot <- unique(subset(out_table_variants, select = c(Sample, Gene, Phenotype, Diplotype, Species)))

#Factor genes
out_table_plot <- out_table_plot %>%
  mutate(Gene = factor(Gene, levels=c("CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8",
                                      "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP2J2",
                                      "CYP3A4", "CYP3A5")))

#Change CYP2C9 Vindija to "N.D." #This is for CYP2C9
if("cyp2c9" %in% geneList){
  out_table_plot$Phenotype[24] <- "Unknown"
  out_table_plot$Diplotype[24] <- "N.D."
}

#Plot all genes/samples in heatmap
fullPlot_snv <- ggplot(out_table_plot, aes(Sample, Gene, fill= Phenotype)) +
  scale_fill_manual(values = colPalette) +
  geom_tile() +
  geom_text(aes(label=Diplotype, fontface="bold"), size = 5) +
  theme_bw() +
  facet_grid(.~Species, scales = "free_x", space = "free") +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "Phenotype\n") +
  scale_y_discrete(limits = rev(levels(out_table_plot$Gene))) +
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

plotPath <- paste(path.output, "full_sge_plot.png", sep = "")
ggsave(plotPath, fullPlot_snv, width = 10, height = 8, units = "in", dpi = 350)




