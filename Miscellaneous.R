##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##******************************************************************************
## 
## Emine Ozsahin
## Analysis of Splicing in weat
## Summer 2020

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#libraries ---
library(qtl2)
library(tidyr)
library(tibble)

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Input Files for genetic map (wheat_gmap.csv) and genotypes (wheat_geno.csv) for qtl mapping 

# Input phenotype files (wheat_pheno.csv) were prepared for each chromosome. 
# Therefore, codes for preparation of phenotype_files are in the rqtl.R script 

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 1- gmap: wheat_gmap.csv ----
# Following codes were used to prepare the genetic map file called wheat_gmap.csv. 
# The file SNP-location2.csv was the genetic map which includes markers (SNPs) and their physical and genetic positions. 

genetic_map <- read.csv("SNP-location2.csv", header = TRUE, stringsAsFactors = FALSE) #modified excel file Lewis send it to me 
genetic_map[1:2,]

genetic_map$Chr <- genetic_map$Chr %>%
  str_replace_all(c("chr" = ""))

wheat_gmap <- data.frame(marker=genetic_map$Marker, chr= genetic_map$Chr, pos=genetic_map$Marker.position..CM..on.genetic.map)
wheat_gmap[1:5,1:3]

write.csv(wheat_gmap, "wheat_gmap.csv", row.names = FALSE)

##******************************************************************************
# 2- geno: wheat_geno.csv ----
# A genotype file was used to prepare the rqtl input genotype file called  the wheat_geno.csv. 
# Genotypes were coded by their nucleotides at the marker positions. 
# These has to be changed to A and B style. 
# Nucleotides in the final file are converted to A if the nucleotide for an individual is identical with the reference nucleotide 
# and B if the individual nucleotide differs at the position of the markers.    
# Change the values 
SNP_file_name <- read.csv("SNP.txt", sep = "\t")
colnames(SNP_file_name)
dim(SNP_file_name)

# Preparation to convert the values
ref <- SNP_file_name[,4]
genotype <- SNP_file_name[,-c(1:4)]
length(colnames(genotype[1,])) #156

# Convert factors to characters
ref <- lapply(ref, as.character)  
genotype <- apply(genotype, 2, as.character)

# 0|0 means that individual is homozygous for the reference sequence at that position
# 1|1 means that individual is homozygous for the alternative sequence at that position
## Change the Nucleotides to 0 or 1: if the individuals have the reference (REF) nucleotide change it to 0, if not 1
SNP_file_name_pre <- data.frame()
for (index in 1:length(ref)) { out <- ifelse(genotype[index,]==ref[index], 0, 1)
  SNP_file_name_pre <- rbind(SNP_file_name_pre, out)
  } 

colnames(SNP_file_name_pre) <- colnames(genotype)
colnames(SNP_file_name_pre) <- gsub("LL2017_", "", colnames(SNP_file_name_pre))
colnames(SNP_file_name_pre) <- gsub("L2017_", "", colnames(SNP_file_name_pre))
SNP_file_name_pre <- cbind(snpid = SNP_file_name$markerID, SNP_file_name_pre)
which(SNP_file_name_pre == "chr1_mA:64141585")

# delete the "_part[0-9]" from snp names    
SNP_file_name_pre$snpid <- gsub("_part[0-9]_", ":", SNP_file_name_pre$snpid) 

write.table(SNP_file_name_pre, "SNP01.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Use the SNP01.txt to change 0 and 1s to A and Bs
geno <- read.csv("SNP01.txt", sep = "\t")
geno <- geno[, -c(2,3,4,5,6)]

geno$markerID <- geno$markerID %>%
  str_replace_all(c("_part[0-9]" = ""))

geno$snpid <- gsub(":", "_", geno$snpid)
colnames(geno) <- gsub("X", "", colnames(geno))
colnames(geno)[1] = "markerID"

cnames<- unlist(geno$markerID)
geno <- as.data.frame(t(geno))
colnames(geno) <- cnames
geno <- geno[-1,]
rnames <- rownames(geno)

geno <- apply(geno, 2, function(x) gsub(1, "A", x))
geno <- apply(geno, 2, function(x) gsub(0, "B", x))
geno <- apply(geno, 2, function(x) gsub(" ", "", x))

geno <- as.data.frame(geno)
rownames(geno)[1] = "id"
geno <- rownames_to_column(geno, var = "id")
colnames(geno) <- c( "id", cnames)

# change NAs to "-"
geno <- as.data.frame(ifelse(apply(geno, 2, is.na) == TRUE, "-", apply(geno, 2, print)))

write.csv(geno, "wheat_geno.csv", row.names = FALSE)

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Annotations- Step 1 assign GeneIDs to the introns

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# first need a gene table which includes gene names and locations, I prepared it from exons_table (I could have extracted the genes from gff3 file)
genes <- unique(exons_table$gene_name)
genes_table <- data.frame()
for (i in 1:length(genes)){
  gene_loc <- exons_table[which(exons_table$gene_name == genes[i]), ]
  gene_start <- gene_loc$start[1]
  gene_end <- gene_loc$end[dim(gene_loc)[1]]
  temp <- data.frame(gene = genes[i], start=gene_start, end =  gene_end, chr = gene_loc$chr)
  genes_table <- rbind(genes_table, temp)}

write.csv(genes_table, "wheat_genes_table.csv")

#grep "gene"  wheat_genes_table.csv > wheat_gene_table_no_noncoding.csv

#sed "s/,/    /g"  wheat_gene_table_no_noncoding.csv |cut -f2 -f3 -f4 -f5 |sort |uniq > wheat_gene_table_no_noncoding_uniq.csv

genes_table <- read.csv("wheat_gene_table_no_noncoding_uniq.csv", sep = "\t")

# Assign the GeneIDs to the introns 
qtl$geneID2 <- 1:length(qtl$clusterID) 
for (i in 1:length(qtl$geneID2)) { 
  region <- genes_table[which(qtl$intron_chr[i] == genes_table$chr), ] 
  gene_region <- region[which(qtl$intron_start[i] > region$start & 
                                qtl$intron_end[i] < region$end & 
                                qtl$intron_start[i] < region$end & 
                                qtl$intron_end[i] > region$start),]
  value <- length(unique(gene_region$gene)) 
  qtl$geneID2[i] <- if (value == 1 )  as.character(gene_region$gene) else 
    if (value == 0 | value > 1 ) {gene_region <- region[which(qtl$intron_start[i] < region$start
                                                              & qtl$intron_end[i] < region$end 
                                                              & qtl$intron_start[i] < region$end
                                                              & qtl$intron_end[i] > region$start),]
    value <- length(unique(gene_region$gene)) 
    qtl$geneID2[i] <- if (value == 1 )  as.character(gene_region$gene) else 
      if (value == 0 | value > 1 ) {gene_region <- region[which(qtl$intron_start[i] > region$start & 
                                                                  qtl$intron_end[i] > region$end &
                                                                  qtl$intron_start[i] < region$end & 
                                                                  qtl$intron_end[i] > region$start ),]
      value <- length(unique(gene_region$gene)) 
      qtl$geneID2[i] <- if (value == 1 )  as.character(gene_region$gene) else value }
    }
  }




##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Fill the missing QTL_types with "trans" if the intron and QTL chromosomes are different.   

sQTL_Results$qtl_type2 <- 1:length(sQTL_Results$qtl_type)
for (i in 1:length(sQTL_Results$qtl_type)) {
  sQTL_Results$qtl_type2[i] <- if (sQTL_Results$intron_chr[i] != sQTL_Results$qtl_chr[i]) "trans" else as.character(sQTL_Results$qtl_type[i])
}

na_qtl_type <- sQTL_Results[which(is.na(sQTL_Results$qtl_type2)),]
head(na_qtl_type[,-c(15,19,26,27)])

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Gene Ontology Graph
# Gene Ontology terms retrieved from http://geneontology.org
# specifically http://pantherdb.org/webservices/go/overrep.jsp
# 

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene names to search for gene ontology
# File obtained from Panther 

BP <- read.csv("BP_pantherChart.txt", sep="\t", header = FALSE)
CC <- read.csv("CC_pantherChart.txt", sep="\t", header = FALSE)
MF <- read.csv("MF_pantherChart.txt", sep="\t", header = FALSE)
MF$type <- c(rep("MF", 5))
CC$type <- c(rep("CC", 10))
BP$type <- c(rep("BP", 10))
GO <- rbind(BP,CC,MF)
colnames(GO) <- c("index", "GO_term", "Number of sGenes", "V4", "V5", "GO class")
GO$GO_term <- gsub(" \\(GO:[0-9]+)", "", GO$GO_term)

library(ggpubr)
png("GO_sQTL.png", units="in", width=10, height=6, res=300)
ggbarplot(GO, x = "GO_term", y = "Number of sGenes",
          fill = "GO class",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          orientation = "horiz",
          lab.size = 0.005
)

dev.off()

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Intron selection

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


library(stringr)
#sQTL <- cbind(sQTL, str_split_fixed(sQTL$lodcolumn, ":", 4))
unique_clusterIDs <- unique(sQTL$clusterID.1)

# check the LOD scores if they are sorted or you may apply this to individual peaks file--> try it

# use intron start position, qtl location and position to filter the introns

sQTL_intron_filter <- data.frame()
for (i in 1:length(unique_clusterIDs)) {
  introns <- sQTL[which(sQTL$clusterID.1 == unique_clusterIDs[i]), ]
  #select_start <- unique(introns[,c(9,18,23)])
  select_start <- unique(introns[,c(9,28)])
  temp <- introns[rownames(select_start),]
  select_end <- unique(temp[,c(10,28)])
  temp <- introns[rownames(select_end),]
  sQTL_intron_filter <- rbind(sQTL_intron_filter, temp)}

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# QTL additive effects

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sqtl_intron_filtered$additive_effect <- 1:length(sqtl_intron_filtered$intron)
for (i in 1:length(sqtl_intron_filtered$intron)) { 
  if (is.na(sqtl_intron_filtered$A.pheno.mean[i])) { sqtl_intron_filtered$additive_effect[i] <- NA } else {
    sqtl_intron_filtered$additive_effect[i] <- if (sqtl_intron_filtered$A.pheno.mean[i] > sqtl_intron_filtered$B.pheno.mean[i])
      (sqtl_intron_filtered$A.pheno.mean[i]-sqtl_intron_filtered$B.pheno.mean[i])/2 else (sqtl_intron_filtered$B.pheno.mean[i] - sqtl_intron_filtered$A.pheno.mean[i])/2 
    }}

write.csv(sqtl_intron_filtered, "sqtl_intron_filtered_additive.csv")

length(which(sqtl_intron_filtered$additive_effect > 1))
# 1642 

additive_bigger_than_one <- sqtl_intron_filtered[which(sqtl_intron_filtered$additive_effect > 1),]

length(unique(additive_bigger_than_one$GeneID)) #690

head(additive_bigger_than_one[which(additive_bigger_than_one$effect == "BB"),])
additive_bigger_than_one <- sqtl_intron_filtered[which(sqtl_intron_filtered$additive_effect > 1),]

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Genes have the same qtl

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sGenes <- sqtl_intron_filtered[grep("Traes", sqtl_intron_filtered$GeneID),]
unique(sGenes[,c(16,21)])

write.csv(sGenes, "sGenes.csv")

# I continued with python on jupyter lab, 
#import pandas as pd
#df_csv = pd.read_csv("sGenes.csv", header=0, index_col=0, quotechar='"', sep=",", na_values=["nan", "-", ".", " ", ""], delimiter=',')
#
##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Venn Diagram for QTL effects
#AA --> RF, BB --> Stettler

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# cis sQTL were grouped for QTL effect
annot <- read.csv("sQTL_Results_Report.csv")
annot<- annot[which(annot$qtl_type == "cis"),]

# Count the QTL regions for each gene
unique_genes <- unique(annot$GeneID)
sGenes_qtl_num <- data.frame()
for (i in 1:length(unique_genes)){
  gene <- unique_genes[i]
  subset <- annot[which(annot$GeneID == gene),]  # subset with the gene take all the rows has similar gene
  num <- length(unique(subset$qtl_loc_2))  # count haw many unique sQTL region 
  temp <- data.frame(gene = gene, num_of_qtl = num) 
  sGenes_qtl_num <- rbind(sGenes_qtl_num, temp)
}

# Extract the gene names only the genes has one sQTL region
genes_only_one_cis <- sGenes_qtl_num$gene[which(sGenes_qtl_num$num_of_qtl == 1)]
# the variable "genes_only_one_cis" contains only the gene names

# Subset the original table with only containing the genes has one sQTL region
genes_only_one_cis_all <- annot[unlist(lapply(genes_only_one_cis, function(x) which(annot$GeneID == x))),]  #2153 11
# Grouped them by their QTL_effect: AA means splicing are higher in Red Fife than Stettler
# BB means splicing are higher than in Stettle than Red Fife
genes_only_one_cis_AA <- genes_only_one_cis_all[which(genes_only_one_cis_all$effect == "AA"),]
genes_only_one_cis_BB <-  genes_only_one_cis_all[which(genes_only_one_cis_all$effect == "BB"),]

# Repeat the same prosedures for trans sQTL
annot <- read.csv("sQTL_Results_Report.csv")
annot <- annot[which(annot$qtl_type == "trans"),]

# Count the QTL regions for each gene
unique_genes<- unique(annot$GeneID)
sGenes_qtl_num <- data.frame()
for (i in 1:length(unique_genes)){
  gene <- unique_genes[i]
  subset <- annot[which(annot$GeneID == gene),]  # subset with the gene take all the rows has similar gene
  num <- length(unique(subset$qtl_loc_2))  # count haw many unique sQTL region 
  temp <- data.frame(gene = gene, num_of_qtl = num) 
  sGenes_qtl_num <- rbind(sGenes_qtl_num, temp)
}

genes_only_one_trans <- sGenes_qtl_num$gene[which(sGenes_qtl_num$num_of_qtl == 1)]
genes_only_one_trans_all <- annot[unlist(lapply(genes_only_one_trans, function(x) which(annot$GeneID == x))),]
genes_only_one_trans_AA <- genes_only_one_trans_all[which(genes_only_one_trans_all$effect == 'AA'),]
genes_only_one_trans_BB <- genes_only_one_trans_all[which(genes_only_one_trans_all$effect == 'BB'),]


length(genes_only_one_cis_AA$GeneID)
length(genes_only_one_cis_BB$GeneID)
length(genes_only_one_trans_AA$GeneID)
length(genes_only_one_trans_BB$GeneID)

# Venn Diagram with genes
library("VennDiagram")
venn.diagram (x = list(x1 = genes_only_one_cis_AA$GeneID, 
                       x2 = genes_only_one_cis_BB$GeneID,
                       x3 = genes_only_one_trans_AA$GeneID, 
                       x4 = genes_only_one_trans_BB$GeneID),
              category.names = c("Red Fife cis-sQTL", 
                                 "Stettler cis-sQTL",
                                 "Red Fife trans-sQTL",
                                 "Stettler trans-sQTL"),
              fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"), 
              main = "sGenes", main.pos = c(0.5, 1.05), main.just = c( 0.5, 1), 
              lwd = 2, lty = "solid", col = "black", 
              alpha = 0.5, rotation.degree = 0, 
              rotation.centre = c(0.5, 0.5), 
              label.col = "black", cex = 1, fontface = "plain", 
              fontfamily = "serif", 
              cat.dist = 0.10, cat.cex = 1, cat.col = "black", 
              cat.fontface = "plain", 
              cat.fontfamily = "serif", cat.prompts = FALSE,
              ext.text = FALSE, euler.d = TRUE, 
              scaled = TRUE, sep.dist = 0.05, 
              offset = 0, height = 6, width = 6, 
              resolution = 1000,
              units = "in", 
              filename = "genes.tiff",
              sub.pos = c( 0.5, 1.05), # bu lazim olmayabilir bir bak sonra
              ext.pos = "", ext.percent = "", ext.line.lwd = "", ext.line.lty = "", description = "",  ext.dist = "", ext.length = "")



# Venn Diagram with introns
venn.diagram (x = list(x1 = genes_only_one_cis_AA$intron, 
                       x2 = genes_only_one_cis_BB$intron,
                       x3 = genes_only_one_trans_AA$intron, 
                       x4 = genes_only_one_trans_BB$intron),
              category.names = c("Red Fife cis-sQTL", 
                                 "Stettler cis-sQTL",
                                 "Red Fife trans-sQTL",
                                 "Stettler trans-sQTL"),
              fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"), 
              main = "Introns", main.pos = c(0.5, 1.05), main.just = c( 0.5, 1), 
              lwd = 2, lty = "solid", col = "black", 
              alpha = 0.5, rotation.degree = 0, 
              rotation.centre = c(0.5, 0.5), 
              label.col = "black", cex = 1, fontface = "plain", 
              fontfamily = "serif", 
              cat.dist = 0.10, cat.cex = 1, cat.col = "black", 
              cat.fontface = "plain", 
              cat.fontfamily = "serif", cat.prompts = FALSE,
              ext.text = FALSE, euler.d = TRUE, 
              scaled = TRUE, sep.dist = 0.05, 
              offset = 0, height = 6, width = 6, 
              resolution = 1000,
              units = "in", 
              filename = "introns.tiff",
              sub.pos = c( 0.5, 1.05), # bu lazim olmayabilir bir bak
              ext.pos = "", ext.percent = "", ext.line.lwd = "", ext.line.lty = "", description = "",  ext.dist = "", ext.length = "")


# Find the intersecting genes
overlap <- get.venn.partitions( x = list("AAcis" = genes_only_one_cis_AA$GeneID, 
                              "BBcis" = genes_only_one_cis_BB$GeneID, 
                              "BBtrans" = genes_only_one_trans_BB$GeneID, 
                              "AAtrans" = genes_only_one_trans_AA$GeneID));

overlap[1:15,1:5]

# AAcis BBcis BBtrans AAtrans                         ..set..
# 1   TRUE  TRUE    TRUE    TRUE     AAcis∩BBcis∩BBtrans∩AAtrans
# 2  FALSE  TRUE    TRUE    TRUE (BBcis∩BBtrans∩AAtrans)∖(AAcis)
# 3   TRUE FALSE    TRUE    TRUE (AAcis∩BBtrans∩AAtrans)∖(BBcis)
# 4  FALSE FALSE    TRUE    TRUE (BBtrans∩AAtrans)∖(AAcis∪BBcis)
# 5   TRUE  TRUE   FALSE    TRUE (AAcis∩BBcis∩AAtrans)∖(BBtrans)
# 6  FALSE  TRUE   FALSE    TRUE (BBcis∩AAtrans)∖(AAcis∪BBtrans) 11 #opposing
# 7   TRUE FALSE   FALSE    TRUE (AAcis∩AAtrans)∖(BBcis∪BBtrans) 13 #opposing
# 8  FALSE FALSE   FALSE    TRUE (AAtrans)∖(AAcis∪BBcis∪BBtrans)
# 9   TRUE  TRUE    TRUE   FALSE (AAcis∩BBcis∩BBtrans)∖(AAtrans)
# 10 FALSE  TRUE    TRUE   FALSE (BBcis∩BBtrans)∖(AAcis∪AAtrans) 19 #reinforcing
# 11  TRUE FALSE    TRUE   FALSE (AAcis∩BBtrans)∖(BBcis∪AAtrans) 9  #reinforcing
# 12 FALSE FALSE    TRUE   FALSE (BBtrans)∖(AAcis∪BBcis∪AAtrans)
# 13  TRUE  TRUE   FALSE   FALSE (AAcis∩BBcis)∖(BBtrans∪AAtrans)
# 14 FALSE  TRUE   FALSE   FALSE (BBcis)∖(AAcis∪BBtrans∪AAtrans)
# 15  TRUE FALSE   FALSE   FALSE (AAcis)∖(BBcis∪BBtrans∪AAtrans)

summary(overlap)

#7, 10, 11, 6

gene_select <- c(as.character(unlist(overlap[6,][6])), as.character(unlist(overlap[7,][6])), as.character(unlist(overlap[10,][6])), as.character(unlist(overlap[11,][6])))
length(gene_select)
gene_select 


# To use in functional analysis
write.csv(gene_select, "gene_select.csv")
#sed 's/,/ /g' gene_select.txt |cut -f2 |sed 's/"//g'
# continue to AmiGO

# subset the variable contains all results 
gene_select_with_lod <- annot[ unlist(lapply(gene_select, function(x) which(annot$GeneID == x))),]
summary(gene_select_with_lod)

inf9 <- gene_select_with_lod[unlist(lapply(gsub("[0-9] ", "", reinf9$x), function(x) which(gene_select_with_lod$gene == x))),]
inf19 <- gene_select_with_lod[unlist(lapply(gsub("[0-9] ", "", reinf19$x), function(x) which(gene_select_with_lod$gene == x))),]
reinforcing <- rbind(inf9, inf19)

opposing <- gene_select_with_lod[unlist(lapply(reinforcing$gene, function(x) which(gene_select_with_lod$gene == x))),]


write.table(as.character(unlist(overlap[6,][6])), "venn11.txt")
#sed 's/ / /g' venn6.txt |cut -f2 |sed 's/"//g'

write.table(as.character(unlist(overlap[7,][6])), "venn13.txt")
#sed 's/ / /g' venn6.txt |cut -f2 |sed 's/"//g'

write.table(as.character(unlist(overlap[10,][6])), "venn19.txt")
#sed 's/ / /g' venn6.txt |cut -f2 |sed 's/"//g'

write.table(as.character(unlist(overlap[11,][6])), "venn9.txt")
#sed 's/ / /g' venn6.txt |cut -f2 |sed 's/"//g'

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##BioMart
library(biomaRt)

# R package biomart do not work with queuing system. It is not connected to the server.
# Interactively it is disconnecting after almost 14000 search (I had 29000 introns to annotate)
# Therefore, I used it by dividing the file into 3. Or you may do this for each chromosomes.

# To work in quening system with Biomart there are other versions like MartShell etc.

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# following for the genes found in differentially alternatively spliced introns analyses 
# Leafcutter produced the variable 'introns'

introns$region <-  1:length(introns$gene)
for (i in 1:length(introns$gene)) { 
  region <- paste(introns$chr[i], introns$start[i], introns$end[i], sep=":") 
  introns$region[[i]] <- region
  }

regions <- introns$region <- gsub("chr", "", introns$region)

mart_plant = useMart(biomart = "plants_mart", 
                     host = "https://plants.ensembl.org", 
                     port = 443)

mart_plant <- useDataset(mart = mart_plant, dataset = "taestivum_eg_gene")

filter <- "chromosomal_region"

attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position" )

bio <- data.frame()
per <- data.frame() # this intermediate is necerrary to prevent having error
for (i in 1:length(reagions)) {temp <- getBM(attributes=attributes, filters=filter, values=reagions[i], mart=mart_plant)
  temp <- if (dim(temp)[1] == 0) data.frame(a=NA,b=NA,c=NA,d=NA) else temp
  colnames(temp) <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
  per <- rbind(temp, per)
  bio <- rbind(per, bio)
  per <- data.frame()
  }

cat("writing the results...")
write.csv(cbind(qtl, bio), "qtl_with_annot.csv")

# my qtl file for all chromosomes.
qtl <- read.csv("all_chrom_qtl_effect.csv", sep=",", header=FALSE)
# intron regions 
reagions <- paste(qtl$V8, qtl$V9, qtl$V10, sep=":")
colnames(qtl) <- c("", "intron","AA","BB","A.pheno.mean","B.pheno.mean","effect","intron_chr","intron_start","intron_end","qtl_chr","qtl_start","qtl_end","qtl_type","lod")


##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##BioCircos
library(BioCircos)
library(RColorBrewer)
library(grDevices)
library(stringr)

##******************************************************************************
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Define the genome
genomeChr = c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Un")

lengthChr = c(594102056/1000000, 689851870/1000000, 495453186/1000000, 780798557/1000000, 801256715/1000000, 651852609/1000000, 
              750843639/1000000, 830829764/1000000, 615552423/1000000, 
              744588157/1000000, 673617499/1000000, 509857067/1000000, 
              709773743/1000000, 713149757/1000000, 566080677/1000000, 
              618079260/1000000, 720988478/1000000, 473592718/1000000, 
              736706236/1000000, 750620385/1000000, 638686055/1000000, 
              480980714/1000000)

names(lengthChr) <- genomeChr

listcsv <- sort(dir(pattern = "*peaks.csv"))

peaks <- list()
for (i in 1:length(listcsv)){
  peaks[[i]] <- read.csv(listcsv[i])
  }

p <- c()
for (i in 1:length(peaks)){
  p <- c(p, peaks[i][[1]]$lod)
  }


tracks = BioCircosTracklist()

# Add one track for each chromosome
for (i in 1:length(genomeChr)){nbBars = length(hist(peaks[i][[1]]$lod)$breaks) - 1
barValues = hist(peaks[i][[1]]$lod)$counts
names(barValues) = rep(i, nbBars)
barColor = colorRampPalette(brewer.pal(8, "YlOrBr"))(length(genomeChr))[i]
tracks = tracks + BioCircosBarTrack(paste0("bars", i), chromosome = genomeChr[i],
                                    starts = cumsum(rep(lengthChr[i]/nbBars, nbBars)) - (lengthChr[i]/nbBars),
                                    ends = cumsum(rep(lengthChr[i]/nbBars, nbBars)), values = barValues, color = barColor,
                                    range = c(min(p), max(p)), maxRadius = 0.53, minRadius = 0.5) }

print("bar tracks are generated")

# Add background
tracks = tracks + BioCircosBackgroundTrack("bars_background", colors = "#2222EE")

print("tracks are generated")

save(tracks,p, file = "BioCircos.RData")

save.image()

#BioCircos(tracks, genomeFillColor = "YlOrBr", genome = as.list(lengthChr), 
#genomeTicksDisplay = F, genomeLabelDy = 0)
#print("Biociros is done")

######################################################################333

### LINK TRACKS FOR BIOCIRCOS 
file <- read.csv("jul3rd_Results/sQTL_Results.csv", header=TRUE)
#dim(file)
#head(file)
#colSums(is.na(file[!is.na(file$qtl_start),]))
#file <- file[-which(is.na(file$qtl_start)),]
#dim(file)

#head(file)
file <- file[unique(file$clusterID.1),]
#dim(file)

chr1 <- file[which(file$intron_chr == "1D"),]
chr2 <- file[which(file$intron_chr == "3D"),]
chr3 <- file[which(file$intron_chr == "2D"),]
#chromosomes 1
links_chr1_intron = c(as.character(chr1$intron_chr))
links_chr1_qtl = c(as.character(chr1$qtl_chr))

# start positions 1
links_chr1_intron_start = chr1$intron_start/1000000
links_chr1_qtl_start = chr1$qtl_start/1000000

# end positions 1
links_chr1_intron_end = chr1$intron_end/1000000
links_chr1_qtl_end = chr1$qtl_end/1000000

#chromosomes 2
links_chr2_intron = c(as.character(chr2$intron_chr))
links_chr2_qtl = c(as.character(chr2$qtl_chr))

# start positions 2
links_chr2_intron_start = chr2$intron_start/1000000
links_chr2_qtl_start = chr2$qtl_start/1000000

# end positions 2
links_chr2_intron_end = chr2$intron_end/1000000
links_chr2_qtl_end = chr2$qtl_end/1000000

#chromosomes 3
links_chr3_intron = c(as.character(chr3$intron_chr))
links_chr3_qtl = c(as.character(chr3$qtl_chr))

# start positions
links_chr3_intron_start = chr3$intron_start/1000000
links_chr3_qtl_start = chr3$qtl_start/1000000

# end positions
links_chr3_intron_end = chr3$intron_end/1000000
links_chr3_qtl_end = chr3$qtl_end/1000000

#tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0,35, maxRadius = 1, borderSize = 0,5, fillColors = "#EEFFEE")  

#length(links_chr1_intron)
#length(links_chr2_intron)
#length(links_chr3_intron)
#colors()
#?rgb() # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
#col_blue2 <- rgb(50, 50, 200, max = 200, alpha = 125, names = "skyblue2")
#col_red <- rgb(200, 50, 50, max = 200, alpha = 125, names = "rosybrown")
#col_green <- rgb(100, 255, 150, max = 255, alpha = 125, names = "darkgreen")
library("yarrr")
#col_green <- yarrr::transparent("darkgreen", trans.val = .7)
#col_blue <- yarrr::transparent("skyblue2", trans.val = .7)
#col_red <- yarrr::transparent("rosybrown", trans.val = .7)

col_red <- yarrr::transparent("#724E2F", trans.val = .7)
col_green <- yarrr::transparent("darkgreen", trans.val = .7)
col_blue <- yarrr::transparent("#D9870F", trans.val = .7)


# SNP track
snps <- read.csv("SNP.txt", sep = "\t")
snps <- snps[,1:3]
#head(snps)
snps$CHROM <- gsub("chr", "", snps$CHROM)
snps$CHROM <- gsub("_part[0-9]", "", snps$CHROM)

# Chromosomes on which the points should be displayed
points_chromosomes = c(as.character(gsub('"', '', snps$CHROM)))
# Coordinates on which the points should be displayed
points_coordinates = snps$POS/1000000
# Values associated with each point, used as radial coordinate 
#   on a scale going to minRadius for the lowest value to maxRadius for the highest value
points_values = c(rep(0:9, length(points_chromosomes)/10))
points_values <- c(points_values, c(1:8))
#length(points_values)
#length(points_chromosomes)

rm(tracklist)
tracklist = BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates, 
                              points_values, colors = col_blue[[1]],minRadius = 1.26, maxRadius = 1.55) #minRadius = 1.27, maxRadius = 2


tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', 
                               links_chr1_intron, links_chr1_intron_start, links_chr1_intron_end,
                               links_chr1_qtl, links_chr1_qtl_start, links_chr1_qtl_end,
                               color = col_blue[[1]], labels = "", 
                               maxRadius = 1)

tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', 
                                           links_chr2_intron, links_chr2_intron_start, links_chr2_intron_end,
                                           links_chr2_qtl, links_chr2_qtl_start, links_chr2_qtl_end,
                                           color = col_red[[1]], labels = "", 
                                           maxRadius = 1)

tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', 
                                           links_chr3_intron, links_chr3_intron_start, links_chr3_intron_end,
                                           links_chr3_qtl, links_chr3_qtl_start, links_chr3_qtl_end,
                                           color = col_green[[1]], labels = "", 
                                           maxRadius = 1)


BioCircos(tracklist, genomeFillColor = c(rep("wheat2", 22)),
          chrPad = 0.02, displayGenomeBorder = FALSE, yChr =  FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = "12pt", 
          genomeLabelDy = 0, genome = as.list(lengthChr), width= '1000px', height= '700px', elementId = "wheat")


BioCircosOutput("wheat", width = "100%", height = "400px")

