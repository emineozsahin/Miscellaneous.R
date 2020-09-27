annot <- read.csv("sQTL_annotated.csv", header = FALSE)

annot_trans<- annot[which(annot$V15 == "trans"),]
unique_introns_trans<- unique(annot_trans$V3)

length(unique_introns_trans)

intron_qtl_num_trans <- data.frame()
for (i in 1:length(unique_introns_trans)){
  intron <- unique_introns_trans[i]
  subset <- annot_trans[which(annot_trans$V3 == intron),]  # subset with the gene take all the rows has similar gene
  num <- length(unique(subset$V20))  # count haw many unique sQTL region 
  temp <- data.frame(intron = intron, num_of_qtl = num) 
  intron_qtl_num_trans <- rbind(intron_qtl_num_trans, temp)
}

dim(intron_qtl_num_trans)
head(intron_qtl_num_trans)

#unique introns
introns_only_one_trans <- intron_qtl_num_trans$intron[which(intron_qtl_num_trans$num_of_qtl == 1)] 
length(introns_only_one_trans)
#Extract the introns which has one trans QTL
introns_only_one_trans_all <- annot_trans[unlist(lapply(introns_only_one_trans, function(x) which(annot_trans$V3 == x))),]
# Group them by their QTL effect 
introns_only_one_trans_AA <- introns_only_one_trans_all[which(introns_only_one_trans_all$V8 == 'AA'),]
introns_only_one_trans_BB <- introns_only_one_trans_all[which(introns_only_one_trans_all$V8 == 'BB'),]
introns_only_one_trans_AA <- as.character(introns_only_one_trans_AA$V3)
introns_only_one_trans_BB <- as.character(introns_only_one_trans_BB$V3)

length(introns_only_one_trans_AA)
length(introns_only_one_trans_BB)

######
annot <- read.csv("sQTL_annotated.csv", header = FALSE)
annot_cis <- annot[which(annot$V15 == 'cis'),]
unique_introns_cis<- unique(annot_cis$V3)

intron_qtl_num_cis <- data.frame()
for (i in 1:length(unique_introns_cis)){
  intron <- unique_introns_cis[i]
  subset <- annot_cis[which(annot_cis$V3 == intron),]  # subset with the gene take all the rows has similar gene
  num <- length(unique(subset$V20))  # count haw many unique sQTL region 
  temp <- data.frame(intron = intron, num_of_qtl = num) 
  intron_qtl_num_cis <- rbind(intron_qtl_num_cis, temp)
}

head(intron_qtl_num_cis)

#unique introns
introns_only_one_cis <- intron_qtl_num_cis$intron[which(intron_qtl_num_cis$num_of_qtl == 1)] 
#Extract the introns which has one trans QTL
introns_only_one_cis_all <- annot_cis[unlist(lapply(introns_only_one_cis, function(x) which(annot_cis$V3 == x))),]
# Group them by their QTL effect 
introns_only_one_cis_AA <- introns_only_one_cis_all[which(introns_only_one_cis_all$V8 == 'AA'),]
introns_only_one_cis_BB <- introns_only_one_cis_all[which(introns_only_one_cis_all$V8 == 'BB'),]
introns_only_one_cis_AA <- as.character(introns_only_one_cis_AA$V3)
introns_only_one_cis_BB <- as.character(introns_only_one_cis_BB$V3)


length(introns_only_one_cis_AA)
length(introns_only_one_cis_BB)


Use venndiagram codes found in miscellaneous.R file
library("VennDiagram")
venn.diagram (x = list(x1 = introns_only_one_cis_AA, 
                       x2 = introns_only_one_cis_BB,
                       x3 = introns_only_one_trans_AA, 
                       x4 = introns_only_one_trans_AA),
              category.names = c("Red Fife cis-sQTL", 
                                 "Stettler cis-sQTL",
                                 "Red Fife trans-sQTL",
                                 "Stettler trans-sQTL"),
              fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"), 
              main = "VennDiagram", main.pos = c(0.5, 1.05), main.just = c( 0.5, 1), 
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
              filename = "a.tiff",
              sub.pos = c( 0.5, 1.05), # bu lazim olmayabilir bir bak
              ext.pos = "", ext.percent = "", ext.line.lwd = "", ext.line.lty = "", description = "",  ext.dist = "", ext.length = "")














