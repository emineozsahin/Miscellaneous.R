#Emine Ozsahin
# July 2020

library(qtl2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

#Example file name: QTL.qqnorm_chr1.RData

load(filename)

peaks <- peaks[order(peaks[,5],decreasing = TRUE),]

peaks <- cbind(peaks, unlist(str_split_fixed(peaks$lodcolumn, ":",4)))

peaks$`1` <- ifelse(peaks$`1`== "21", "7D",
                    ifelse(peaks$`1`== "20", "6D",
                           ifelse(peaks$`1`== "19", "5D",
                                  ifelse(peaks$`1`== "18", "4D",
                                         ifelse(peaks$`1`== "17", "3D",
                                                ifelse(peaks$`1`== "16", "2D",
                                                       ifelse(peaks$`1`== "15", "1D",
                                                              ifelse(peaks$`1`== "14", "7B",
                                                                     ifelse(peaks$`1`== "13", "6B",
                                                                            ifelse(peaks$`1`== "12", "5B",
                                                                                   ifelse(peaks$`1`== "11", "4B",
                                                                                          ifelse(peaks$`1`== "10", "3B",
                                                                                                 ifelse(peaks$`1`== "9", "2B",
                                                                                                        ifelse(peaks$`1`== "8", "1B",
                                                                                                               ifelse(peaks$`1`== "7", "7A",
                                                                                                                      ifelse(peaks$`1`== "6", "6A",
                                                                                                                             ifelse(peaks$`1`== "5", "5A",
                                                                                                                                    ifelse(peaks$`1`== "4", "4A",
                                                                                                                                           ifelse(peaks$`1`== "3", "3A",
                                                                                                                                                  ifelse(peaks$`1`== "2", "2A",
                                                                                                                                                         ifelse(peaks$`1`== "1", "1A",
                                                                                                                                                                ifelse(peaks$`1`== "22", "Un", NA))))))))))))))))))))))
















cat("finding qtl effects...")
qtl_effects <- data.frame()
for (i in 1:dim(peaks)[1]) { g <- as.data.frame(maxmarg(pr, map, chr=peaks$chr[i], pos=peaks$pos[i], return_char=TRUE))
colnames(g) <- "genotype"
c2eff <- scan1coef(pr[,peaks$`1`[i]], wheat$pheno[,i])
#g <- cbind(g, phenotype)
#g <- g[!is.na(g$genotype),]
#A <- g[g == "AA", ]
#B <- g[g == "BB", ]
#Amean <- mean(A$phenotype)
#Bmean <- mean(B$phenotype)
Amean <- mean(c2eff[,1])
Bmean <- mean(c2eff[,2])
#Alen <- length(A$genotype)
#Blen <- length(B$genotype)
chi_sq <- chisq.test(c2eff[,1],c2eff[,2], simulate.p.value = TRUE)
e <- if(Amean == "NaN" | Bmean == "NaN") "NA" else if(Amean > Bmean) "AA" else "BB"
temp_gen <- data.frame(intron=peaks$lodcolumn[i], AA=Alen, BB=Blen, A.pheno.mean=Amean, B.pheno.mean=Bmean, effect=e)
qtl_effects <- rbind(qtl_effects, temp_gen)}


cat("qtl effects are generated..", "\n")
map <- read.csv("SNP-location2.csv")

map$Chr <- gsub("chr", "", map$Chr)

#finding the qtl type
cat("finding qtl positions on phsical map...", "\n")

colnames(peaks)[8:11] <- c("intron_chr", "intron_start", "intron_end", "cluster")

qtl_pos_on_physical_map <- data.frame(stringsAsFactors = FALSE)
for (i in 1:length(peaks$chr)) {
  pos <- map[map$Marker.position..CM..on.genetic.map > peaks$ci_lo[i]  & map$Marker.position..CM..on.genetic.map < peaks$ci_hi[i] & map$Chr == peaks$chr[i],]
  pos <- pos[order(pos$Marker.position.on.physical.map),]
  pos <- pos[(!is.na(pos$Marker.position..CM..on.genetic.map)),]
  temp <- data.frame(qtl_pos_on_physical_map=paste0(pos$Marker.position.on.physical.map[1], ":", tail(pos, n=1)$Marker.position.on.physical.map))
  qtl_pos_on_physical_map <- rbind(qtl_pos_on_physical_map, temp)
}

f <- cbind(qtl_effects, peaks[,8:11], qtl_pos_on_physical_map, peaks$chr)

f <- cbind(f, unlist(str_split_fixed(f$qtl_pos_on_physical_map, ":",2)))

colnames(f)[12:14] <- c("qtl_chr", "qtl_start", "qtl_end") 

as.num = function(x, na.strings = "NA") {
  stopifnot(is.character(x))
  na = x %in% na.strings
  x[na] = 0
  x = as.numeric(x)
  x[na] = NA_real_
  x
}

f$intron_start <- as.numeric(as.character(f$intron_start))
f$intron_end <- as.numeric(as.character(f$intron_end))

f$qtl_start <- as.num(as.character(f$qtl_start, na.strings = "NA"))
f$qtl_start <- f$qtl_start-1000000

f$qtl_end <- as.num(as.character(f$qtl_end, na.strings = "NA"))
f$qtl_end <- f$qtl_end+1000000

f$t <- 1:length(f$intron)
for (i in 1:length(f$intron)) { f$t[i] <- if(!is.na(f$qtl_end[i]) & !is.na(f$qtl_end[i]) & f$qtl_chr[i] == f$intron_chr[i] & f$intron_start[i] > f$qtl_start[i] & f$intron_start[i] < f$qtl_end[i] & f$intron_end[i] > f$qtl_start[i] & f$intron_end[i] < f$qtl_end[i] )  "cis" else "trans" }

f$t <- as.character(f$t)
f$qtl_type <- f$t
f$qtl_type <- as.character(f$qtl_type)
for (i in 1:length(f$t)) {f$qtl_type[i] <- if(is.na(f$qtl_end[i])) f$qtl_end[i] else f$t[i]}

f <- f[,-which(colnames(f) == "t")]

f <- cbind(f, peaks$lod)
colnames(f)[which(colnames(f) =="peaks$lod")] <- "lod"

write.csv(f[,-c(10,11)], paste0("chr", peaks$`1`, "_qtl_effects_corrected_p_Val.csv"))

#save(c2eff, file = paste0(filename, "c2eff.RData"))
#save.image()
