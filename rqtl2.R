# Emine Ozsahin June 2020
library("qtl2")
library("stringr")
library("tibble")

args <- commandArgs(trailingOnly = TRUE)
basename <- args[1]

#Example file name: QTL_perind.counts_chrname.gz.qqnorm_chr1.gz

# Prepare the phenotype file
# 3- pheno: wheat_pheno.csv ----

chr <- read.csv(basename, sep = "\t", stringsAsFactors = FALSE)

colnames(chr) <- colnames(chr) %>%
  str_replace_all(c("_05f2_sort4ReadGroupMarkDuplUniqMap.bam" = "",
                    "merged" = "",
                    ".bam" = "",
                    "ID" = "id")) 

chr$id <- chr$id %>%
  str_replace_all(c("_NA" = "")) 

cnames <- unlist(chr$id)

chr <- as.data.frame(t(chr[,-c(1,2,3)]))

colnames(chr) <- cnames

chr <- chr[-1, ]

rownames(chr) <- gsub("LL2017_", "", rownames(chr))
rownames(chr) <- gsub("L2017_", "", rownames(chr))

## rearrange the ind names

geno <- read.csv("wheat_geno.csv")

chr <- chr[as.character(unlist(geno$id)),]

chr <- rownames_to_column(chr, var = "id")

basename <- gsub("_perind.counts_chrname.gz", "", basename)
basename <- gsub(".gz", "", basename)
write.csv(chr, paste("wheat_", basename, "_pheno.csv", sep = ""), row.names = FALSE)

#write control file

write_control_file(paste("wheat_", basename, "_pheno.yaml", sep = ""), 
	comments=c("Emine Ozsahin June 2020","Wheat Alternative Splicing"), 
	crosstype="dh", 
	geno_file="wheat_geno.csv", 
	gmap_file="wheat_gmap.csv", 
	pheno_file=paste("wheat_", basename, "_pheno.csv", sep = ""), 
	geno_codes=c(A=1L, B=2L), 
	alleles=c("A","B"),
	na.strings=c("-","NA"), overwrite=TRUE)

## qtl  ---

#Calculating genotype probabilities
wheat <- read_cross2(paste("wheat_", basename, "_pheno.yaml", sep = ""))

map <- insert_pseudomarkers(map=wheat$gmap, step=1)
print ("map is done")

pr <- calc_genoprob(cross=wheat, map=map, error_prob=0.002)

#Performing a genome scan
out <- scan1(genoprobs = pr,  pheno = wheat$pheno, cores=8)
print("Scan is done.")

#Performing a permutation test ----
print ("permutation is starting...")
operm <- scan1perm(genoprobs = pr, pheno = wheat$pheno, n_perm = 1000, cores = 8)
print ("permutation is done.")

thr = summary(operm)

# Finding LOD peaks ----

peaks <- find_peaks(scan1_output = out, map = map, threshold = thr, prob = 0.95, expand2markers = FALSE, cores = 8)

write.csv(peaks[order(peaks[,5],decreasing = TRUE),], paste("wheat_", basename, "_peaks.csv", sep = ""))

high <- lapply(peaks$lodcolumn, function(x) which(colnames(out) == x))

ymx <- maxlod(out)

png(paste0(basename, "qtl_effect.png"), units="in", width=5, height=5, res=300)
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
lapply(high, function(x) plot(out, map, lodcolumn=x, ylim=c(0, ymx*1.02), col="slateblue", add=TRUE))
dev.off()

save(wheat, map, pr, out, operm, thr, peaks, file = paste(basename, ".RData", sep=""))
 
save.image()
