library(BioCircos)
library(RColorBrewer)
library(grDevices)
library(stringr)

genomeChr = c("1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Un")

lengthChr = c(594102056/1000000, 689851870/1000000, 495453186/1000000, 780798557/1000000, 801256715/1000000, 651852609/1000000, 
              750843639/1000000, 830829764/1000000, 615552423/1000000, 
              744588157/1000000, 673617499/1000000, 509857067/1000000, 
              709773743/1000000, 713149757/1000000, 566080677/1000000, 
              618079260/1000000, 720988478/1000000, 473592718/1000000, 
              736706236/1000000, 750620385/1000000, 638686055/1000000, 
              480980714/1000000)

names(lengthChr) <- genomeChr

file <- read.csv("jul3rd_Results/sQTL_Results.csv", header=TRUE)
dim(file)
head(file)
colSums(is.na(file[!is.na(file$qtl_start),]))
#file <- file[-which(is.na(file$qtl_start)),]
dim(file)

head(file)
file <- file[unique(file$clusterID.1),]
dim(file)

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

length(links_chr1_intron)
length(links_chr2_intron)
length(links_chr3_intron)
colors()
?rgb() # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
#col_blue <- rgb(50, 50, 200, max = 200, alpha = 125, names = "skyblue2")
#col_red <- rgb(200, 50, 50, max = 200, alpha = 125, names = "rosybrown")
#col_green <- rgb(100, 255, 150, max = 255, alpha = 125, names = "darkgreen")
library("yarrr")
col_green <- yarrr::transparent("darkgreen", trans.val = .7)
col_blue <- yarrr::transparent("skyblue2", trans.val = .7)
col_red <- yarrr::transparent("rosybrown", trans.val = .7)

col_red <- yarrr::transparent("#724E2F", trans.val = .7)
col_green <- yarrr::transparent("darkgreen", trans.val = .7)
col_blue <- yarrr::transparent("#D9870F", trans.val = .7)


rm(tracklist)

tracklist = BioCircosLinkTrack('myLinkTrack', 
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
          genomeLabelDy = 0, genome = as.list(lengthChr))



















+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
listcsv <- sort(dir(pattern = "*peaks.csv"))
listcsv <- listcsv[-1]

peaks <- list()
for (i in 1:length(listcsv)){
  peaks[[i]] <- read.csv(listcsv[i])
}


library(stringr)

for (i in 1:length(peaks)) {peaks[[i]]
  
  
}



lapply(peaks, function(x) str_split_fixed(x$lodcolumn, ":", 4))

colnames(file)[which(colnames(file) == "V9")] <- c("intron_chr")
colnames(file)[which(colnames(file) == "V10")] <- c("intron_start")
colnames(file)[which(colnames(file) == "V11")] <- c("intron_end")

colnames(file)[which(colnames(file) == "V12")] <- c("qtl_chr")
colnames(file)[which(colnames(file) == "V13")] <- c("qtl_start")
colnames(file)[which(colnames(file) == "V14")] <- c("qtl_end")
