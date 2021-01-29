library("dplyr")
library("stats")
library("bedr")
library("tidyr")
library("data.table")
library("ggplot2")

setwd("/path/to/epigenetic/mod/.bed/file")

#Before running this script, please check whether the current wd has the files

#hist_HiC_df <- bedr(engine = "bedtools",params = c("-wo"), input = list(x = "HiC_fullgen_ygi.bed", y = "GM12878_H3K27ac.bed"), method = "intersect")
#This will intersect the HiC data with the epigenetic bed file data
cat("** Using Bedtools (R-package) to form 1MB fragments of Epigenetic data ** \n \n")
hist_HiC_df <- bedr(engine = "bedtools", params=c("-wo"), input = list(a="HiC_all_chr_loci.bed", b=file.choose(new = FALSE)), method = "intersect") 

head(hist_HiC_df)
dim(hist_HiC_df)
hist_HiC_df$V11 <- as.numeric(hist_HiC_df$V11)
hist_HiC_df$V17 <- as.numeric(hist_HiC_df$V17)
#Calculating signal
cat("** Calculating Signal value for every 1MB fragment of Epigenetic data ** \n \n")
Signal <- hist_HiC_df$V11*hist_HiC_df$V17/1000000
hist_HiC_df <- mutate(hist_HiC_df,Signal)
#signal_sum <- aggregate(Signal~start+end+V4+V8+V9, data=hist_HiC_df, FUN=sum)
#signal_sum <- as.data.table(hist_HiC_df)[, sum(Signal), by = .(start, end)]
cat("** Summing up the signal value for every 1MB fragment ** \n \n")
ab <- group_by(hist_HiC_df,chr,start,end) 
signal_sum <- summarize(ab,sum_Signal = sum(Signal))
signal_sum_df <- as.data.frame(signal_sum)
head(signal_sum_df)
dim(signal_sum_df)
chr_start <- as.numeric(signal_sum_df$start)
length(chr_start)
chr_start <- format(chr_start + 1,scientific = FALSE)
signal_sum_df <- format(mutate(signal_sum_df,chr_start),scientific = FALSE)
drops <- c("start")
signal_sum_df <- signal_sum_df[ , !(names(signal_sum_df) %in% drops)]
signal_sum_df <- signal_sum_df %>% 
  rename(
    start = chr_start
  )
cat("** Generating final dataframe containing 1MB fragments with Epigenetic Signal value ** \n \n")
signal_sum_df <- signal_sum_df[c("chr","start","end","sum_Signal")]
head(signal_sum_df)
#signal_sum_df <- as.matrix(signal_sum_df)
#signal_sum_df <- signal_sum_df[, 2] + 1
#signal_sum_df <- signal_sum_df$start+1
cat("Processing for integrating ChIP Seq with HiC Seq: \n")
cat("A) Taking all permutations and combinations for comparisons between Chromosome 1-22. \n")
cat("B) Subsetting Epigenetics dataframe for extracting info of each Chromosome. \n")
cat("C) Changing the format in Epigenetics dataframe to match HiC data format. \n")
cat("D) Taking HiC matrices and extarcting rows and columns with corresponding Epigenetic data. \n")
cat("E) Calculating Epigenetic Signal diffrence between 2 interacting fragments, i.e. D(i,j) matrix. \n")
cat("F) Flattening the HiC matrix & the D(i,j) matrix into vector form. \n")
cat("G) Storing Vector forms of HiC & D(i,j) matrix for all Chromosome combinations. \n \n")
cat("Printing Chromosome Combinations: \n")

all_Dij <- c()
all_HiC <- c()
all_Dij_med <- c()
all_HiC_med <- c()
setwd("/path/to/HiC/matrix/files")
for (p in 1:21) {
  for (q in (p+1):22) {
    
    #Chromosome A- histone modification dataset
    chrNo=paste("chr",p,sep = "")
    # get subset of chrA histone modification data
    hist_mod_chromA=subset(signal_sum_df,signal_sum_df$chr==chrNo) 
    #convert sum_signal in hist_mod_chromA into a numeric form for use in calculating D(i.j) matrix
    #later
    hist_mod_chromA$sum_Signal <- as.numeric(hist_mod_chromA$sum_Signal)
    chr_A <- paste0(hist_mod_chromA$chr,":",hist_mod_chromA$start,"-",hist_mod_chromA$end)
    chr_A <- gsub(" ","",chr_A)
    #Chromosome B- histone modification dataset
    chrNo=paste("chr",q,sep = "")
    #get subset of chrB histone modification data
    hist_mod_chromB=subset(signal_sum_df,signal_sum_df$chr==chrNo)
    #convert sum_signal in hist_mod_chromB into a numeric form for use in calculating D(i.j) matrix
    #later
    hist_mod_chromB$sum_Signal <- as.numeric(hist_mod_chromB$sum_Signal)
    chr_B <- paste0(hist_mod_chromB$chr,".",hist_mod_chromB$start,".",hist_mod_chromB$end)
    chr_B <- gsub(" ","",chr_B)
    #ordering the peaks from smallest position along chromosome
    #hist_mod_chromB=hist_mod_chromB[order(hist_mod_chromB$start),]
    
    file <- paste0("HIC_gm06690_chr",p,"_","chr",q,"_","1000000_","pearson",".txt")
    whole_matrix <- readLines(file)
    ignore <- whole_matrix[-c(1:1)]
    chrAB_HiC <- read.csv(textConnection(ignore), sep="\t", stringsAsFactors = FALSE, row.names = 1)
    #chrAB_HiC %>% remove_rownames %>% column_to_rownames(var="HIC_bin1.hg18.chr22.1.999999")
    drops <- c("X.1")
    chrAB_HiC <- chrAB_HiC[ , !(names(chrAB_HiC) %in% drops)]
    chromoA <- rownames(chrAB_HiC)
    #chromoA <- gsub("/|",".",chromoA)
    chromoB <- names(chrAB_HiC)
    
    
    chrAB_HiC <- chrAB_HiC[grep(paste(chr_A,collapse="|"), chromoA, value = TRUE),]
    chrAB_HiC <- t(chrAB_HiC)  
    #chrAB_HiC <- chrAB_HiC[,which(chromoB %in% chr_B)]
    chrAB_HiC <- chrAB_HiC[grep(paste(chr_B,collapse="|"), chromoB, value = TRUE),]
    chrAB_HiC <- t(chrAB_HiC)
    chrAB_HiC <- as.matrix(chrAB_HiC)
    chrAB_HiC <- matrix(chrAB_HiC,ncol = ncol(chrAB_HiC),dimnames = NULL)
    
    #chrAB_HiC[grep("", rownames(chrAB_HiC)), ]
    #chrAB_HiC <- subset(chrAB_HiC)
    #D(i.j) matrix to find the difference in the signal from 2 interacting fragments
    #D(i,j) <- log10(hist_mod_chromA$sum_Signal/hist_mod_chromB$sum_Signal)
    
    Difference_mat <- matrix(0 , nrow = nrow(hist_mod_chromA), ncol = nrow(hist_mod_chromB))  #initialize matrix with 0's
    for (x in 1:nrow(hist_mod_chromA)) {
      for (y in 1:nrow(hist_mod_chromB)) {
        Difference_mat[x,y] <- abs(log10(hist_mod_chromA[x,"sum_Signal"]/hist_mod_chromB[y,"sum_Signal"]))
      }
    }
    
    #Flattening the HiC and D(i,j) matrices
    Dij_vec <- as.vector(Difference_mat)
    chrAB_HiC_vec <- as.vector(chrAB_HiC)
    #Dij_med <- median(Dij_vec)
    #HiC_med <- median(chrAB_HiC_vec)
    #all_Dij_med <- append(all_Dij_med,Dij_med)
    #all_HiC_med <- append(all_HiC_med,HiC_med)
    all_Dij <- append(all_Dij,Dij_vec)
    all_HiC <- append(all_HiC,chrAB_HiC_vec)
    
    print(paste0("Chromosome A: ", p))
    print(paste0("Chromosome B: ", q))
    
    
  }
  
}
cat("\n \n")
cat("Generating plots... \n \n")

minimum=-0.32
maximum=0.57
HiC_D=data.frame(HiC_distance=all_HiC,D=all_Dij)
HiC_D=subset(HiC_D,HiC_distance > minimum)
HiC_D=subset(HiC_D,HiC_distance < maximum)
intervals=18.0
intervalWidth=(maximum-minimum)/intervals
HiC_D=mutate(HiC_D,Grouped_HiC_distance=ceiling((HiC_distance-minimum)/intervalWidth))
intervalNames=seq((minimum+intervalWidth/2.0),(maximum-intervalWidth/2.0),intervalWidth)

boxplot(D~Grouped_HiC_distance,HiC_D,outline = FALSE, boxlty = 1,
        whisklty = 1, staplelty = 1,las=2,names=round(intervalNames,2),
        xlab = "Spatial proximity", ylab = "Histone modification distance D(i,j)",
        main="name_of_modification",cex.lab=1.3, cex.main=1.5)





