#' Axel Rodriguez Perez
#' Final GWAS script
#' Edited from Chris Van Hout lab scripts
#' nov 27, 2023

# 1. READ IN THE DATA ----------------------------------------------------------
# Loading required libraries
library(RCurl)
library(downloader)
library(snpStats)
library(SNPRelate)
library(gdsfmt)
library(dplyr)
library(doParallel)

# Reading files
bedFile<-"/home/axel/Documents/ENES/5to/StatisticsModels/lab/GWAS_data.bed"
bimFile<-"/home/axel/Documents/ENES/5to/StatisticsModels/lab/GWAS_data.bim"
famFile<-"/home/axel/Documents/ENES/5to/StatisticsModels/lab/GWAS_data.fam"

data <- read.plink(bedFile,bimFile,famFile)

clinicalFile <- "/home/axel/Documents/ENES/5to/StatisticsModels/final_GWAS/Clinical_IBD.csv"
clinical <- read.csv(clinicalFile, colClasses = c("character", rep("factor",2), rep("numeric", 2)))

# 2. UNDERSTANDING THE DATA ----------------------------------------------------
#' Here, we use the PennCath cohort data, comming from a GWAS study of Coronary Artery 
#' Disease and cardiovascular risk factors based at the University of Pennsylvania 
#' Medical Center. 
#' De-identified data used in this project has a total n = 1401 individuals with 
#' genotype information across 500,000 single-nucleotide polymorphisms (SNPs). 
#' Corresponding clinical data, including age, sex, high-density lipoprotein (HDL)-cholesterol (HDL-C), 
#' low-density lipoprotein (LDL)-cholesterol, triglycerides (TGs) and CAD status 
#' are available as well. 
#' 
#' The BIM file contains information on each SNP; the BED file,contains every SNP in the study, 
#'  as well as the genotype at that snp for each individual; and the FAM file contains the participant identification information.
bim <- data$map
head(bim)

bed <- data$genotypes
head(bed)

fam <- data$fam
head(fam)

#' As we can observe, we have information of 1,401 individuals across 500,000 SNPs.
#' From those individuals, 937 are males and 464 are females (sum(fam$sex == 1) 
#' for males and sum(fam$sex == 2) for females).

## 2.1 The Clinical Data -------------------------------------------------------
#' The clinical file contains information on things like additional covariates 
#' for each individual and other “outcome variables”. The rows of this 
#' file represent each subject and the columns correspond to the available 
#' covariates and phenotypes.
rownames(clinical) <- clinical$FamID
print(head(clinical))
#' From this, we have 1,281 complete observations and can observe that the age 
#' of range of study participants is from 22 to 87 (summary(clinical$age)).
#' In this project, our trait of interest is Inflamatory Bowel Disease (IBD). 
#' We have 367 cases of positive IBD (sum(clinical$IBD_1Control_2Case == 2)).

# 3. PROCESSING OF GENOTYPE DATA -----------------------------------------------
genotype <- data$genotypes # SnpMatrix
genoMap <- data$map # Table with SNPs info
colnames(genoMap)<-c("chr", "SNP", "gen.dist", "position", "A1", "A2")

snpsum.col <- col.summary(genotype) # Summary statistics of each SNP
head(snpsum.col)

# 4. QUALITY CONTROL -----------------------------------------------------------
## 4.1. Call Rate
# For a given SNP, its call rate is defined as the proportion of inidividuals in 
# our study for which the corresponding SNP information is not missing. For 
# example, a call rate of 95% for a certain SNP means that 95% of the individuals 
# have data for this SNP.
call <- 0.95
keep <- with(snpsum.col, (!is.na(Call.rate) & Call.rate >= call))
head(keep)

keep[is.na(keep)] <- FALSE

cat(ncol(genotype)-sum(keep),"SNPs will be removed due to low call rate.\n")

genotype <- genotype[,keep]
snpsum.col <- snpsum.col[keep,]

print(genotype)
summary(genotype)
# After this, we have 445417 remaining SNPs

## 4.2 Minor allele frequency (MAF) --------------------------------------------
# Inadequate power to infer a statiscally significant relationship between the 
# SNP and the trait under study is the result of a large degree of homogeneity at
# a given SNP across study participants. This occurs when we have a very small 
# minor allele frequency (MAF), which means the majority of the individuals have
# two copies of the same major alleles. We remove SNPs whose minor allele fequency is 
# less than 1%.
minor <- 0.01

keep1 <- with(snpsum.col, (!is.na(MAF) & MAF > minor) )
keep1[is.na(keep1)] <- FALSE

cat(ncol(genotype)-sum(keep1),"SNPs will be removed due to low MAF .\n"  )

# After knowing the SNPs needed to be removed, we subset the genotype and SNP 
# summary data for the SNPs that pass MAF criteria.

genotype <- genotype[,keep1] 
snpsum.col <- snpsum.col[keep1,] 

snpsum.row <- row.summary(genotype) # <-- summary statistics per individuals
head(snpsum.row)

snpsum.col<-col.summary(genotype)  # per SNP
head(snpsum.col)

MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype) # Removing NAs from genotype info

## 4.3 Heterozygosity ----------------------------------------------------------
# Heterozygosity refers to the presence of each of two alleles at a given SNP 
# within an individual. This is expected under Hardy-Weinberg Equilibrium (HWE) 
# to occur with probability of 2p(1-p) where p is the dominant allele frequency 
# at that SNP (assuming a bi-allelic SNP). Excess or deficient heterozygosity 
# across typed SNPs within an individual is commonly considered as an indication 
# of inbreeding between relatives, but also serves as a screen for sample 
# contamination or poor sample quality.

hetExp <- callmatrix %*% (2*MAF*(1-MAF)) #Note that %*% is the 'multiply_by_matrix' function from plyr
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

par(mfrow=c(1,2))
# Plot of heterozygosity
plot(hetObs, hetExp,
     xlab = "hetExp", ylab = "hetOBS",
     main = "Heterozygosity",
     col = rgb(0, 0, 0, 0.5), pch = 16, cex = 1.5)
abline(0,1, col = "red", lty = 2)
grid()

# Inbreeding coefficient histogram
hist(snpsum.row$hetF, breaks = 30, col = rgb(0.2, 0.5, 0.8, 0.7), border = "white", 
     main = "Ibreeding coefficient", xlab = "Value", ylab = "Frequency")
lines(density(snpsum.row$hetF), col = "darkblue", lwd = 2) # density curve
abline(v = mean(snpsum.row$hetF), col = "red", lwd = 2, lty = 2) # vertical line for the mean

head(snpsum.row)

hetcutoff <- 0.055

sampleuse <- with(snpsum.row, abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE  
cat(nrow(genotype)-sum(sampleuse), "individuals will be removed due to low sample call rate or unusual inbreeding coefficient.\n")

# After identifying the individuals to be removed, we subset the genotype and 
# clinical data set using information obtained above to remove the individuals 
# that do not pass the criteria from our sample.

genotype <- genotype[sampleuse,]
clinical<- clinical[rownames(genotype), ]
snpsum.row<-row.summary(genotype)

## 4.4 Call Rate ---------------------------------------------------------------
# To filter by call rate, we remove individuals that are missing genotype data 
# across more than a pre-defined percentage of the typed SNPs. This proportion of 
# missingness across SNPs is the sample call rate and we apply a threshold of 5%, 
# meaning that we drop those individuals missing genotype data for more than 5% 
# of the typed SNPs. A new reduced dimesion SnpMatrix gentotype is created from this filter.

sampcall <- 0.95 

sampleuse1 <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall)
sampleuse1[is.na(sampleuse1)] <- FALSE 
cat(nrow(genotype)-sum(sampleuse1), "individuals will be removed due to low sample call rate.\n")

genotype <- genotype[sampleuse1,]
clinical<- clinical[ rownames(genotype), ]
snpsum.row<-row.summary(genotype)

## 4.5 Hardy-Weinberg Equilibrium (HWE) ----------------------------------------
#' HWE departure is common signal across SNPs with genotyping errors. Deviations 
#' from HWE within individuals with disease (cases) can be due to association 
#' between genotypes and disease status.
hardy <- 10^-6
controls <- as.character(clinical[ clinical$IBD_1Control_2Case == 1, 'FamID' ]) 
snpsum.colControls <- col.summary(genotype[controls,])
HWEuse <- with(snpsum.colControls, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colControls)
HWEuse[is.na(HWEuse)] <- FALSE          
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n") 
genotype <- genotype[,HWEuse]
print(genotype)

## 4.6 Cryptic relatedness, duplicates and gender identity ---------------------
# Population-based cohort studies often only concerns unrelated individuals and 
# the generalized linear modeling approach assumes independence across individuals. 
# In regional cohort studies (e.g. hospital-based cohort studies) of complex 
# diseases, however, individuals from the same family could be unintentionally 
# recruited. To mitigate this problem, we thus employ Identity-by-descent (IBD) 
# analysis which is a common measurement of relatedness and duplication. We set 
# LD threshold to be 20% and 
# kinship threshold to be 5%, which means that we would remove any SNPs whose 
# relevant statistics are higher than our stated threshold values.
ld.thresh <- 0.2
kin.thresh <- 0.05

snpgdsBED2GDS(bedFile, famFile, bimFile,"GWAS_data.gds") 

genofile <- snpgdsOpen("GWAS_data.gds", readonly = FALSE)
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          
                          sample.id = geno.sample.ids, 
                          
                          snp.id = colnames(genotype)) 

snpset.ibd <- unlist(snpSUB, use.names=FALSE)

cat(length(snpset.ibd),"will be used in IBD analysis\n")

ibd <- snpgdsIBDMoM(genofile, kinship=TRUE, sample.id = geno.sample.ids, snp.id = snpset.ibd,  num.thread = 1)

ibdcoeff <- snpgdsIBDSelection(ibd)    
head(ibdcoeff)

kins<-subset(ibdcoeff,ibdcoeff$kinship>0.05) #for context, the kinship coefficient for 1st cousins is 0.0625

ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]
related.samples <- NULL

while ( nrow(ibdcoeff) > 0 ) {
  sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing sample", as.character(rm.sample), 'too closely related to', 
      sample.counts[1, 'freq'],'other samples.\n')
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n") 
print(genotype)
### 4.6.1 IBD Plot -------------------------------------------------------------
par(mfrow=c(1,1))
plot(kins$k0, kins$k1, 
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "IBD0", ylab = "IBD1",
     main = "IBD Analysis",
     col = rgb(0, 0, 0, 0.5), pch = 16, cex = 1.5)
abline(h = 0.5, v = 0.5, col = "red", lty = 2)
grid()

## 4.7 PCA ---------------------------------------------------------------------
# Principle components (PCs) from multi-dimensional scaling (MDS) is one approach 
# to visualizing and classifying #individuals into ancestry groups based on their 
# observed genetics make-up.
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=2)

pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10], stringsAsFactors = FALSE)

colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))

pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    
                    PC2 = pca$eigenvect[,2],
                    PC3 = pca$eigenvect[,3],
                    stringsAsFactors = FALSE)
par(mfrow=c(2,2))
# PCA 1 vs PCA 2
plot(pctab$PC2, pctab$PC1, 
     xlab = "Principal Component 2", ylab = "Principal Component 1", 
     main = "Ancestry Plot PCA1 vs PCA2", pch = 16, col = rgb(0, 0, 0, 0.5))
grid()
axis(side = 1, at = seq(-0.2, 0.2, 0.1), labels = seq(-0.2, 0.2, 0.1))
axis(side = 2, at = seq(-0.2, 0.2, 0.1), labels = seq(-0.2, 0.2, 0.1))

# PCA 1 vs PCA 3
plot(pctab$PC3, pctab$PC1, 
     xlab = "Principal Component 3", ylab = "Principal Component 1", 
     main = "Ancestry Plot PCA1 vs PCA3", pch = 16, col = rgb(0, 0, 0, 0.5))
grid()
axis(side = 1, at = seq(-0.2, 0.2, 0.1), labels = seq(-0.2, 0.2, 0.1))
axis(side = 2, at = seq(-0.2, 0.2, 0.1), labels = seq(-0.2, 0.2, 0.1))

# PCA 2 vs PCA 3
plot(pctab$PC3, pctab$PC2, 
     xlab = "Principal Component 3", ylab = "Principal Component 2", 
     main = "Ancestry Plot PCA2 vs PCA3", pch = 16, col = rgb(0, 0, 0, 0.5))
grid()
axis(side = 1, at = seq(-0.2, 0.2, 0.1), labels = seq(-0.2, 0.2, 0.1))
axis(side = 2, at = seq(-0.2, 0.2, 0.1), labels = seq(-0.2, 0.2, 0.1))

plot(pca$eigenval,ylab="Eigenvalue",xlab="Index of Eigenvector",xlim=c(1,20),pch=20,main="Scree plot")
# sum(pca$eigenval[1:3])

# since we did not calculate all PCs, just the first few, we can't calculate the 
# proportion of variance captured by each PC, which is the eigenvalue / trace, 
# where the trace is the sum of all eigenvalues. Instead, we simply plot the eigenvalue.

closefn.gds(genofile)

# 5. RUNNING  A GWAS ------------------------------------------------------------
genodata <- genotype

# Print the number of SNPs to be checked
cat(paste(ncol(genodata), "Thus far, SNPs are included in analysis.\n"))


# create text file for GWAA output to be written to
columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")

# the output file handle is written in the mac style, with '~' referring to the current users home directory. 
write.table(t(columns), "~/Documents/ENES/5to/StatisticsModels/final_GWAS/GWAA.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

## 5.1 Manipulating Phenotype Data ---------------------------------------------
traitName<-'IBD_1Control_2Case' # The trait
clinicalSub<-clinical[,c('FamID','sex','age',traitName)]
clinicalSub <- clinicalSub %>% rename(prePhenotype = traitName)
clinicalSub$phenotype <- ifelse(clinicalSub$prePhenotype == 1, 0, ifelse(clinicalSub$prePhenotype == 2, 1, clinicalSub$prePhenotype))  # Originally controls = 1 & cases = 2, now controls = 0 & cases = 1
clinicalSub$prePhenotype <- NULL
pcsSub<-pcs[,c('FamID','pc1','pc2')] #add additional PCs here and merge in the next step

phenodata <- merge(clinicalSub,pcsSub,by="FamID")     
phenodata<-phenodata[!is.na(phenodata$phenotype),]
genodata <- genodata[as.character(phenodata$FamID),]

cat(nrow(genodata), "samples included in analysis thus far.\n")

## 5.2 The GWAA function
# The function itself has four major elements. Given the data described above, 
# it (1) determines how many cores to use/which method to employ for parallel 
# processing (depending on specified cores and operating system), (2) converts 
# the genotype data into the necessary format for an additive model, (3) fits a 
# linear model for each SNP and covariates (one group at a time) and then (4) 
# writes the SNP coefficient details (the estimate, standard error, test statistic 
# and p value) from each model out to the GWAA.txt file. Read through the function 
# to get a general idea of the format in regards to these four steps. Further 
# details regarding the function elements are provided below.

GWAA <- function(genodata, phenodata, filename = NULL, append = FALSE, workers = getOption("mc.cores", 2L), flip = FALSE, select.snps = NULL, hosts = (detectCores()/2)) {
  if (is.null(hosts)) {
    cl <- makeCluster(workers)
  } else {
    cl <- makeCluster(hosts, "PSOCK")
  }
  show(cl)
  registerDoParallel(cl)
  # Function that will change which allele is counted (major or minor)
  flip.matrix <- function(x) {
    zero2 <- which(x == 0)
    two0 <- which(x == 2)
    x[zero2] <- 2
    x[two0] <- 0
    return(x)
  }
  foreach(part = 1:nSplits) %do% {
    genoNum <- as(genodata[, snp.start[part]:snp.stop[part]], "numeric")
    # flip.matrix function employed
    if (isTRUE(flip)) 
      genoNum <- flip.matrix(genoNum)
    rsVec <- colnames(genoNum)
    res <- foreach(snp.name = rsVec, .combine = "rbind") %dopar% {
      result <- summary(glm(phenotype ~ . - FamID, family = binomial, data = cbind(phenodata, snp = genoNum[, snp.name])))  #GLM family must match the continuously distributed trait family=gaussian, or a binary trait family=binomial
      result$coefficients["snp", ]
    } 
    write.table(cbind(rsVec, res), filename, append = TRUE, quote = FALSE, 
                
                col.names = FALSE, row.names = FALSE)
    cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 
                100 * part/nSplits))
  }
  stopCluster(cl)
  return(print("Done."))
} #end GWAA

nSNPs <- ncol(genodata)
nSplits <- 20
genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

start <- Sys.time()
GWAA(genodata, phenodata, filename="~/Documents/ENES/5to/StatisticsModels/final_GWAS/GWAA.txt") #note this file has 5 columns; SNP,Estimate,Std.Error,t-value,Pr(>|t|) 
end <- Sys.time()
print(end-start) #output approximately how long the analysis took to compute

gwas<-read.table("~/Documents/ENES/5to/StatisticsModels/final_GWAS/GWAA.txt",sep=" ",head=T)
# Annotate variant information for any SNP that exists in the gwas output. 
# left_join annotates any SNP in the gwas data.frame with information from genoMap. 
# right_join annotates any SNP in the genoMap data.frame with information from gwas. 
gwasAnnotation<- gwas %>% left_join(genoMap, by = "SNP")

# 6. PLOTS ---------------------------------------------------------------------
## 6.1 Manhattan Plot
par(mfrow=c(1,1))
library('qqman')
if (exists("genoMap")) {
  gwas<-read.table("~/Documents/ENES/5to/StatisticsModels/final_GWAS/GWAA.txt",sep=" ",head=T)
  gwasAnnotated<- gwas %>% left_join(genoMap, by = "SNP") #Annotate variant 
  # information for any SNP that exists in the gwas output. left_join annotates 
  # any SNP in the gwas data.frame with information from genoMap. right_join 
  # annotates any SNP in the genoMap data.frame with information from gwas. 
  # note that base R 'merge' threw an error that the vectors in the merge were 
  # too long, which they shouldn't be, but this works as expected
  
  significance_level <- 0.05/nSNPs # Bonferroni correction
  snps_of_interest <- gwasAnnotated$SNP[gwasAnnotated$p.value < significance_level]
  snps_of_interest_extended <- gwasAnnotated$SNP[gwasAnnotated$p.value < 1.0e-05]
  manhattan(gwasAnnotated,chr="chr",bp="position",p="p.value",snp="SNP", highlight = snps_of_interest, annotatePval = 1.0e-05, annotateTop = FALSE, col = c("blue", "orange"))
} else {
  print("Quitting. Have you run Lab.3 or loaded lab3GWAS_save.RData? Appropriate datasets do not exist.")
}

## 6.2 QQ Plot -----------------------------------------------------------------
# Calculation of lambdaGC (genomic control)
p_values <- gwasAnnotated$p.value
# Calculate observed lambda (genomic control)
lambdaGC <- qchisq(1-median(p_values),1)/qchisq(.5,1)
cat("Lambda (Genomic Control):", lambdaGC, "\n")

### 6.2.1 QQ Plot Function -----------------------------------------------------------
#' QQ plot From Hoffman et al Bioinformatics 2013
#' QQ plot and lambda_GC optimized for large datasets.
#' 
#' EXAMPLE 
#' p = runif(1e6)
#' QQ_plot(p)
#' 
#' @param p_values vector, matrix or list of p-values
#' @param col colors corresponding to the number of columns in matrix, or entries in the list
#' @param main title
#' @param pch pch
#' @param errors show 95\% confidence interval
#' @param lambda calculate and show genomic control lambda.  Lambda_GC is calcualted using the 'median' method on p-values > p_thresh.
#' @param p_thresh Lambda_GC is calcualted using the 'median' method on p-values > p_thresh.
#' @param showNames show column names or list keys in the legend
#' @param ylim ylim 
#' @param xlim xlim
#' @param plot make a plot.  If FALSE, returns lamda_GC values without making plot
#' @param new make a new plot.  If FALSE, overlays QQ over current plot
#' @param box.lty box line type
#' @param collapse combine entries in matrix or list into a single vector
#' @param ... other arguments
#' 
#' 
#' # get lambda_GC values without making plot
#' lambda = QQ_plot(p, plot=FALSE)
#'
#' @export 
QQ_plot = function(p_values, col=rainbow(min(length(p_values), ncol(p_values))), main="", pch=20, errors=TRUE, lambda=TRUE, p_thresh = 1e-7, showNames=FALSE, ylim=NULL, xlim=NULL, plot=TRUE,new=TRUE, box.lty=par("lty"), collapse=FALSE,...){
  
  if( collapse ){
    p_values = as.vector(unlist(p_values, use.names=FALSE))
  }
  
  # convert array, vector or matrix into list
  if( ! is.list(p_values) ){
    names(p_values) = c()
    keys = colnames(p_values)
    
    # if there is am empty name
    if( "" %in% keys ){
      keys[which("" %in% keys)] = "NULL"
    }
    p_values = as.matrix(p_values)
    p_values_list = list()
    
    for(i in 1:ncol(p_values) ){
      p_values_list[[i]] = p_values[,i]
    }
    names(p_values_list) = keys
    p_values = p_values_list
    rm( p_values_list )
  }
  
  rge = range(p_values, na.rm=TRUE)
  
  if( rge[1] < 0 || rge[2] > 1 ){
    stop("p_values outside of range [0,1]")
  }
  
  # assign names to list entries if they don't exist
  if( is.null( names( p_values ) ) ){
    names( p_values ) = 1:length( p_values )
  }
  
  # set pch values if not defined
  if( is.null( pch ) ){
    pch = rep(20, length(p_values))
  }
  
  if( length(pch) == 1){
    pch = rep(pch, length(p_values))
  }
  
  p_values = as.list(p_values)
  # Set the x and y ranges of the plot to the largest such values 
  # encountered in the data 
  ry = 0; rx = 0
  
  for( key in names( p_values ) ){
    # remove NA values
    p_values[[key]] = p_values[[key]][which( ! is.na( p_values[[key]] ))]
    j = which(p_values[[key]] == 0)
    if( length(j) > 0){
      p_values[[key]][j] = min(p_values[[key]][-j])
    }
    #ry = max(ry,range( -log10(p_values[[key]]), finite=TRUE) )
    ry = max(ry, -log10(min(p_values[[key]])) )
    rx = max(rx, -log10( 1/(length(p_values[[key]])+1)  )) 
  }
  
  if(!is.null(ylim)){
    ry = max(ylim)
  }
  
  if(!is.null(xlim)){
    rx = max(xlim)
  }
  
  r = max(rx, ry)
  xlab = expression(-log[10]('expected p-value'))
  ylab = expression(-log[10]('observed p-value'))
  
  if( plot && new ){
    # make initial plot with proper ranges
    plot(1, type='n', las=1, pch=20, xlim=c(0, rx), ylim=c(0, ry), xlab=xlab, ylab=ylab, main=main,...)
    abline(0, 1)
  }
  
  lambda_values = c()
  se = c()
  
  # Plots points for each p-value set
  # Since 90% of p-values should be < .1 and 99% < 0.01, plotting all of these points is
  # time consuming and takes up space, even though all the points appear on top of each other
  # Therefore, thin the p-values at the beginning and increases density for smaller p-values
  # This allows saving QQ plots as a PDF.....with out thinning, a PDF can be VERY large
  i = 1
  for( key in names( p_values ) ){
    observed = sort(p_values[[key]], decreasing=TRUE)
    # Calculated lambda_GC directly from sorted observed values
    # This is MUCH faster than using estlambda() from GenABEL
    # Compute only for p-values < 1e-6
    obs = sort(p_values[[key]][which(p_values[[key]] > p_thresh)], decreasing=TRUE)
    
    lambda_values[i] = qchisq(median( obs), 1, lower.tail = FALSE) / qchisq(0.5, 1)
    # option for a regression based approach, making sure that the observed pvalue distribution has a limit on the minimum
    
    if( plot ){
      # if a p-value is exactly zero, set it equal to the smallest nonzero p-value
      j = which( observed == 0)
      
      if( length(j) > 0){
        observed[j] = min(observed[-j])
      }
      
      expected = (length(observed):1) / (length(observed)+1) 
      
      p = length(expected)
      
      if( p < 1e6 ){
        # Thin p-values near 1, and increase density for smaller p-values 
        intervals = ceiling(c(1, 0.2*p, 0.4*p, 0.7*p, 0.9*p, 0.95*p, 0.99*p, p))
        scaling = c(28, 200, 800, 1000, 3000, 5000, p)
      }else{
        intervals = ceiling(c(1, 0.2*p, 0.4*p, 0.7*p, 0.9*p, 0.95*p, 0.99*p, 0.995*p, 0.999*p, p))
        scaling = c(28, 200, 800, 1000, 3000, 5000, 10000, 50000, p)
      }
      
      for( j in 1:length(scaling) ){
        k = seq(intervals[j], intervals[j+1], intervals[j+1]/scaling[j])
        points(-log10(expected[k]), -log10(observed[k]),  col=col[i], pch=pch[i])#,...)
      }
    }
    i = i + 1
  }
  
  # Plot lambda values, if desired
  if( plot && lambda ){
    # Calculated lambda using estlambda() from GenABLE
    # As of Sept24, 2013, I calculated it much faster using the sorted observed p-values
    
    #result = sapply( p_values, estlambda, plot=FALSE,method=method) # Better call GenABEL before uncommenting this
    #result = matrix(unlist(result), ncol=2, byrow=TRUE)
    #lambda_values = result[,1]
    #se = result[,2]
    
    if( ! showNames ){
      namesStrings = rep('', length(names( p_values )) )
    }else{
      namesStrings = paste( names( p_values ), ": ", sep='')
    }
    
    legend("topleft", legend = paste(namesStrings, format(lambda_values, digits = 4, nsmall = 3)), col=col, pch=15, pt.cex=1.5, title=expression(lambda['GC']), box.lty=box.lty)
  }
  
  # Standard errors 
  # From Supplementary information from Listgarten et al. 2010. PNAS 
  if( plot && errors ){
    error_quantile = 0.95
    # Using M:1 plots too many points near zero that are not discernable 
    # Reduce the number of points near zero 
    #plot(1, type='n', las=1, pch=20, xlim=c(0, rx), ylim=c(0, ry))
    M = length(p_values[[key]])
    #alpha = M:1
    alpha = seq(M, M/10 + 1, length.out=1000)
    
    if( M/10 > 1) alpha = append(alpha, seq(M/10, M/100 + 1, length.out=1000))
    if( M/100 > 1) alpha = append(alpha, seq(M/100, M/1000 + 1, length.out=1000))
    if( M/1000 > 1) alpha = append(alpha, seq(M/1000, M/10000 + 1, length.out=10000))
    alpha = append(alpha, seq(min(alpha), 1, length.out=10000))
    
    alpha = round(alpha)
    beta = M - alpha + 1
    
    x_top = qbeta(error_quantile, alpha, beta)
    x_bot = qbeta(1-error_quantile, alpha, beta)
    
    lines( -log10(alpha/(M+1)), -log10(x_top))
    lines( -log10(alpha/(M+1)), -log10(x_bot))
  }
  
  return( invisible(lambda_values) )
}



### 6.2.3 QQplot ---------------------------------------------------------------
QQ_plot(p_values)
