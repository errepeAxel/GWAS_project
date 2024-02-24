rm(list=ls()) # remove any data objects or functions to avoid collisions or misspecification
getwd()
load("lab1_MyGWAS_save.RData")
setwd("~/Documents/ENES/5to/StatisticsModels/final_GWAS/")

#BiocManager::install("snpStats",force=TRUE) #run only once to install
#BiocManager::install("SNPRelate",force=TRUE) #run only once to install

library(snpStats)
library(SNPRelate)
library(gdsfmt)
library(plyr)

snpsum.row <- row.summary(genotype) # <-- por individuos
head(snpsum.row)
##       Call.rate    Certain.calls Heterozygosity
## 10002 0.9826528             1      0.3281344
## 10004 0.9891708             1      0.3231777
## 10005 0.9918277             1      0.3239852
## 10007 0.9860113             1      0.3241203
## 10008 0.9824172             1      0.3231505
## 10009 0.9912570             1      0.3205117

snpsum.col<-col.summary(genotype)  # por SNP
head(snpsum.col)
##            Calls Call.rate Certain.calls       RAF        MAF         P.AA       P.AB      P.BB      z.HWE
## rs4579145   1400 0.9992862             1 0.8235714 0.17642857 0.0285714286 0.29571429 0.6757143  0.6580953
## rs10125738  1401 1.0000000             1 0.9768023 0.02319772 0.0007137759 0.04496788 0.9543183 -0.2901317
## rs888263    1382 0.9864383             1 0.6649783 0.33502171 0.1041968162 0.46164978 0.4341534  1.3420757
## rs7639361   1396 0.9964311             1 0.6840974 0.31590258 0.0988538682 0.43409742 0.4670487  0.1626160
## rs675257    1391 0.9928622             1 0.7203451 0.27965492 0.0877066858 0.38389648 0.5283968 -1.7587956
## rs10836873  1398 0.9978587             1 0.9223891 0.07761087 0.0042918455 0.14663805 0.8490701  0.9044091

MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)

#-------------------------- Heterozygosity -------------------------------------

# Heterozygosity refers to the presence of each of two alleles at a given SNP 
# within an individual. This is expected under Hardy-Weinberg Equilibrium (HWE) 
# to occur with probability of 2p(1-p) where p is the dominant allele frequency 
# at that SNP (assuming a bi-allelic SNP). Excess or deficient heterozygosity 
# across typed SNPs within an individual is commonly considered as an indication 
# of inbreeding between relatives, but also serves as a screen for sample 
# contamination or poor sample quality.

# In this lab, heterozygosity F statistic calculation is carried out with the form, 
# |F|=(1-O/E), where O is observed proportion of heterozygous genotypes for a 
# given sample and E is the expected proportion of heterozygous genotypes for a 
# given sample based on the minor allele frequency across all non-missing SNPs for 
# a given sample. We remove any sample with inbreeding coefficient |F| > 0.10, 
# which might indicate autozygosity (which is not inherently a problem, but perhaps 
# unusual) or bad quality sample. We also remove any subjects with NA values.

hetExp <- callmatrix %*% (2*MAF*(1-MAF)) #Note that %*% is the 'multiply_by_matrix' function from plyr
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

plot(hetObs,hetExp);abline(0,1)
hist(snpsum.row$hetF)
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

#------------------------------- Call rate -------------------------------------

# To filter by call rate, we remove individuals that are missing genotype data 
# across more than a pre-defined percentage of the typed SNPs. This proportion of 
# missingness across SNPs is the sample call rate and we apply a threshold of 5%, 
# meaning that we drop those individuals missing genotype data for more than 5% 
# of the typed SNPs. A new reduced dimesion SnpMatrix gentotype is created from this filter.

sampcall <- 0.90   

sampleuse1 <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall)
sampleuse1[is.na(sampleuse1)] <- FALSE 
cat(nrow(genotype)-sum(sampleuse1), "individuals will be removed due to low sample call rate.\n")

genotype <- genotype[sampleuse1,]
clinical<- clinical[ rownames(genotype), ]
snpsum.row<-row.summary(genotype)

#------------------------ Hardy-Weinberg Equilibrium --------------------------
#how many variants are omitted comparing CAD controls only vs. the whole sample?

hardy <- 10^-6
controls <- as.character(clinical[ clinical$CAD==0, 'FamID' ]) 

snpsum.colControls <- col.summary(genotype[controls,])

HWEuse <- with(snpsum.colControls, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )

rm(snpsum.colControls)

HWEuse[is.na(HWEuse)] <- FALSE          

cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n") 

genotype <- genotype[,HWEuse]

print(genotype)

#----------- Cryptic relatedness, duplicates and gender identity ----------------

# Population-based cohort studies often only concerns unrelated individuals and 
# the generalized linear modeling approach assumes independence across individuals. 
# In regional cohort studies (e.g. hospital-based cohort studies) of complex 
# diseases, however, individuals from the same family could be unintentionally 
# recruited. To mitigate this problem, we thus employ Identity-by-descent (IBD) 
# analysis which is a common measurement of relatedness and duplication. Note that 
# gender identity can also be checked at this stage to confirm self-reported 
# gender consistency on X and Y chromosomes. With this data, we do not check gender 
# identity due to unavailability of sex chromosomes information.

# To filter on relatedness criteria, we use the SNPRelate package to perform 
# identity-by-descent (IBD) analysis. This package requires that the data be 
# transformed into a GDS format file. IBD analysis is performed on only a subset 
# of SNPs that are in linkage equilibrium by iteratively removing adjacent SNPs 
# that exceed an linkage disequilibrium (LD) threshold in a sliding window using 
# the snpgdsLDpruning() function.In this lab, we set LD threshold to be 20% and 
# kinship threshold to be 10%, which means that we would remove any SNPs whose 
# relevant statistics are higher than our stated threshold values. Then we used 
# snpgdsBED2GDS to create a GDS file of the data set from .bim, .bed and .fam files.

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

# print(genotype)


plot(kins$k0,kins$k1,xlim=c(0,1),ylim=c(0,1),xlab="IBD0",ylab="IBD1")
#plot(ibdcoeff$k0,ibdcoeff$k1,xlim=c(0,1),ylim=c(0,1),xlab="IBD0",ylab="IBD1") #this is ALL PAIRS of individuals, mostly minimally related (unrelated)

#-------------------------- PCA ------------------------------------------------

# Principle components (PCs) from multi-dimensional scaling (MDS) is one approach 
# to visualizing and classifying #individuals into ancestry groups based on their 
# observed genetics make-up. There are two reasons for doing this. First, 
# self-reported race/ethnicity can differ from clusters of individuals that are 
# based on genetics infomation. Second, the #presence of an individual not 
# appearing to fall within a racial/ethnic cluster may be suggestive of a sample error.

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=2)

pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10], stringsAsFactors = FALSE)

colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))

pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    
                    PC2 = pca$eigenvect[,2],    
                    stringsAsFactors = FALSE)

plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", main = "Ancestry Plot")

hist(pca$eigenvect[,1])
plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", main = "Ancestry Plot", pch = ".")

plot(pca$eigenval,ylab="Eigenvalue",xlab="Index of Eigenvector",xlim=c(1,20),pch=20,main="Scree plot")

# since we did not calculate all PCs, just the first few, we can't calculate the 
# proportion of variance captured by each PC, which is the eigenvalue / trace, 
# where the trace is the sum of all eigenvalues. Instead, we simply plot the eigenvalue.

closefn.gds(genofile)

save(genotype, snpsum.col, snpsum.row, genofile, clinical, pcs, genoMap, file= "lab2_MyGWAS_save.RData")
