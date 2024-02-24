# Here, we use a function we’ve written GWAA() to conduct the analysis. To prepare
# for it, we need to load the data saved from last module and install necessary 
# packages like GenABEL and doParallel. The GWAA function requires two arguments. 
# The genodata argument should specify the entire genotype data object in SnpMatrix 
# format. The phenodata argument should be a data frame with a column of sample 
# IDs, corresponding to the row names of genodata, followed by columns for the 
# outcome variable and covariates to be included in the model. In the following 
# code, we assign the genodata variable and write out a table called GWAA.txt for
# the results to later be written to.

rm(list=ls()) #remove any data objects or functions to avoid collisions or misspecification

load("lab2_MyGWAS_save.RData")

#BiocManager::install("doParallel",force=TRUE) #run only once to install

library(snpStats)
library(dplyr)
library(doParallel)

genodata <- genotype

# Print the number of SNPs to be checked
cat(paste(ncol(genodata), "Thus far, SNPs are included in analysis.\n"))


# create text file for GWAA output to be written to
columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")

# the output file handle is written in the mac style, with '~' referring to the current users home directory. 
write.table(t(columns), "~/Documents/ENES/5to/StatisticsModels/final_GWAS/GWAA.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

#---------------- Manipulating Phenotype Data ----------------------------------
# Next we focus on the phenodata argument. First we assign phenoSub to be 
# phenodata, which is created from merging the clinical and pcs data. We then 
# create the phenotype value (outcome variable) by using a normal 
# rank-transformation of the hdl variable, which we show is necessary later on. 
# Next we remove anyone that doesn’t have HDL data available, and subset the 
# genotype data for the remaining individuals. Lastly, we want to remove variables 
# that won’t be used in analysis, e.g. the remaining outcome variables and the 
# untransformed version of HDL (or whatever trait is specified by traitName).

traitName<-'IBD_1Control_2Case'

# Note that the linear model as coded takes all variables in the phenodata data.frame 
# as covariates, i.e. age sex and any pcs. One may control which covariates are 
# included in the phenodata data.frame, or by changing the glm in the GWAA 
# function to explicitly account for covariates

clinicalSub<-clinical[,c('FamID','sex','age',traitName)]
clinicalSub <- clinicalSub %>% rename(prePhenotype = traitName)
clinicalSub$phenotype <- ifelse(clinicalSub$prePhenotype == 1, 0, ifelse(clinicalSub$prePhenotype == 2, 1, clinicalSub$prePhenotype))  # Originally controls = 1 & cases = 2, now controls = 0 & cases = 1
#clinicalSub$phenotype <- rntransform(as.numeric(clinicalSub$prePhenotype, family="gaussian")) #TRANSFORMATION!!!
clinicalSub$prePhenotype <- NULL
pcsSub<-pcs[,c('FamID','pc1','pc2')] #add additional PCs here and merge in the next step

phenodata <- merge(clinicalSub,pcsSub,by="FamID")     
phenodata<-phenodata[!is.na(phenodata$phenotype),]
genodata <- genodata[as.character(phenodata$FamID),]

cat(nrow(genodata), "samples included in analysis thus far.\n")
#print(head(phenodata))

par(mfrow=c(1,2))
hist(as.numeric(clinical[,c(traitName)]), main=NULL, xlab=paste(traitName))
hist(phenodata$phenotype, main=NULL, xlab=paste("Transformed ",traitName))

#------------------------ The GWAA function ------------------------------------
# The function itself has four major elements. Given the data described above, 
# it (1) determines how many cores to use/which method to employ for parallel 
# processing (depending on specified cores and operating system), (2) converts 
# the genotype data into the necessary format for an additive model, (3) fits a 
# linear model for each SNP and covariates (one group at a time) and then (4) 
# writes the SNP coefficient details (the estimate, standard error, test statistic 
# and p value) from each model out to the GWAA.txt file. Read through the function 
# to get a general idea of the format in regards to these four steps. Further 
# details regarding the function elements are provided below.

# note that the number of hosts in the cluster defaults to half of those available 
# on the current machine, the default mc.cores=2

# note that each core uses about 250Mb memory to analyze all autosomes ~350k variants, 
# or 2Gb total.

# Low memory computers may want to reduce the number of hosts to 2 or even 1 to 
# limit use of RAM. (and/or close your memory hungry web browser while you run this analysis!)

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

# perhaps subset the data to make it run faster. Change the commented out lines below to control full GWAS vs. single chromosome
SNPs_sub <- genoMap$SNP[genoMap$chr==15] #this took about 30 seconds on 4 2.3 GHz i9 cores to run the GWAA function, Nov 14 2022
# SNPs_sub <- genoMap$SNP[genoMap$chr<=22] #this took considerably longer, about 20 minutes #updated Nov 14 2022

genodata_sub <- genodata[,colnames(genodata)%in%SNPs_sub]

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

# note that base R 'merge' threw an error that the vectors in the merge were too 
# long, which they shouldn't be, but this works as expected
gwasAnnotation<- gwas %>% left_join(genoMap, by = "SNP")

save(phenodata, genoMap, genodata_sub, GWAA, nSNPs, gwasAnnotation, file = "lab3MyGWAS_save.RData")
# You have just run a Genome Wide Association Study!