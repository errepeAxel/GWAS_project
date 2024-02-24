#------------------- Read in the data ------------------------------------------
# In order to break down the structure of the different data files, we will 
# first read in the .bim, .bed and .fam files together using the command read.plink() 
# found in the snpStats package provided by Bioconductor. We will also be utilizing 
# the R package downloader in order to read in the data directly from the web and 
# write local versions to use. Please install the downloader package using the 
# command library(downloader) before preceeding.

rm(list=ls()) #remove any data objects or functions to avoid collisions or misspecification

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.17") 
#install.packages("RCurl")
#install.packages("downloader")
#install.packages("snpStats")
#BiocManager::install("snpStats")

library(RCurl)
library(downloader)
library(snpStats)

# read in data in two steps. 
## ------------------------------ step 1 --------------------------------------

bedFile<-"/home/axel/Documents/ENES/5to/StatisticsModels/lab/GWAS_data.bed"
bimFile<-"/home/axel/Documents/ENES/5to/StatisticsModels/lab/GWAS_data.bim"
famFile<-"/home/axel/Documents/ENES/5to/StatisticsModels/lab/GWAS_data.fam"

data <- read.plink(bedFile,bimFile,famFile)

clinicalFile <- "/home/axel/Documents/ENES/5to/StatisticsModels/final_GWAS/Clinical_IBD.csv"
clinical <- read.csv(clinicalFile, colClasses = c("character", rep("factor",2), rep("numeric", 2)))
##--------------------------- step 2 --------------------------------------------
# read in the local files using read.plink() in snpStats

# what class is "data"
class(data)
## [1] "list"

# how many elements are in "data"
length(data)
## [1] 3

# You will first notice that when read in this way, the data is of class “list” 
# and has 3 elements in it. Each element corresponds to the three kinds of data 
# that were read into R. We will next explore each component as each serves a 
# different purpose in analysis and will be needed in the upcoming labs.

##------------------------------- Exercise ------------------------------------
# Before we look at each piece of data in detail, try looking at the general 
# structure of it using the function str(). You’ll notice that there is information 
# on all three of the objects in the list. What does each element of the list 
# correspond to (e.g. what is it called and what appears to be in it)? Of what 
# ``class" is the genotypes object?

str(data$fam)
#'data.frame':	1401 obs. of  6 variables:
#  $ pedigree: int  10002 10004 10005 10007 10008 10009 10010 10011 10012 10013 ...
#$ member  : int  1 1 1 1 1 1 1 1 1 1 ...
#$ father  : logi  NA NA NA NA NA NA ...
#$ mother  : logi  NA NA NA NA NA NA ...
#$ sex     : int  1 2 1 1 1 1 1 2 1 2 ...
#$ affected: int  1 1 2 1 2 2 2 1 2 NA ...

str(data$map)
#'data.frame':	500000 obs. of  6 variables:
#$ chromosome: int  4 6 9 18 3 17 19 11 10 1 ...
#$ snp.name  : chr  "rs4579145" "rs2768995" "rs10125738" "rs888263" ...
#$ cM        : logi  NA NA NA NA NA NA ...
#$ position  : int  78303781 6911315 101066136 12750253 80138205 69675940 54268380 37859556 100361583 201052199 ...
#$ allele.1  : chr  "T" "C" "A" "G" ...
#$ allele.2  : chr  "C" "T" "G" "A" ...

str(data$genotypes)
#Formal class 'SnpMatrix' [package "snpStats"] with 1 slot
#..@ .Data: raw [1:1401, 1:500000] 03 02 03 03 ...
#.. ..- attr(*, "dimnames")=List of 2
#.. .. ..$ : chr [1:1401] "10002" "10004" "10005" "10007" ...
#.. .. ..$ : chr [1:500000] "rs4579145" "rs2768995" "rs10125738" "rs888263" ...
#..$ dim     : int [1:2] 1401 500000
#..$ dimnames:List of 2
#.. ..$ : chr [1:1401] "10002" "10004" "10005" "10007" ...
#.. ..$ : chr [1:500000] "rs4579145" "rs2768995" "rs10125738" "rs888263" ...

#------------------------------ The .bim file ----------------------------------
# The .bim file contains information on each SNP and can be seen by calling the 
# map element in the data object (note that .bim and .map file formats are similar).

bim <- data$map

head(bim)

##            chromosome   snp.name cM  position allele.1 allele.2
## rs4579145           4  rs4579145  0  78303781        T        C
## rs2768995           6  rs2768995  0   6911315        C        T
## rs10125738          9 rs10125738  0 101066136        A        G
## rs888263           18   rs888263  0  12750253        G        A
## rs7639361           3  rs7639361  0  80138205        A        C
## rs2430512          17  rs2430512  0  69675940        T        C

# As you can see from the first six observations displayed above, this file has a 
# row for each SNP and six columns describing the SNPs. The columns correspond to 
# chromosome number, RS number (SNP identifier), genetic distance, the position on 
# the chromosome and the two alleles, respectively. Note that in most contexts, 
# “allele 1” is also referred to as the “minor allele” and “allele 2” refers to 
# the “major allele” at each SNP. It is also important to note that SNP names and 
# locations can change based on the genetic build your data is based on. Here, 
# the RS numbers and positions are from build Hg19.


##------------------------------- Exercise -------------------------------------
# Start investigating the data at hand. How many SNPs are present? Execute the 
# command table(bim$chromosome). Which chromosomes are represented? Which has the 
# most SNPs in this data?

table(bim$chromosome)

# How many SNPs are present?
sum(table(bim$chromosome))
## There are 500,000.

# Which chromosomes are represented? 
## Everyone without the X and Y.

# Which has the most SNP in this data? 
which.max(table(bim$chromosome))
## El cromosoma 2

#------------------------------ The .bed file ----------------------------------
# The raw version of the .bed file contains a binary version of the genotype data. 
# This is the largest of the three files because it contains every SNP in the study, 
# as well as the genotype at that snp for each individual. When we read this file 
# into R, however, it’s automatically translated into the genotype information from 
# the binary code. R then store this object as a SnpMatrix object that is critical 
# to SNP analysis. This SnpMatrix corresponds to the genotypes element in the data object.

bed <- data$genotypes
head(data$genotypes)
bed

## A SnpMatrix with  1401 rows and  500000 columns
## Row names:  10002 ... 11596 <-- Individuals
## Col names:  rs4579145 ... rs946221 <-- SNPs


##-------------------------------- Exercise ------------------------------------
# What do the columns in this matrix represent? What do the rows represent? How many 
# SNPs are in this matrix? How many individuals are present? What does the element 
# corresponding to the ith row and jth column represent?

# ANSWER The columns of this matrix represents info about each SNP.
# The rows represents individuals
# There are 500,000 SNPs and 1401 individuals.
# The element corresponding to the ith row and jth column represents the genotype 
# of the ith-individual at their jth-snp.

#------------------------------ The .fam file ----------------------------------
# Finally, the .fam file contains the participant identification information, 
# including a row for each individual and six columns, corresponding to “Family 
# ID Number”, “Sample ID Number”,“Paternal ID Number”, “Maternal ID Number”, “Sex”, 
# and “Phenotype”. The .fam file can be called by referring to the .fam element 
# of the data object.

fam <- data$fam

head(fam)

##       FamID IndID PatID MatID sex phenotype
## 10002    10002      1      0      0   1        1
## 10004    10004      1      0      0   2        1
## 10005    10005      1      0      0   1        2
## 10007    10007      1      0      0   1        1
## 10008    10008      1      0      0   1        2
## 10009    10009      1      0      0   1        2

# Note that not all of these columns contain information depending the nature of
# the data collection process and study design. In this case, we also have a 
# separate clinical file with several outcomes and covariates that will be read 
# into R in the next section.

##---------------------------------- Exercise ----------------------------------
# How many observations are there in the fam file? How many males are there 
# (assuming males are coded as 1)? How many females? How many levels of the 
# phenotype variable are there? What does the "pedigree" column correspond to? 
# Does it look familar (where else have you seen it)?

# How many observations are there in the fam file?
length(fam[,1]) # ANSWER 1401
# How many males are there?
male <- fam$sex == 1
sum(male) # ANSWER 937
# How many levels of the phenotype variable are there?
as.factor(fam$affected) # ANSWER Just two: 1 and 2
#  What does the "pedigree" column correspond to?
# ANSWER To the family they belong

#----------------------------- The clinical data -------------------------------
# Finally, in some instances you may have access to a separate phenotype file that 
# contains information on things like additional covariates for each individual 
# and other “outcome variables” you may be interested in (e.g. something other 
# than affected or variable that codes affected differently). The rows of this 
# file represent each subject and the columns correspond to the available 
# covariates and phenotypes. We read this data in separately as it is a simple .txt file.

rownames(clinical) <- clinical$FamID
print(head(clinical))

## ID    FamID CAD    sex age tg  hdl  ldl IBD_1Control_2Case
## 10002 10002   1   1  60 <NA> <NA> <NA>                  1
## 10004 10004   1   2  50   55   23   75                  1
## 10005 10005   1   1  55  105   37   69                  1
## 10007 10007   1   1  52  314   54  108                  2
## 10008 10008   1   1  58  161   40   94                  2
## 10009 10009   1   1  59  171   46   92                  1

##----------------------------- Exercise ---------------------------------------
# Explore the clinical data. How many complete observations are there? What is 
# the age range of study participants? Which column of this data set links it to
# all the others? How many people are CAD+ (CAD=1)?

# How many complete observations are there?
sum(complete.cases(clinical)) # ANSWER 1281
# What is the age range of study participants?
summary(clinical$age) # ANSWER From 22 to 87
# Which column of this data set links it to all the others? 
## ANSWER The 'FamID' column
# How many people are CAD+ (CAD=1)?
cad_positive <- clinical$CAD == 1
sum(cad_positive) # ANSWER 933

#------------------------------------------------------------------------------
genotype <- data$genotypes
genoMap <- data$map
colnames(genoMap)<-c("chr", "SNP", "gen.dist", "position", "A1", "A2")

snpsum.col <- col.summary(genotype)
help(col.summary)
head(snpsum.col)

#-------------------------------- Call rate ------------------------------------
# For a given SNP, its call rate is defined as the proportion of inidividuals in 
# our study for which the corresponding SNP information is not missing. For 
# example, a call rate of 95% for a certain SNP means that 95% of the individuals 
# have data for this SNP. Thus to filter by call rate, we look at the column of 
# Call.rate in our dataframe and choose the SNPs that has percentage of missing 
# observations below our chosen threshold. In this module we set a call rate 
# threshold to be 0.95 and then select for SNPs that pass the call rate.

call <- 0.95
keep <- with(snpsum.col, (!is.na(Call.rate) & Call.rate >= call))
head(keep)

keep[is.na(keep)] <- FALSE

cat(ncol(genotype)-sum(keep),"SNPs will be removed due to low call rate.\n")

genotype <- genotype[,keep]
snpsum.col <- snpsum.col[keep,]

print(genotype)
summary(genotype)

#----------------------------- Minor allele frequency --------------------------
# Inadequate power to infer a statiscally significant relationship between the 
# SNP and the trait under study is the result of a large degree of homogeneity at
# a given SNP across study participants. This occurs when we have a very small 
# minor allele frequency (MAF), which means the majority of the individuals have
# two copies of the same major alleles. Minor allele frequency (also known as 
# variant allele frequency) is defined as frequency of the less common allele at 
# a variable site. In this module, we remove SNPs whose minor allele fequency is 
# less than 1%. In other cases, particularly when sample size is small, we may 
# apply a cut point of 5%.

minor <- 0.01

keep1 <- with(snpsum.col, (!is.na(MAF) & MAF > minor) )
keep1[is.na(keep1)] <- FALSE

cat(ncol(genotype)-sum(keep1),"SNPs will be removed due to low MAF .\n"  )

# After knowing the SNPs needed to be removed, we subset the genotype and SNP 
# summary data for the SNPs that pass MAF criteria.

genotype <- genotype[,keep1] 
snpsum.col <- snpsum.col[keep1,] 

#----------------------- Saving work for following labs… -----------------------
# Now that all the data is read into R, we need to save the objects for the rest 
# of the labs. Note that some of the objects must first be renamed for the purpose 
# of future labs. You can rename and save the necessary data by running the following commands:

save(genotype,genoMap,snpsum.col,genoMap,bedFile,bimFile,famFile, clinical, file= "lab1_MyGWAS_save.RData")
