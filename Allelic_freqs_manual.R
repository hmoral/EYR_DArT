#! /usr/bin/Rscript

#########################################
## Calculate allelic frequencies       ##
## from a SNP data (0,1,2) were        ##
## 0: Reference homozygous             ##
## 1: Alternative homozygous           ##
## 2: Heterozygous                     ##
#########################################

# writen by Hernan E. Morales
# please cite Morales, H. E., et al. (2016). "Mitochondrial-nuclear interactions maintain a deep mitochondrial split in the face of nuclear gene flow." bioRxiv 095596. doi: https://doi.org/10.1101/09559.
# last modified 2016-07-22

# USAGE
# This script has to be used "manually"
# It does NOT take any input from the command line

##################
# Load libraries #
##################
library(data.table)

##################
# set directory  #
##################
setwd("./")

##################################
# Read input parameters and data #
##################################
# Table in 012 format with headers with marker columns and sample names
# The table has n columns with individual names and n rows with marker names
SNPtable<- fread("SNP.TABLE.txt", header=T, sep=",", data.table=F)
# fread adds a prefix "X" to header names that start with number, if this is the case read the first row of the file again (header) and assign back to column names
# a work-around this problem is using "read.table" with check.names=FALSE, but would be very slow for large tables
x.1<- fread("SNP.TABLE.txt", header=F, sep=",",nrows = 1)
names(SNPtable)<- x.1

# Population file: text file with each individual assigned to one population. Individual names should be the same as in "SNPtable"
POPtable<- read.table("POP.FILE.txt")

# Make vector containing genotypes only (i.e. remove extra columns, i.e. Indv names, mapping information, etc) by subsetting data
geno<- SNPtable[,N:ncol(SNPtable)] # Change "N" with the number of column where genotypes start in SNPtable
# Add column with the name of the loci from SNPtable
row.names(geno)<- SNPtable$SNP_Code

# Dummy set with 1000 loci for code testin, uncomment to use
#geno<- geno[1:1000,]

## assign all non {0,1,2} to NA
geno[(geno=="NA")] <- NA # Or whatever character for missing data was used e.g. -9
# Transpose the geno object
geno<- data.frame(t(geno))
# Assign a new column  with the name of the individuals
geno$Code<- row.names(geno)
# merge the sample info table ABOVE with the geno object to assign populations to individuals
geno<- merge(geno,POPtable[c("Code", "POP")], by = "Code")
# remove first column of geno (Indv names)
geno<- geno[,-1]
# Make a vector with list of pops
POPS<- names(table(geno$POP))
# Make an empty frequence table
FREQ.TABLE.2<- data.frame(Locus= names(geno[,1:(ncol(geno)-1)]))
# loop through pops to make allelic frequency tables
for (i in POPS){
geno.pop<- subset(geno, geno$POP == i)
FREQ.TABLE.1<- data.frame(NULL)
# Loop thorugh markers to estimate frequency per marker
for (j in c(1:(ncol(geno.pop)-1))){
geno.pop.1<- data.frame(geno.pop[,j])
n0 <- apply(geno.pop.1==0,1,sum,na.rm=T)
n1 <- apply(geno.pop.1==1,1,sum,na.rm=T)
n2 <- apply(geno.pop.1==2,1,sum,na.rm=T)
n <- n0 + n1 + n2
p <- mean(((2*n0)+n2)/(2*n),na.rm = T)
q <- mean(1 - p, na.rm = T)
F.table<- data.frame(Locus = j, p = p, q = q)
names(F.table)<- c("Locus",paste("p",i,sep=""),paste("q",i,sep=""))
FREQ.TABLE.1<- rbind(FREQ.TABLE.1, F.table)
}
FREQ.TABLE.2<- cbind(FREQ.TABLE.2, FREQ.TABLE.1)
}
# FREQ.TABLE.2 is the final table with all allelic frequencies for all markers per population
head(FREQ.TABLE.2)

# check that correlations between pops are as expected
plot(FREQ.TABLE.2$POP1,FREQ.TABLE.2$POP2)

# Save table with frequencies
write.table(FREQ.TABLE.2, "FILE.txt",sep="\t",quote=F, row.names = F)
