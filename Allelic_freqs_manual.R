#! /usr/bin/Rscript

#########################################
## Calculate allelic frequencies       ##
## from a SNP data (0,1,2) were        ##
## 0: Reference homozygous             ##
## 1: Alternative homozygous           ##
## 2: Heterozygous                     ##
#########################################

# writen by Hernan E. Morales
# please cite Morales et al. in preparation
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
# Read in the 012 table, the whole table and first row for name of samples
# The table has n columns with individual names and n rows with marker names
x<- fread("SNP.TABLE.txt", header=T, sep=",")
x<- data.frame(x)
# fread complains if individual names start with number (it adds a "X."), if this is the case read the first row of the file again (header) and ssign it as names
x.1<- fread("../../SNP_Tables/SNP_Table_DArT_hetero_error.012.transectSAMPLES.csv", header=F, sep=",",nrows = 1)
names(x)<- x.1

# Read info file with samples populations, with same individual names as "x" (SNP.TABLE.txt)
z<- read.table("FILE.txt")

# Make vector containing genotypes only (i.e. remove extra columns, i.e. Indv names, mapping information, etc) by subsetting data
geno<- x[,12:ncol(x)] # Change "12" with the number of column where genotypes start in "x" (SNP.TABLE.txt)
# Add column with the name of the loci
row.names(geno)<- x$SNP_Code

# Dummy set with 1000 loci for test code
#geno<- geno[1:1000,]

## assign all non {0,1,2} to NA
geno[(geno=="NA")] <- NA
# Transpose the geno object
geno<- data.frame(t(geno))
# Assign a new column  with the name of the individuals
geno$Code<- row.names(geno)
# merge the sample info table ABOVE with the geno object to assign populations to individuals
geno<- merge(geno,z[c("Code", "POP")], by = "Code")
# remove first column of geno (names)
geno<- geno[,-1]
# Make a vector with list of pops
POPS<- names(table(geno$POP))
# Make an emtpty frequence table
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
