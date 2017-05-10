###########################################
## Define genomic regions of a given size #
## and counts annotated genes given a     #
## GO term, the second part of the script #
#  produce a random distribution of counts#
###########################################

# writen by Hernan E. Morales
# please cite Morales, H. E., et al. (2016). "Mitochondrial-nuclear interactions maintain a deep mitochondrial split in the face of nuclear gene flow." bioRxiv 095596. doi: https://doi.org/10.1101/09559.
# last modified 2017-01-22

# USAGE
# This script has to be used "manually"
# It does NOT take any input from the command line

##################
# Load libraries #
##################
 library(data.table)
library(ggplot2)
library(car)
##################
# set directory  #
##################
setwd("./")

##################################
# Read input parameters and data #
##################################

# Load table with all the annotated genes (e.g. GTF or GFF3 files from ENSEMBL)
# in this case I used the Zebra Finch genome
ZF_ALL<- fread("FILE", sep="\t", header=T)
ZF_ALL<- data.frame(ZF_ALL)
ZF_ALL<- ZF_ALL[!duplicated(ZF_ALL$Ensembl.Gene.ID),] # Eliminate dupicated entries
CHR_Table<- table(ZF_ALL$Chromosome.Name) # Confirm all chromsomes were included
# Edit chromosome names, delete unkown genomic region and the mitochondrial genome
ZF_ALL<- ZF_ALL[- grep("random", ZF_ALL$Chromosome.Name),] 
ZF_ALL<- ZF_ALL[- grep("Un", ZF_ALL$Chromosome.Name),]
ZF_ALL<- ZF_ALL[- grep("MT", ZF_ALL$Chromosome.Name),]
CHROM<- table(ZF_ALL$Chromosome.Name) # Confirm only desired chromsomes were included
print(CHROM)

# Load files of annotated genes (Ensembl.Gene.ID) with specific function
# For example, in this case I obtained:
# OXPHOS protein-coding of Zebra Finch in KEGG (GO term: 0006119, oxidative phosphorylation)
# All mitonuclear annotated genes in Zebra Finch from ENSEMBL (GO term: 0005739 , mitochondrion)
OXPHOS<- read.table("FILE", sep="\t", header=T)
OXPHOS<- unique(OXPHOS) # delete duplicated entries
OXPHOS<- data.frame(OXPHOS$Ensembl.Gene.ID) # Keep Ensembl Gene IDs only
mito_GO<- read.table("FILE", sep="\t", header=T)
mito_GO<- unique(mito_GO$Ensembl.Gene.ID) # delete duplicated entries
mito_GO<- data.frame(mito_GO$Ensembl.Gene.ID)  # Keep Ensembl Gene IDs only
# make sure both lists have the same column name

############# Count genes in areas of interest (e.g. divergent peaks) ##################

# Filter all genes by chromsomome (e.g. 1A, 4 or Z) and genomic coordinates
# When selecting genomic coordinates consider the rate of LD decay to capture linked genes (7.8 Kb in EYR, see Morales et al. 2016)

Chr.Table<- subset(ZF_ALL, ZF_ALL$Chromosome.Name == XX)
Chr.Table<- subset(Chr.Table, Chr.Table$GeneStart >= (XXXX - 7.8 Kb) & Chr.Table$GeneEnds <= (XXXX + 7.8 Kb))

# Count genes
OXPHOS_count<- nrow(merge(Chr.Table, OXPHOS, by = "Ensembl.Gene.ID"))
mito_GO_count<- nrow(merge(Chr.Table, mito_GO, by = "Ensembl.Gene.ID"))

############# Random count  ##################

# Define size of the divergent cluster
Start= min(Chr.Table$Gene.Start..bp.)
End= max(Chr.Table$Gene.End..bp.)
Peak_size<- (End - Start)

PEAK<- data.frame(PEAK$Ensembl.Gene.ID) # Rename divergent cluster frame
names(PEAK)<- "Ensembl.Gene.ID"
PEAK<- unique(PEAK)
GENES_PEAK<- nrow(PEAK) # Count total number of gene in cluster

CHR_Size<- data.frame(CHR=NULL,Size=NULL) # create empty data frame

# Loop goes through each chromosome, make sure that the chromosome is big enough to hold a cluster of gene of at least "Peak_size" lenght
# Define 1000 random genomic regions of "Peak_size" lenght and count OXPHOS and mitonuclear genes

for (i in names(CHROM)){
  x<- subset(ZF_ALL, ZF_ALL$Chromosome.Name == i)
  Range<- ((max(x$Gene.End..bp.)-min(x$Gene.Start..bp.)))/1000000 
  Table<- data.frame(CHR=i,Size=Range)
  CHR_Size<- rbind(CHR_Size,Table)
}
CHR_Large<- subset(CHR_Size, CHR_Size$Size >= Peak_size)

Match_Table<- data.frame(CHR="PEAK", START= Start, END= End, GENES= GENES_PEAK, OXPHOS=PEAK_OXPH_Count, mito=PEAK_Mito_Count)
for (i in CHR_Large$CHR){
  x<- subset(ZF_ALL, ZF_ALL$Chromosome.Name == i)
  Top<- max(x$Gene.End..bp.) - Peak_size
  Bottom<- min(x$Gene.Start..bp.)
  for (j in c(1:1000)){
    Start<- sample(Bottom:Top,1,replace=T)
    End<- Start + Peak_size
    y<- subset(x, x$Gene.Start..bp. >= Start)
    y<- subset(y, y$Gene.End..bp. <= End)
    StarT<- min(y$Gene.Start..bp.)
    EnD<- min(y$Gene.End..bp.)
    Range<- EnD-StarT
    Genes_PEAK<- nrow(y)
    OXPH_Count<- nrow(merge(y,OXPHOS, by="Ensembl.Gene.ID"))
    Mito_Count<- nrow(merge(y,mito_GO, by="Ensembl.Gene.ID"))
    z<- data.frame(CHR=i, START= StarT, END= EnD, GENES= Genes_PEAK, OXPHOS=OXPH_Count, mito=Mito_Count)
    Match_Table<- rbind(Match_Table,z)
  }
}

z<- z[!duplicated(Match_Table[c("START", "END", "GENES")]),] # Eliminate duplicated entries
# eliminate counts that fall within the "Peak_size" lenght by HAND (depends on were the are of interest is)
# "z" is the final table, save and plot distributions. To test if the count of your are of interest is greater than the random distribution is a t.test

# E.g. Plot random distrubutions
ggplot(z, aes(x=mito)) + geom_density(color = "black", fill = "black", alpha = 0.4)  +
  xlab("Genes with mitochondrial GO term") +
  ggtitle("Density distribution of mitochondrial GO terms\n accross random genomic regions")

ggplot(z, aes(x=OXPHOS)) + geom_histogram(color = "black", fill = "black", alpha = 0.4, binwidth=1, right = T) +
  xlab("Genes with OXPHOS GO term") + xlim(-1,5) +
  ggtitle("Histogram of OXPHOS GO terms\n accross  random genomic regions")
#
