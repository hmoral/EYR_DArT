############################################
## Calculate allelic frequencies of target #
## and random loci, plots distributions    #
## and significance of test vs random      #
###########################################

# writen by Hernan E. Morales
# please cite Morales et al. in preparation
# last modified 2016-07-22

# USAGE
# This script has to be used "manually"
# It does NOT take any input from the command line

##################
# Load libraries #
##################
library(ggplot2)
##################
# set directory  #
##################
setwd("./")

##################################
# Read input parameters and data #
##################################
# Read in allelic freq table
FREQ.TABLE.2<- read.table("./AlleleFreq.FILE.txt",header=T)

# Read inforation of Fst outliers loci (i.e. file with a q-value, Fst or simmilar to define outlier and background loci)
Marker.List<- read.table("FILE.BayeSc.txt",header=T)
names(Marker.List)[1]<- "Locus" # Make sure heade is the same as in FREQ.TABLE.2


### For Overlap outliers
Overlap.Outliers<- data.frame(Set = NULL, EYR_A = NULL, EYR_B = NULL) # make a empty table to store allelic freq correlations
SUBSET<- subset(Marker.List, Marker.List$OUTLIER.South == "YES" & Marker.List$OUTLIER.North == "YES") # make a list of outlier markers common across the entire species range
SELECTION<- as.character(SUBSET$Locus)
OUTLIERS <- FREQ.TABLE.2[FREQ.TABLE.2$Locus %in% SELECTION,] # Subset allelic freq table
EYR_A.cor<- cor(OUTLIERS$pEYR_A_Nth,OUTLIERS$pEYR_A_Sth) # Meassure correlation in the SOUTH (or whatever comparison) for all outlier markers
EYR_B.cor<- cor(OUTLIERS$pEYR_B_Nth,OUTLIERS$pEYR_B_Sth) # Meassure correlation in the NORTH (or whatever comparison) for all outlier markers
OUT.table<- data.frame(Set = "Outlier", EYR_A = EYR_A.cor, EYR_B = EYR_B.cor) #write table stating that these are Outlier loci
Overlap.Outliers<- rbind(Overlap.Outliers,OUT.table) # combine table
RANDOM.table<- FREQ.TABLE.2[!FREQ.TABLE.2$Locus %in% SELECTION,] # remove these markers from freq table
# loop through all other markers (background) and subset as many random sets as possible of the same size as "OUTLIERS"
# repeat the same correlation as above until there are not enough markers left to make a new random set
repeat{
  SEL<- sample(RANDOM.table$Locus, length(SELECTION))
  RANDOM.SEL<- RANDOM.table[RANDOM.table$Locus %in% SEL,]
  EYR_A.cor<- cor(RANDOM.SEL$pEYR_A_Nth,RANDOM.SEL$pEYR_A_Sth)
  EYR_B.cor<- cor(RANDOM.SEL$pEYR_B_Nth,RANDOM.SEL$pEYR_B_Sth)
  RND.table<- data.frame(Set = "Random", EYR_A = EYR_A.cor, EYR_B = EYR_B.cor)
  Overlap.Outliers<- rbind(Overlap.Outliers,RND.table)
  RANDOM.table<- RANDOM.table[!RANDOM.table$Locus %in% SEL,]
  if((nrow(RANDOM.table) > length(SELECTION)) != TRUE){
    break
  }
}
# Overlap.Outliers is the frequency of outlier loci and random loci
write.table(Overlap.Outliers, ".FILE.txt", sep="\t", quote=F, row.names = F)

# make a t.test, check that distribution is normal, etc and plot the distribution for each comparison
plot(density(Overlap.Outliers$EYR_A));
shapiro.test(Overlap.Outliers$EYR_A)
qqnorm(Overlap.Outliers$EYR_A);qqline(Overlap.Outliers$EYR_A, col = 2)
T.test<- t.test(Overlap.Outliers$EYR_A[2:nrow(Overlap.Outliers)], mu = Overlap.Outliers$EYR_A[1],
                alternative = c("less"))
ggplot(Overlap.Outliers, aes(x=EYR_A)) + geom_histogram(binwidth = 0.03) +
  annotate("text", x=Overlap.Outliers$EYR_A[1], y=4, label= "Outliers") +
  ggtitle(paste("Overlap.Outliers vs. background (N =", (nrow(Overlap.Outliers)-1), ")",
                "\n", "pval = ", round(T.test$p.value, digits = 6)))

T.test<- t.test(Overlap.Outliers$EYR_B[2:nrow(Overlap.Outliers)], mu = Overlap.Outliers$EYR_B[1],
                alternative = c("less"))
ggplot(Overlap.Outliers, aes(x=EYR_B)) + geom_histogram(binwidth = 0.03) +
  annotate("text", x=Overlap.Outliers$EYR_B[1], y=4, label= "Outliers") +
  ggtitle(paste("Overlap.Outliers vs. background (N =", (nrow(Overlap.Outliers)-1), ")",
                "\n", "pval = ", round(T.test$p.value, digits = 6)))
