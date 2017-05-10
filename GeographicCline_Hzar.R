#! /usr/bin/Rscript
#########################################
## Models geographic clines for SNP    ##
## with HZar using 3 cline models:     ##
## ModelI: fixed parameters no tails   ##
## ModelII: free parameters no tails   ##
## ModelIII: free parameters both tails##
## Neutral model: no cline             ##
#########################################

# based on Derryberry et al. 2014 Mol. Ecol Res
# adapted by Hernan Morales (hern.moral@gmail.com)
# please cite:
# Derryberry, E. P., et al. (2014). "HZAR: hybrid zone analysis using an R software package." Molecular ecology resources 14(3): 652-663
# and if applies also: Morales, H. E., et al. (2016). "Mitochondrial-nuclear interactions maintain a deep mitochondrial split in the face of nuclear gene flow." bioRxiv 095596. doi: https://doi.org/10.1101/09559.
# last modified 2016-08-19

# USAGE
# REQUIRED
# R-libraries doMC, hzar, data.table

# OUTPUT, directories with these names should be created
# - Best Models - list object for the best model for each SNP
# - Result Model - list object for all models for each SNP
# - Result Plots - pdf plots with comparison of the 4 models and best models (if best model = neutral, no plot is produced), object for each of these plots is also saved for each SNP
# - Result tables:
        # AIC_table: for each SNP (rows) AIC score for each  model and name of best model
        # Model_convergence_table: three tables, one for each model. Convergence of each paramter, TRUE/FALSE
        # Model.RESULT.Table: three tables, one for each model. Parameter estimates and CI's
        # Check4ModelIncluded: sanity check that all models were included
# - TRACE plots: three plots, one for each model

##################
# Load libraries #
##################
library(hzar)
library(data.table)
library(doMC)

# Set directory
setwd("./")

# Check if output directory exist and if not create one
ifelse(!dir.exists("BestModels"), dir.create("BestModels"), FALSE)
ifelse(!dir.exists("ResultPlots"), dir.create("ResultPlots"), FALSE)
ifelse(!dir.exists("ResultTables"), dir.create("ResultTables"), FALSE)

# Reads input parameters
args<-commandArgs(trailingOnly=T)

# Checks if enough input parameters were given
if(length(args)<4){
  stop("Usage: Hzar_cline_fitting.R filename.suffix MarkerList transect chainLenght ncores \n\t\t- filename.suffix is a file with SNPs, format: 0= homo, 0.5=hetero, 1=homo. Mandatory columns also include Distance of each individual NOTE my file has many more columns and the code is set for SNP to start in the 9th column \n\t\t- MarkerList file with markers to be included \n\t\t- transect: either NORTH or SOUTH (case sensitive)\n\t\t- chainLenght for the MCMC (e.g. 1e5) \n\t\t- ncores is the number of cores allowed for parallel computation") 
}

# Create empty tables
modelI.Convergence.TABLE<- c("SNP","Center", "width", "model.LL")
modelII.Convergence.TABLE<- c("SNP","Center", "width","pMin","pMax", "model.LL")
modelIII.Convergence.TABLE<- c("SNP","Center", "width","pMin","pMax","deltaL","tauL","deltaR","tauR", "model.LL")
FourModels.Table<- c("SNP","nullModel","modelI","modelII","modelIII")
AIC.Table<- data.frame(nullModel=NULL,modelI=NULL,modelII=NULL,modelIII=NULL,Min.AIC=NULL,SNP=NULL)
modelI.Table<-data.frame(Model=NULL,SNP=NULL,center=NULL,width=NULL,center2LLLow=NULL,center2LLHigh=NULL,width2LLLow=NULL,width2LLHigh=NULL)
modelII.Table<-data.frame(Model=NULL,SNP=NULL,center=NULL,width=NULL,pMin=NULL,pMax=NULL,center2LLLow=NULL,center2LLHigh=NULL,width2LLLow=NULL,width2LLHigh=NULL,pMin2LLLow=NULL,pMin2LLHigh=NULL,pMax2LLLow=NULL,pMax2LLHig=NULL)
modelIII.Table<-data.frame(Model=NULL,SNP=NULL,center=NULL,width=NULL,pMin=NULL,pMax=NULL,deltaL=NULL,tauLdeltaR=NULL,tauR=NULL,center2LLLow=NULL,center2LLHigh=NULL,width2LLLow=NULL,width2LLHigh=NULL,pMin2LLLow=NULL,pMin2LLHigh=NULL,pMax2LLLow=NULL,pMax2LLHig=NULL,deltaL2LLLow=NULL,deltaL2LLHigh=NULL,tauL2LLLow=NULL,tauL2LLHigh=NULL,deltaR2LLLow=NULL,deltaR2LLHigh=NULL,tauR2LLLow=NULL,tauR2LLHigh=NULL)


## A typical chain length.  This value is the default setting in the package.
chainLength=as.numeric(args[4])   
# Specify which transect
TRANSECT<- paste(args[3],"/",sep="")

#Load data
MAIN_data<- fread(args[1],header=T,sep="\t")
MAIN_data<- data.frame(MAIN_data)

#Subset selected markers
# names of all markers
MarkerNames<- names(MAIN_data)
# Read marker list and subset markers that are included in the set
MARKERS<- read.table(args[2],header=T)
MARKERS<- as.character(MARKERS[,1])
MARKERS<- MarkerNames[MarkerNames %in% MARKERS]


# Default number of random starting parameters
if(is.na(args[5])){args[5]=2}

# Start loop thorugh markers
for (i in MARKERS){
## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528))


if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
registerDoMC(as.numeric(args[5]))
} else {
  ## Use foreach in sequential mode
registerDoSEQ();
}

## ## Picking an allele for a locus
KEPP<- names(MAIN_data)[1:8]
DATA<- data.frame(MAIN_data[c(KEPP,i)])
COL<- ncol(DATA)
names(DATA)[COL]<-"SNP"
DATA$Sample<- 1
## Blank out space in memory to hold molecular analysis
if(length(apropos("^mkn$",ignore.case=FALSE)) == 0 ||
   !is.list(mkn) ) mkn <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
mkn$AdaA <- list();
## Space to hold the observed data
mkn$AdaA$obs <- list();
## Space to hold the models to fit
mkn$AdaA$models <- list();
## Space to hold the compiled fit requests
mkn$AdaA$fitRs <- list();
## Space to hold the output data chains
mkn$AdaA$runs <- list();
## Space to hold the analysed data
mkn$AdaA$analysis <- list();



## Locus Ada, Allele A from Brumfield et al 2001
mkn$AdaA$obs <-
  hzar.doMolecularData1DPops(DATA$Distance,
                             as.numeric(as.character(DATA$SNP)),
                             DATA$Sample);
#hzar.plot.obsData(mkn$AdaA$obs)

## Make a helper function
mkn.loadAdaAmodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep="."))
  mkn$AdaA$models[[id]] <<- hzar.makeCline1DFreq(mkn$AdaA$obs, scaling, tails)

mkn.loadAdaAmodel("fixed","none","modelI");
mkn.loadAdaAmodel("free" ,"none","modelII");
mkn.loadAdaAmodel("free" ,"both","modelIII");

## Check the default settings
print(mkn$AdaA$models)

## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between 0 and 570 km.
MINdistance<- min(DATA$Distance)-30
MAXdistance<- max(DATA$Distance)+30
mkn$AdaA$models <- sapply(mkn$AdaA$models,
                          hzar.model.addBoxReq,
                          MINdistance , MAXdistance,
                          simplify=FALSE)

## Check the updated settings
print(mkn$AdaA$models)

## Compile each of the models to prepare for fitting
mkn$AdaA$fitRs$init <- sapply(mkn$AdaA$models,
                              hzar.first.fitRequest.old.ML,
                              obsData=mkn$AdaA$obs,
                              verbose=FALSE,
                              simplify=FALSE)

## Update the settings for the fitter if desired.
mkn$AdaA$fitRs$init$modelI$mcmcParam$chainLength <-
  chainLength;                          #1e5
mkn$AdaA$fitRs$init$modelI$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
mkn$AdaA$fitRs$init$modelI$mcmcParam$seed[[1]] <-
  mainSeed$A

mkn$AdaA$fitRs$init$modelII$mcmcParam$chainLength <-
  chainLength;                          #1e5
mkn$AdaA$fitRs$init$modelII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
mkn$AdaA$fitRs$init$modelII$mcmcParam$seed[[1]] <-
  mainSeed$B 

mkn$AdaA$fitRs$init$modelIII$mcmcParam$chainLength <-
  chainLength;                          #1e5
mkn$AdaA$fitRs$init$modelIII$mcmcParam$burnin <-
  chainLength %/% 10;                   #1e4
mkn$AdaA$fitRs$init$modelIII$mcmcParam$seed[[1]] <-
  mainSeed$C 

## Check fit request settings
#print(mkn$AdaA$fitRs$init)

## Run initial chains for all models
mkn$AdaA$runs$init <- list()
mkn$AdaA$runs$init$modelI <-
  hzar.doFit(mkn$AdaA$fitRs$init$modelI)
mkn$AdaA$runs$init$modelII <-
  hzar.doFit(mkn$AdaA$fitRs$init$modelII)
mkn$AdaA$runs$init$modelIII <-
  hzar.doFit(mkn$AdaA$fitRs$init$modelIII)



#Store intial chain resulting parameters
modelI.center<-mkn$AdaA$runs$init$modelI$modelParam$upper$center
modelII.center<-mkn$AdaA$runs$init$modelII$modelParam$upper$center
modelIII.center<-mkn$AdaA$runs$init$modelIII$modelParam$upper$center
modelI.width<-mkn$AdaA$runs$init$modelI$modelParam$upper$width
modelII.width<-mkn$AdaA$runs$init$modelII$modelParam$upper$width
modelIII.width<-mkn$AdaA$runs$init$modelIII$modelParam$upper$width
modelII.pMin<-mkn$AdaA$runs$init$modelII$modelParam$upper$pMin
modelIII.pMin<-mkn$AdaA$runs$init$modelIII$modelParam$upper$pMin
modelII.pMax<-mkn$AdaA$runs$init$modelII$modelParam$upper$pMax
modelIII.pMax<-mkn$AdaA$runs$init$modelIII$modelParam$upper$pMax
modelIII.deltaR<-mkn$AdaA$runs$init$modelIII$modelParam$upper$deltaR
modelIII.deltaL<-mkn$AdaA$runs$init$modelIII$modelParam$upper$deltaL
modelIII.tauR<-mkn$AdaA$runs$init$modelIII$modelParam$upper$tauR
modelIII.tauL<-mkn$AdaA$runs$init$modelIII$modelParam$upper$tauL

## Compile a new set of fit requests using the initial chains 
mkn$AdaA$fitRs$chains <-lapply(mkn$AdaA$runs$init,hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
mkn$AdaA$fitRs$chains <-
  hzar.multiFitRequest(mkn$AdaA$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit
## runif(9,-30,600) center for modelI, modelII, modelIII
for (j in c(1:9)){
mkn$AdaA$fitRs$chains[[j]]$modelParam$init["center"]= runif(1,(modelI.center-(modelI.center*0.6)),(modelI.center+(modelI.center*0.6)))
}
## runif(9,0,630) width for modelI, modelII, modelIII
for (j in c(1:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["width"]= runif(1,(modelI.width-(modelI.width*0.6)),(modelI.width+(modelI.width*0.6)))
}
## runif(6,0,1) pMin for modelII, modelIII
for (j in c(4:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["pMin"]= runif(1,(modelII.pMin-(modelII.pMin*0.6)),(modelII.pMin+(modelII.pMin*0.6)))
}
## runif(6,0,1) pMax for modelII, modelIII 
for (j in c(4:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["pMax"]= runif(1,(modelII.pMax-(modelII.pMax*0.6)),(modelII.pMax+(modelII.pMax*0.6)))
}
## runif(3,0,630) deltaR for modelIII
for (j in c(7:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["deltaR"]= runif(1,(modelIII.deltaR-(modelIII.deltaR*0.6)),(modelIII.deltaR+(modelIII.deltaR*0.6)))
}
## runif(3,0,630) deltaL for modelIII
for (j in c(7:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["deltaL"]= runif(1,(modelIII.deltaL-(modelIII.deltaL*0.6)),(modelIII.deltaL+(modelIII.deltaL*0.6)))
}
## runif(3,0,630) tauL for modelIII
for (j in c(7:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["tauL"]= runif(1,(modelIII.tauL-(modelIII.tauL*0.6)),(modelIII.tauL+(modelIII.tauL*0.6)))
}
## runif(3,0,1) tauR for modelIII
for (j in c(7:9)){
  mkn$AdaA$fitRs$chains[[j]]$modelParam$init["tauR"]= runif(1,(modelIII.tauR-(modelIII.tauR*0.6)),(modelIII.tauR+(modelIII.tauR*0.6)))
}

## Go ahead and run a chain of 3 runs for every fit request
mkn$AdaA$runs$chains <-  hzar.doChain.multi(mkn$AdaA$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did modelI converge?
modelI.Convergence<- summary(do.call(mcmc.list,
                lapply(mkn$AdaA$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
modelI.a<- (modelI.Convergence$statistics[1]*0.25>modelI.Convergence$statistics[4])
modelI.b<- (modelI.Convergence$statistics[2]*0.25>modelI.Convergence$statistics[5])
modelI.c<- (modelI.Convergence$statistics[3]*0.25>modelI.Convergence$statistics[6])
modelI.Convergence<- c(i,TRUE == c(modelI.a, modelI.b, modelI.c))
modelI.Convergence.TABLE<-rbind(modelI.Convergence.TABLE,modelI.Convergence) 
pdf(paste("./TracePlots/",TRANSECT,"ModelI_",i,".pdf",sep=""))
plot(do.call(mcmc.list,
             lapply(mkn$AdaA$runs$chains[1:3],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )
dev.off()

## Yes it did.

## Did modelII converge?
modelII.Convergence<-summary(do.call(mcmc.list,
                lapply(mkn$AdaA$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

modelII.a<- (modelII.Convergence$statistics[1]*0.25>modelII.Convergence$statistics[6])
modelII.b<- (modelII.Convergence$statistics[2]*0.25>modelII.Convergence$statistics[7])
modelII.c<- (modelII.Convergence$statistics[3]*0.25>modelII.Convergence$statistics[8])
modelII.d<- (modelII.Convergence$statistics[4]*0.25>modelII.Convergence$statistics[9])
modelII.e<- (modelII.Convergence$statistics[5]*0.25>modelII.Convergence$statistics[10])
modelII.Convergence<- c(i,TRUE == c(modelII.a, modelII.b, modelII.c,modelII.d,modelII.e))
modelII.Convergence.TABLE<-rbind(modelII.Convergence.TABLE,modelII.Convergence) 
pdf(paste("./TracePlots/",TRANSECT,"ModelII_",i,".pdf",sep=""))
plot(do.call(mcmc.list,
             lapply(mkn$AdaA$runs$chains[4:6],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )
dev.off()

## Did modelIII converge?
modelIII.Convergence<-summary(do.call(mcmc.list,
                                      lapply(mkn$AdaA$runs$chains[7:9],
                                             function(x) hzar.mcmc.bindLL(x[[3]]) )) )
modelIII.a<- (modelIII.Convergence$statistics[1]*0.25>modelIII.Convergence$statistics[10])
modelIII.b<- (modelIII.Convergence$statistics[2]*0.25>modelIII.Convergence$statistics[11])
modelIII.c<- (modelIII.Convergence$statistics[3]*0.25>modelIII.Convergence$statistics[12])
modelIII.d<- (modelIII.Convergence$statistics[4]*0.25>modelIII.Convergence$statistics[13])
modelIII.e<- (modelIII.Convergence$statistics[5]*0.25>modelIII.Convergence$statistics[14])
modelIII.f<- (modelIII.Convergence$statistics[6]*0.25>modelIII.Convergence$statistics[15])
modelIII.g<- (modelIII.Convergence$statistics[7]*0.25>modelIII.Convergence$statistics[16])
modelIII.h<- (modelIII.Convergence$statistics[8]*0.25>modelIII.Convergence$statistics[17])
modelIII.i<- (modelIII.Convergence$statistics[9]*0.25>modelIII.Convergence$statistics[18])
modelIII.Convergence<- c(i,TRUE == c(modelIII.a, modelIII.b, modelIII.c,modelIII.d,modelIII.e,modelIII.f, modelIII.g, modelIII.h,modelIII.i))
modelIII.Convergence.TABLE<-rbind(modelIII.Convergence.TABLE,modelIII.Convergence) 
pdf(paste("./TracePlots/",TRANSECT,"modelIII_",i,".pdf",sep=""))
plot(do.call(mcmc.list,
             lapply(mkn$AdaA$runs$chains[7:9],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )
dev.off()


## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele
## frequency independent of distance along cline) to include in
## analysis.
mkn$AdaA$analysis$initDGs <- list(
  nullModel =  hzar.dataGroup.null(mkn$AdaA$obs))

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
mkn$AdaA$analysis$initDGs$modelI <-
  hzar.dataGroup.add(mkn$AdaA$runs$init$modelI)
mkn$AdaA$analysis$initDGs$modelII <-
  hzar.dataGroup.add(mkn$AdaA$runs$init$modelII)
mkn$AdaA$analysis$initDGs$modelIII <-
  hzar.dataGroup.add(mkn$AdaA$runs$init$modelIII)

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (nullModel, modelI,
## modelII, modelIII).
mkn$AdaA$analysis$oDG <-
  hzar.make.obsDataGroup(mkn$AdaA$analysis$initDGs)
mkn$AdaA$analysis$oDG <-
  hzar.copyModelLabels(mkn$AdaA$analysis$initDGs,
                       mkn$AdaA$analysis$oDG)

## Convert all 27 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
mkn$AdaA$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn$AdaA$runs$chains,
                                hzar.dataGroup.add),
                         mkn$AdaA$analysis$oDG);


## Check to make sure that there are only four hzar.dataGroup
## objects named nullModel, modelI, modelII, and modelIII in the
## hzar.obsDataGroup object.
FourModels<-(summary(mkn$AdaA$analysis$oDG$data.groups))
FourModels<-c(i,rownames(FourModels))
FourModels.Table<- rbind(FourModels.Table,FourModels)

## Compare the 3 cline models to the null model graphically
hzar.plot.cline(mkn$AdaA$analysis$oDG,main=paste("3 models",i))
assign(paste("AllModelCline_", i,sep=""), recordPlot())
#print(AllModelCline_ak2.A)

#Save object with all models
assign(paste("ALL_models_results",i,sep=""),mkn$AdaA$analysis$oDG)
All_models<- get(paste("ALL_models_results",i,sep=""))
saveRDS(All_models, paste("./ResultModels/",TRANSECT,"ALL_models_results",i,sep=""), ascii=TRUE)

## Do model selection based on the AICc scores
mkn$AdaA$analysis$AICcTable<-hzar.AICc.hzar.obsDataGroup(mkn$AdaA$analysis$oDG)
# build report table of AIC scores
AIC.t<-mkn$AdaA$analysis$AICcTable
AIC.t<- data.frame(t(AIC.t))
Min.AIC <- t(sapply(seq(nrow(AIC.t)), function(i) {
  j <- which.min(AIC.t[i,])
  colnames(AIC.t)[j]
}))
AIC.t$Min.AIC<- Min.AIC
AIC.t$SNP<- i
AIC.Table<- rbind(AIC.Table,AIC.t)

## Print out the model with the minimum AICc score
mkn$AdaA$analysis$model.name <-
        rownames(mkn$AdaA$analysis$AICcTable
        )[[ which.min(mkn$AdaA$analysis$AICcTable$AICc )]]
Model_Name<- mkn$AdaA$analysis$model.name

## Extract the hzar.dataGroup object for the selected model
mkn$AdaA$analysis$model.selected <-
  mkn$AdaA$analysis$oDG$data.groups[[mkn$AdaA$analysis$model.name]]
x<-mkn$AdaA$analysis$model.selected
assign(paste("BestModel_",i,sep=""), x)
save(x,file = paste("./BestModels/",TRANSECT,i,"_",Model_Name,sep=""))


## Plot best model
if(Model_Name != "nullModel"){
hzar.plot.fzCline(mkn$AdaA$analysis$model.selected, main= paste("BestModel",i))
assign(paste("Cline_bestModel_", i,sep=""), recordPlot())}
#print(Cline_bestModel1)

#Make table of models summary and save each model
#Model 1
modelI.Estimates<- hzar.get.ML.cline(mkn$AdaA$analysis$oDG$data.groups$modelI)
modelI.Estimates<- data.frame(modelI.Estimates$param.free)
modelI.CI<- data.frame(hzar.getLLCutParam(mkn$AdaA$analysis$oDG$data.groups$modelI,
                                          names(mkn$AdaA$analysis$oDG$data.groups$modelI$data.param)))
modelI.T<- data.frame(Model= "modelI", SNP= i)
modelI.T<- cbind(modelI.T,modelI.Estimates,modelI.CI)
modelI.Table<- rbind(modelI.Table,modelI.T)
assign(paste("modeI_results_SNP_",i,sep=""),mkn$AdaA$analysis$oDG$data.groups$modelI)
Data_modelI<- get(paste("modeI_results_SNP_",i,sep=""))
saveRDS(Data_modelI, paste("./ResultModels/",TRANSECT,"modelI_results_SNP_",i,sep=""), ascii=TRUE)
#Model 2
modelII.Estimates<- hzar.get.ML.cline(mkn$AdaA$analysis$oDG$data.groups$modelII)
modelII.Estimates<- data.frame(modelII.Estimates$param.free)
modelII.CI<- data.frame(hzar.getLLCutParam(mkn$AdaA$analysis$oDG$data.groups$modelII,
                                           names(mkn$AdaA$analysis$oDG$data.groups$modelII$data.param)))
modelII.T<- data.frame(Model= "modelII", SNP= i)
modelII.T<- cbind(modelII.T,modelII.Estimates,modelII.CI)
modelII.Table<- rbind(modelII.Table,modelII.T)
assign(paste("modeI_results_SNP_",i,sep=""),mkn$AdaA$analysis$oDG$data.groups$modelII)
Data_modelII<- get(paste("modeI_results_SNP_",i,sep=""))
saveRDS(Data_modelII, paste("./ResultModels/",TRANSECT,"modelII_results_SNP_",i,sep=""), ascii=TRUE)
#Model 3
modelIII.Estimates<- hzar.get.ML.cline(mkn$AdaA$analysis$oDG$data.groups$modelIII)
modelIII.Estimates<- data.frame(modelIII.Estimates$param.free)
modelIII.CI<- data.frame(hzar.getLLCutParam(mkn$AdaA$analysis$oDG$data.groups$modelIII,
                                            names(mkn$AdaA$analysis$oDG$data.groups$modelIII$data.param)))
modelIII.T<- data.frame(Model= "modelIII", SNP= i)
modelIII.T<- cbind(modelIII.T,modelIII.Estimates,modelIII.CI)
modelIII.Table<- rbind(modelIII.Table,modelIII.T)
assign(paste("modeI_results_SNP_",i,sep=""),mkn$AdaA$analysis$oDG$data.groups$modelIII)
Data_modelIII<- get(paste("modeI_results_SNP_",i,sep=""))
saveRDS(Data_modelIII, paste("./ResultModels/",TRANSECT,"modelIII_results_SNP_",i,sep=""), ascii=TRUE)
}

# Write out all summary tables
write.table(modelI.Convergence.TABLE,paste("./ResultTables/",TRANSECT,"ModelI_convergence_table.txt",sep=""),row.names = F,col.names = F,quote = F,sep="\t")
write.table(modelII.Convergence.TABLE,paste("./ResultTables/",TRANSECT,"ModelII_convergence_table.txt",sep=""),row.names = F,col.names = F,quote = F,sep="\t")
write.table(modelIII.Convergence.TABLE,paste("./ResultTables/",TRANSECT,"ModelIII_convergence_table.txt",sep=""),row.names = F,col.names = F,quote = F,sep="\t")
write.table(FourModels.Table,paste("./ResultTables/",TRANSECT,"Check4ModelIncluded.txt",sep=""),row.names = F,col.names = F,quote = F,sep="\t")
write.table(AIC.Table,paste("./ResultTables/",TRANSECT,"AIC_table.txt",sep=""),row.names = F,col.names = T,quote = F,sep="\t")
write.table(modelI.Table,paste("./ResultTables/",TRANSECT,"ModelI.RESULT.Table.txt",sep=""),row.names = F,col.names = T,quote = F,sep="\t")
write.table(modelII.Table,paste("./ResultTables/",TRANSECT,"ModelII.RESULT.Table.txt",sep=""),row.names = F,col.names = T,quote = F,sep="\t")
write.table(modelIII.Table,paste("./ResultTables/",TRANSECT,"ModelIII.RESULT.Table.txt",sep=""),row.names = F,col.names = T,quote = F,sep="\t")
