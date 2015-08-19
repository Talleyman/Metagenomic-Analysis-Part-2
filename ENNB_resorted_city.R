#Resort the metagenomic data, this time by city
#First, get the data set up
setwd("~/GitHub/Metagenomic-Analysis-Part-2")
source("TwoStage_Package_Code.R")
#Grab all cleaned data
path<-as.list(dir(pattern=".*cleaned.*"))
my.data<-lapply(path, function(x) read.csv(x, header=T))
factors<-read.csv("phyllo_factors.csv",header=T)
rownames(factors)<-factors$treatment<-NULL

for (i in 1:length(my.data)){
  my.data[[i]][,19]<-NA
  colnames(my.data[[i]])[19]<-"PHYLLO30"
  my.data[[i]]$PHYLLO30<-round((my.data[[i]]$PHYLLO28+my.data[[i]]$PHYLLO29)/2)
}

#Impute the new PHYLLO30 sample into the phenotype file
factors[18,]<-NA
factors[,1]<-sapply(factors[,1],as.character)
factors[18,]<-c("PHYLLO30", "HF")
factors[,1:2]<-sapply(factors[,1:2],as.factor)

#Resort the samples according to which samples come from which city for the analysis
for (i in 1:length(my.data)){
  my.data[[i]]<-my.data[[i]][,c(1,2,15,16,17:19,3,8,12,9,13,14,4,6,10,5,7,11)]
}

#NOTE: These filenames are used in this order 
#since this is their order in the original folder (sorted alphabetically)
filenames<-list("ECID", "GOFuncID","GOProcID","PFamID","Phyllo_Genus","TIGRFamID")
for (i in 1:length(my.data)){
  try(TwoStage_Package(my.data[[i]], factors, paste(filenames[[i]],"citysort","sigtest","csv", sep="."), 1))
  try(TwoStage_Package(my.data[[i]], factors, paste(filenames[[i]],"citysort","sigtest-method2","csv",sep="."), 2))
}