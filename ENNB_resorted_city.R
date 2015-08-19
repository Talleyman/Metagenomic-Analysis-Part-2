#Resort the metagenomic data, this time by city
#First, get the data set up
setwd("~/GitHub/Metagenomic-Analysis-Part-2")
source("TwoStage_Package_Code.R")
#Grab all cleaned data
path<-as.list(dir(pattern=".*cleaned.*"))
my.data<-lapply(path, function(x) read.csv(x, header=T))
factors<-read.csv("phyllo_factors.csv",header=T)

for (i in 1:length(my.data)){
  my.data[[i]][,19]<-NA
  colnames(my.data[[i]])[19]<-"PHYLLO30"
  my.data[[i]]$PHYLLO30<-round((my.data[[i]]$PHYLLO28+my.data[[i]]$PHYLLO29)/2)
}

#Impute the new PHYLLO30 sample into the phenotype file
factors[18,]<-NA
factors[,1]<-sapply(factors[,1],as.character)
factors[18,]<-c("PHYLLO30", "HF", "drought")
factors[,1:3]<-sapply(factors[,1:3],as.factor)

rownames(factors)<-NULL
CA.obs<-which(factors$city=="CA")
CA.group<-factors[CA.obs,]
HF.obs<-which(factors$city=="HF")
HF.group<-factors[HF.obs,]
DE.obs<-which(factors$city=="DE")
DE.group<-factors[DE.obs,]
rownames(CA.group)<-CA.group$city<-rownames(DE.group)<-DE.group$city<-rownames(HF.group)<-HF.group$city<-NULL

CA.list<-DE.list<-HF.list<-list()
for (i in 1:length(my.data)){
  CA.list[[i]]<-my.data[[i]][,c(1,CA.obs+1)]
  DE.list[[i]]<-my.data[[i]][,c(1,DE.obs+1)]
  HF.list[[i]]<-my.data[[i]][,c(1,HF.obs+1)]
}

#rearrange based on treatment (HF is conveniently already in order, so no need to rearrange it)
DE.list<-lapply(DE.list, function(x) x[,c(1,2,3,5,4,6,7)])
CA.list<-lapply(CA.list, function(x) x[,c(1,2,4,6,3,5,7)])

#NOTE: These filenames are used in this order 
#since this is their order in the original folder (sorted alphabetically)
filenames<-list("ECID", "GOFuncID","GOProcID","PFamID","Phyllo_Genus","TIGRFamID")
for (i in 1:length(my.data)){
  try(TwoStage_Package(HF.list[[i]], HF.group, paste(filenames[[i]],"analysis-HF","sigtest","csv", sep="."), 1))
  try(TwoStage_Package(HF.list[[i]], HF.group, paste(filenames[[i]],"analysis-HF","sigtest-method2","csv",sep="."), 2))
  try(TwoStage_Package(DE.list[[i]], DE.group, paste(filenames[[i]],"analysis-DE","sigtest","csv", sep="."), 1))
  try(TwoStage_Package(DE.list[[i]], DE.group, paste(filenames[[i]],"analysis-DE","sigtest-method2","csv",sep="."), 2))
  try(TwoStage_Package(CA.list[[i]], CA.group, paste(filenames[[i]],"analysis-CA","sigtest","csv", sep="."), 1))
  try(TwoStage_Package(CA.list[[i]], CA.group, paste(filenames[[i]],"analysis-CA","sigtest-method2","csv",sep="."), 2))
}