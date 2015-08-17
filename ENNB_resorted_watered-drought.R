#Data set up
setwd("~/GitHub/Metagenomic-Analysis-Part-2")
source("TwoStage_Package_Code.R")
#Grab all cleaned data
path<-as.list(dir(pattern=".*cleaned.*"))
my.data<-lapply(path, function(x) read.csv(x, header=T))
#Remove relevant rows, row names, and the city column so the phenotype file will run properly
factors<-read.csv("phyllo_factors.csv",header=T)[-c(1,2, 9:12),]
rownames(factors)<-factors$city<-NULL
removeExtras<-function(data){
  newdata<-data[,-c(2:3)]
  return(newdata)
}
my.data<-lapply(my.data, removeExtras)
#Resort data according by watered and drought groups, respectively
for (i in 1:length(my.data)){
  my.data[[i]]<-my.data[[i]][,c(1,2,4,6,8,10,13,14,3,5,7,9,11,12,15:16)]
}
filenames<-list("ECID", "GOFuncID","GOProcID","PFamID","Phyllo_Genus","TIGRFamID")
remove(path)
for (i in 1:length(my.data)){
  try(TwoStage_Package(my.data[[i]], factors, paste(filenames[[i]],"resorted_sigtest","csv", sep="."), 1))
  try(TwoStage_Package(my.data[[i]], factors, paste(filenames[[i]],"resorted_sigtest-method2","csv",sep="."), 2))
}