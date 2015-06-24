#Elastic net analysis with entries 9 and 10 removed
#Data set up
setwd("~/GitHub/Metagenomic-Analysis-Part-2")
source("TwoStage_Package_Code.R")
#Grab all cleaned data
path<-as.list(dir(pattern=".*cleaned.*"))
my.data<-lapply(path, function(x) read.csv(x, header=T))
#Remove relevant rows, row names, and the city column so the phenotype file will run properly
factors<-read.csv("phyllo_factors.csv",header=T)[-c(1,2, 9:12),]
rownames(factors)<-factors$city<-NULL

#Function for removing the PHYLLO09 and PHYLLO10 entries
removeExtras<-function(data){
  newdata<-data[,-c(2:3)]
  return(newdata)
}
#NOTE: These filenames are used in this order 
#since this is their order in the original folder (sorted alphabetically)
filenames<-list("ECID", "GOFuncID","GOProcID","PFamID","Phyllo_Genus","TIGRFamID")
my.new.data<-lapply(my.data, removeExtras)
#Dump extraneous variables for saving a bit of memory
remove(my.data)
remove(path)
for (i in 1:length(my.new.data)){
  try(TwoStage_Package(my.new.data[[i]], factors, paste(filenames[[i]],"sigtest","csv", sep="."), 1))
  try(TwoStage_Package(my.new.data[[i]], factors, paste(filenames[[i]],"sigtest-method2","csv",sep="."), 2))
}