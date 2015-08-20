#Calculating the trimmed means, sorted by treatment and separated by city
setwd("~/GitHub/Metagenomic-Analysis-Part-2")
path<-as.list(dir(pattern=".*cleaned.*"))
my.data<-lapply(path, function(x) read.csv(x, header=T))
#Remove relevant rows, row names, and the city column so the phenotype file will run properly
factors<-read.csv("phyllo_factors.csv",header=T)
factors.treat<-factors[-c(1,2),]
rownames(factors.treat)<-rownames(factors)<-factors.treat$city<-NULL

#Trimmed mean function (from the edgeR paper)
#Function borrowed from the Two_Stage Package.R file
TMMNorm = function(X, Y, TMM.option){
  require(edgeR)
  factors = calcNormFactors(X,method="TMM") #Calculate normaization factors
  if(TMM.option==1){
    eff.lib.size = colSums(X)*factors;
    ref.lib.size = mean(eff.lib.size); #Use the mean of the effective library sizes as a reference library size 
    X.output = sweep(X,MARGIN=2,eff.lib.size,"/")*ref.lib.size; #Normalized read counts
  } else if(TMM.option==2){ #Regenerate counts with a common dispersion.
    group = as.factor(Y);
    X.tmp = DGEList(counts=X, group=group)
    X.tmp = calcNormFactors(X.tmp)
    X.tmp = estimateCommonDisp(X.tmp)
    X.output = as.matrix(X.tmp$pseudo.counts)
  } else {
    stop("TMM.option must be either 1 or 2!")
  }
  return(X.output)
}

#Function for removing the PHYLLO09 and PHYLLO10 entries
removeExtras<-function(data){
  newdata<-data[,-c(2:3)]
  return(newdata)
}
treat.data<-lapply(my.data, removeExtras)
for (i in 1:length(my.data)){
  treat.data[[i]]<-treat.data[[i]][,c(1,2,4,6,8,10,13,14,3,5,7,9,11,12,15:16)]
}

#Save annotation to be bound to the data frame later
rows.list<-lapply(my.data, function(x) return(x$X))
treat.data<-lapply(treat.data, function(x) return (x[,-1]))
trimmed.means.norm<-new.data<-list()
#Run the trimmed mean function with the normalization method
for (i in 1:length(my.data)){
  trimmed.means.norm[[i]]<-TMMNorm(treat.data[[i]], factors, 1)
  new.data[[i]]<-cbind(rows.list[[i]], trimmed.means.norm[[i]])
  write.csv(new.data[[i]], file=paste(unlist(strsplit(path[[i]], "-"))[1], "trimmed-means", "csv", sep="."), row.names=F)
}