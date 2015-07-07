setwd("~/GitHub/Metagenomic-Analysis-Part-2")
#Grab all cleaned data
path<-as.list(dir(pattern=".*cleaned.*"))
my.data<-lapply(path, function(x) read.csv(x, header=T))
#Remove relevant rows, row names, and the city column so the phenotype file will run properly
factors<-read.csv("phyllo_factors.csv",header=T)[-c(1,2, 9:12),]
rownames(factors)<-factors$city<-NULL

#Trimmed mean function (from the edgeR paper)
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
#Save annotation to be bound to the data frame later
rows.list<-lapply(my.data, function(x) return(x$X))
my.data<-lapply(my.data, function(x) return (x[,-1]))
trimmed.means.norm<-new.data<-list()
#Run the trimmed mean function with the normalization method
for (i in 1:length(my.data)){
  trimmed.means.norm[[i]]<-TMMNorm(my.data[[i]], factors, 1)
  new.data[[i]]<-cbind(rows.list[[i]], trimmed.means.norm[[i]])
  write.csv(new.data[[i]], file=paste(unlist(strsplit(path[[i]], "-"))[1], "trimmed-means", "csv", sep="."), row.names=F)
}