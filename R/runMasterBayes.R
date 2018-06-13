
setwd("~/Projects/parentage_review/results")
library(MasterBayes)
source("../R/parentageConverteR.R")
gty.files<-list.files(pattern="genotypes.txt")

mb<-lapply(gty.files[6:length(gty.files)],function(gfile){
  print(gfile)
  G<-cervus2coancestryG(gfile)
  GdP<-GdataPed(G=G, categories=NULL)
  P<-cervus2coancestryP(G$id,sexes=c(Males="MAL",Females="FEM",Offspring="OFF"))
  colnames(P)[colnames(P)=="IDs"]<-"ids"
  #assuming we know the mothers - need to get family IDs
  mom.off<-read.delim(gfile)[,1:3]
  mom.off<-mom.off[which(!is.na(mom.off$Mom)),1:2]
  P<-merge(P,mom.off,by.x="ids",by.y="ID",all.x = TRUE) #assign moms
  res1<-expression(varPed(x="Mom",gender="Female",relational="OFFSPRING",restrict="=="))#offspring match to females
  res2<-expression(varPed(x="offspring",gender=NULL,relational=FALSE,restrict=0)) #remove offspring from parental pop
  PdP<-PdataPed(formula=list(res1,res2),data=P,USsire=TRUE)
  model<-MCMCped(PdP=PdP, GdP=GdP, verbose=TRUE,burnin=50000,nitt=55000)
  saveRDS(model,paste("MasterBayes/",gsub("genotypes.txt","mod.RDS",gfile),sep=""))
})