#' @author Sarah P. Flanagan


#' Convert a data.frame formatted for CERVUS to input for sequoia
#' @param gt.name The file name of the genotypes file
#' @param sep An optional character that is the separator for the genotype file (default is tab delimited)
#' @param start.col The column where the first allele is found. Default is 2.
#' @param header A boolean variable indicating whether the file has a header row (default is TRUE)
#' @return GenoM A matrix with all of the genotypes in 0,1,2 format with -9 as missing values
#' @export
#' @example
#' GenoM<-cervus2sequoia(gt.name="parentsim_genotypes.txt",start.col=4)
cervus2sequoia<-function(gt.name,sep='\t',start.col=2,header = TRUE){
  gt<-read.delim(gt.name,sep=sep,header=header)
  GenoM<-as.matrix(do.call(cbind,lapply(seq(start.col,ncol(gt),2),function(coln,gt){
    dat<-gt[,c(coln,coln+1)]
    ref<-names(which.max(table(apply(dat,2,factor))))
    out<-apply(dat,1,function(ind,ref){
      if(ind[1]==ind[2]){
        if(ind[2]==ref){ out<-2 }
        if(ind[2]==0){ out<--9 }
        if(ind[2]!=ref & ind[2] != 0 ){ out<-0 }
      }else{
        out<-1
      }
      return(out)
    },ref)
  },gt)))
  rownames(GenoM)<-gt[,1]
  return(GenoM)
}

#' Generate a Life history dataframe for sequoia from IDs
#' @param GenoM.ID A list of individual IDs
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals by sex and age
#' @return LifeHistDat A data.frame with ID, sex, and birth year (BY)
#' @export
#' @example
#' LifeHistData<-generate_seq_lifehist(rownames(GenoM),sexes=c(Males="MAL",Females="FEM",Offspring="OFF"))
generate_seq_lifehist<-function(GenoM.ID,sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
  LifeHistDat<-data.frame(ID=GenoM.ID,sex=NA,BY=1)
  LifeHistDat$sex[grep(sexes["Males"],LifeHistDat$ID)]<-as.numeric(2)
  LifeHistDat$sex[grep(sexes["Females"],LifeHistDat$ID)]<-as.numeric(1)
  LifeHistDat$BY[grep(sexes["Offspring"],LifeHistDat$ID)]<-2
  return(LifeHistDat)
}
