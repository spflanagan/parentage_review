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


#' Convert a CERVUS file into input for FaMoz
#' @param gty A data.frame of genotypes in CERVUS format
#' @param filename The name of the file to export the FaMoz file to.
#' @param first.allele The first column with an allele (locus 1, alele 1). Default is 4.
#' @param miss The character specifying missing data. Default is 0.  
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @export
cervus2famoz<-function(gty,filename,first.allele=4,miss=0,sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
  nloci<-(ncol(gty)-(first.allele-1))/2
  #write loci
  suppressWarnings(write.table(nloci,filename,quote=FALSE,col.names = FALSE,row.names = FALSE))
  frqs<-lapply(seq(first.allele,ncol(gty),2),function(loc1,gty){
    locus<-gty[,c(loc1,loc1+1)]
    frq<-unlist(locus[locus!=miss])
    dat<-data.frame(t(data.frame(table(frq)/sum(table(frq)))$Freq))
    row.names(dat)<-ncol(dat)
    suppressWarnings(write.table(dat,filename,append = TRUE,row.names=TRUE,col.names = FALSE,quote=FALSE))
  },gty=gty)
  #how many parents?
  ns<-c(length(c(grep(sexes["Males"],gty$ID),grep(sexes["Females"],gty$ID))),
        length(grep(sexes["Offspring"],gty$ID)))
  suppressWarnings(write.table(t(ns),filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
  suppressWarnings(write.table('',filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
  #parent genotypes
  gty[gty==miss]<--5
  parents<-as.matrix(gty[c(grep(sexes["Males"],gty$ID),grep(sexes["Females"],gty$ID)),
                         first.allele:ncol(gty)])
  rownames(parents)<-NULL
  suppressWarnings(write.table(parents,filename,append = TRUE,row.names=TRUE,col.names = FALSE,quote=FALSE))
  suppressWarnings(write.table('',filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
  #offspring genotypes
  offspring<-as.matrix(gty[grep(sexes["Offspring"],gty$ID),first.allele:ncol(gty)])
  rownames(offspring)<-NULL
  suppressWarnings(write.table(offspring,filename,append = TRUE,row.names=TRUE,col.names = FALSE,quote=FALSE))
  
}

#' Convert a CERVUS file into input for SOLOMON
#' @param gty A data.frame of genotypes in CERVUS format
#' @param filename The name of the file to export the SOLOMON files to.
#' @param first.allele The first column with an allele (locus 1, alele 1). Default is 4.
#' @param miss The character specifying missing data. Default is 0.  
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @param output A boolean variable designating whether adults should be output together (TRUE) or separately as moms and dads (FALSE)
#' @export
cervus2solomon<-function(gty,filename,first.allele=4,miss=0,sexes=c(Males="MAL",Females="FEM",Offspring="OFF"),
                         together=FALSE,knownDads=FALSE,knownMoms=FALSE){
  if(miss != 0){
    gty[gty==miss]<-0
  }
  if(isTRUE(together)){
    #parent genotypes
    parents<-gty[c(grep(sexes["Males"],gty$ID),grep(sexes["Females"],gty$ID)),c(1,first.allele:ncol(gty))]
    rownames(parents)<-NULL
    write.table(parents,paste(filename,"parents.txt",sep=""),row.names=FALSE,col.names = TRUE,quote=FALSE,sep='\t')
  }else{ #they're separate
    #dad genotypes
    dads<-gty[grep(sexes["Males"],gty$ID),c(1,first.allele:ncol(gty))]
    rownames(dads)<-NULL
    write.table(dads,paste(filename,"dads.txt",sep=""),row.names=FALSE,col.names = TRUE,quote=FALSE,sep='\t')
    #mom genotypes
    moms<-gty[grep(sexes["Females"],gty$ID),c(1,first.allele:ncol(gty))]
    rownames(moms)<-NULL
    write.table(moms,paste(filename,"moms.txt",sep=""),row.names=FALSE,col.names = TRUE,quote=FALSE,sep='\t')
  }
  #offspring genotypes
  offspring<-gty[grep(sexes["Offspring"],gty$ID),c(1,first.allele:ncol(gty))]
  rownames(offspring)<-NULL
  write.table(offspring,paste(filename,"offspring.txt",sep=""),row.names=FALSE,col.names = TRUE,quote=FALSE,sep='\t')
}

#' Convert a data.frame formatted for CERVUS to input for sequoia
#' @param gt.name The file name of the genotypes file
#' @param sep An optional character that is the separator for the genotype file (default is tab delimited)
#' @param start.col The column where the first allele is found. Default is 2.
#' @param header A boolean variable indicating whether the file has a header row (default is TRUE)
#' @return gty A matrix with all of the genotypes in 0,1,2 format with -9 as missing values
#' @export
cervus2colonyG<-function(gt.name,sep='\t',start.col=4,header = TRUE){
  gt<-read.delim(gt.name,sep=sep,header=header)
  gty<-gt[,c(1,start.col:ncol(gt))]
  colnames(gty)[1]<-"id"
  return(gty)
}

#' Convert a data.frame formatted for CERVUS to input for sequoia
#' @param ids A vector or data.frame with all of the individual IDs. If data.frame is provided, then all of the columns are retained in the final dataset, and the IDs need to be in the first column.
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @return A data.frame with columns ID, offspring, and sex
#' @export
cervus2colonyP<-function(ids,sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
  if(class(ids) == "data.frame"){
    dat<-ids
    colnames(dat)[1]<-"IDs"
  }else{
    dat<-data.frame(IDs=ids,stringsAsFactors = FALSE)
  }
  dat$sex<-NA
  dat$sex[grep(sexes["Males"],dat$IDs)]<-"Male"
  dat$sex[grep(sexes["Females"],dat$IDs)]<-"Female"
  dat$offspring<-0
  dat$offspring[grep(sexes["Offspring"],dat$IDs)]<-1
  return(dat)
}

#' Convert a CERVUS file into PLINK format
#' @param gty A data.frame of genotypes in CERVUS format
#' @export
cervus2plink<-function(file,out.dir="./",first.allele=4,sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
  gty<-read.delim(file)
  snp.names<-colnames(gty)[first.allele:ncol(gty)]
  gty$FamID<-gty$ID
  gty$sex<-NA
  gty$sex[grep(sexes["Males"],gty$sex)]<-1
  gty$sex[grep(sexes["Females"],gty$sex)]<-2
  gty$Phenotype<--9
  ped.name<-paste(out.dir,gsub("_genotypes.txt",".ped",file),sep="")
  map.name<-paste(out.dir,gsub("_genotypes.txt",".map",file),sep="")
  ped<-gty[,c("FamID","ID","Dad","Mom","sex","Phenotype",snp.names)]
  map<-data.frame(Chr=0,snp=snp.names[seq(first.allele,ncol(gty),2)],
                  distance=0,bp=0)
  write.table(ped,ped.name,row.names = FALSE,col.names=TRUE,quote=TRUE,sep='\t')
  write.table(map,map.name,row.names = FALSE,col.names=FALSE,quote=TRUE,sep='\t')
  return(ped)
}

#' Convert a CERVUS file into input for PedApp
#' @param gty A data.frame of genotypes in CERVUS format
#' @param filename The name of the file to export the FaMoz file to.
#' @param first.allele The first column with an allele (locus 1, alele 1). Default is 4.
#' @export
cervus2pedapp<-function(gty,filename,first.allele=4){
  #Columns 1-2 is pedigree structure matrix
  #Column 3 is age (generation)
  #Column 4 is output flag (for 1 allows output)
  #Columns 5 - x are genotype data
  #Comma delimited
  last.col<-ncol(gty)
  gty$age<-gsub("(\\w{3}).*","\\1",gty$ID)
  gty$age[gty$age=="OFF"]<-0
  gty$age[gty$age!=0]<-1
  gty$sex<-gsub("(\\w{3}).*","\\1",gty$ID)
  gty$sex[gty$sex=="OFF"]<-0
  gty$output<-1
  new.gty<-cbind(gty[,c("ID","age","sex","output")],gty[,first.allele:last.col])
  write.table(new.gty,filename,row.names = FALSE,col.names=FALSE,
              sep=",",quote=FALSE)
}


#' Convert a CERVUS file into input for FaMoz
#' @param gty A data.frame of genotypes in CERVUS format
#' @param filename The name of the file to export the FaMoz file to.
#' @param first.allele The first column with an allele (locus 1, alele 1). Default is 4.
#' @param miss The character specifying missing data. Default is 0. 
#' @param type Default is "snps". For microsatellites it should be ‘M’ followed by the length of repeat-unit.
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @export
cervus2parfex<-function(gty,filename,first.allele=4,miss=0,type="snps",sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
  #get marker names
  marker.names<-c("Marker",colnames(gty)[seq(first.allele,ncol(gty))])
  marker.names[grep("b",marker.names)]<-""
  suppressWarnings(write.table(as.data.frame(t(marker.names)),filename,quote=FALSE,col.names = FALSE,row.names = FALSE,sep='\t'))
  #marker types
  marker.types<-marker.names
  marker.types[1]<-"MarkerType"
  marker.types[seq(2,length(marker.types),2)]<-type
  suppressWarnings(write.table(as.data.frame(t(marker.types)),filename,quote=FALSE,col.names = FALSE,row.names = FALSE,append=TRUE,sep='\t'))
  #fix the genotypes
  gty[gty==miss]<-"?"
  bps<-sample(c("A","C","G","T"),2)
  gty[gty==1]<-bps[1]
  gty[gty==2]<-bps[2]
  #Offspring
  suppressWarnings(write.table("Offspring",filename,quote=FALSE,col.names = FALSE,row.names = FALSE,append=TRUE,sep='\t'))
  offspring<-as.matrix(gty[grep(sexes["Offspring"],gty$ID),first.allele:ncol(gty)])
  rownames(offspring)<-gty[grep(sexes["Offspring"],gty$ID),"ID"]
  suppressWarnings(write.table(offspring,filename,append = TRUE,row.names=TRUE,col.names = FALSE,quote=FALSE,sep='\t'))
  #Parents
  suppressWarnings(write.table("Parents",filename,quote=FALSE,col.names = FALSE,row.names = FALSE,append=TRUE,sep='\t'))
  parents<-as.matrix(gty[c(grep(sexes["Males"],gty$ID),grep(sexes["Females"],gty$ID)),
                         first.allele:ncol(gty)])
  rownames(parents)<-gty[c(grep(sexes["Males"],gty$ID),grep(sexes["Females"],gty$ID)),"ID"]
  suppressWarnings(write.table(parents,filename,append = TRUE,row.names=TRUE,col.names = FALSE,quote=FALSE,sep='\t'))
}
