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

#' Create Marker Types and Error Rates file for colony
#' @param out.name The file name of the output file
#' @param marker.names A list of marker names
#' @param codom 1 if the markers are codominant, 2 if they are dominant. Either one value for all loci or a vector containing values for each.
#' @param ado Allelic dropout rates. Either one value for all loci or a vector containing values for each.
#' @param ge Genotyping error rates. Either one value for all loci or a vector containing values for each.
#' @export
colonyMarkerTypesError<-function(out.name,marker.names,codom=0,ado=0.01,ge=0.001){
  if(length(ado)==1){ #same info for all markers
    ado<-rep(ado,length(marker.names))
  }
  if(length(codom)==1){ #same info for all markers
    codom<-rep(codom,length(marker.names))
  }
  if(length(ge)==1){ #same info for all markers
    ge<-rep(ge,length(marker.names))
  }
  write.table(data.frame(rbind(marker.names,codom,ado,ge)),
              out.name,col.names = FALSE,row.names=FALSE,quote=FALSE)
  return(invisible(data.frame(rbind(marker.names,codom,ado,ge))))
}

#' Convert a data.frame formatted for CERVUS to input for colony - the offspring genotypes
#' @param gty The genotypes data frame
#' @param out.name The output name
#' @param sub.string A common string in IDs to mark individuals as part of the subset
#' @param first.allele The column where the first allele is found. Default is 4.
#' @return off.gty A matrix with all of the genotypes in 0,1,2 format with -9 as missing values
#' @export
cervus2colonySUB<-function(gty,out.name,sub.string,first.allele=4){
  sub.gty<-gty[grep(sub.string,gty[,1]),c(1,first.allele:ncol(gty))]
  write.table(sub.gty,out.name,
              col.names=FALSE,row.names=FALSE,quote=FALSE)
  return(invisible(sub.gty))
}


#' Convert a data.frame formatted for CERVUS to input for coancestry
#' @param gty.name The file name of the genotypes file
#' @param dir The directory for output files
#' @param first.allele The column where the first allele is found. Default is 4.
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @param known.mom The column with known moms. Use NA to not include.
#' @param known.dad The column with known dads. Use NA to not include.
#' @param codom 1 if the markers are codominant, 2 if they are dominant. Either one value for all loci or a vector containing values for each.
#' @param ado Allelic dropout rates. Either one value for all loci or a vector containing values for each.
#' @param ge Genotyping error rates. Either one value for all loci or a vector containing values for each.
#' @return gty A matrix with all of the genotypes in 0,1,2 format with -9 as missing values
#' @export
cervus2colony<-function(gty.file,dir="./",first.allele=4, sexes=c(Males="MAL",Females="FEM",Offspring="OFF"),
                        known.mom=2,known.dad=1,codom=0,ado=0.01,ge=0.001){
  gty<-read.delim(gty.file)
  basename<-paste(dir,gsub("_genotypes.txt","",gty.file),sep="")
  # Marker Types and Error Rates
  colonyMarkerTypesError(paste(basename,"_TypesError.txt",sep=""),
                         marker.names =
                           colnames(gty)[seq(first.allele,ncol(gty),2)],
                         codom = codom,ado=ado,ge=ge)
  # Genotypes
  cervus2colonySUB(gty = gty,out.name = paste(basename,"_offspring.txt",sep=""),sexes["Offspring"],
                   first.allele=first.allele)
  cervus2colonySUB(gty,paste(basename,"_males.txt",sep=""),sexes["Males"],
                   first.allele=first.allele)
  cervus2colonySUB(gty,paste(basename,"_females.txt",sep=""),sexes["Females"],
                   first.allele=first.allele)
  # Known parent pairs
  if(!is.na(known.mom)){
    known.pairs<-data.frame(Off=gty$ID[grep(sexes["Offspring"],gty$ID)],
                       Moms=gty[grep(sexes["Offspring"],gty$ID),known.mom],
                       stringsAsFactors = FALSE)
    known.pairs<-known.pairs[known.pairs$Moms %in% gty$ID,]
    if (file.exists(paste(basename,"_knownmoms.txt"))) file.remove(paste(basename,"_knownmoms.txt"))
    invisible(do.call(rbind,lapply(unique(known.pairs$Moms),function(mom,allbabes,out.name){
      babies<-allbabes[allbabes[,2]==mom,1]
      write.table(t(c(as.character(mom),as.character(babies))),
                  out.name,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
    },allbabes=known.pairs,out.name=paste(basename,"_knownmoms.txt"))))
  }
  if(!is.na(known.dad)){
    known.pairs<-data.frame(Off=gty$ID[grep(sexes["Offspring"],gty$ID)],
                            Dads=gty[grep(sexes["Offspring"],gty$ID),known.dad],
                            stringsAsFactors = FALSE)
    known.pairs<-known.pairs[known.pairs$Dads %in% gty$ID,]
    if (file.exists(paste(basename,"_knowndads.txt"))) file.remove(paste(basename,"_knowndads.txt"))
    invisible(do.call(rbind,lapply(unique(known.pairs$Dads),function(dad,allbabes,out.name){
      babies<-allbabes[allbabes[,2]==dad,1]
      write.table(t(c(as.character(dad),as.character(babies))),
                  out.name,col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
    },allbabes=known.pairs,out.name=paste(basename,"_knowndads.txt"))))
  }
  invisible(gty.file)
}


#' Convert a data.frame formatted for CERVUS to input for coancestry
#' @param gt.name The file name of the genotypes file
#' @param sep An optional character that is the separator for the genotype file (default is tab delimited)
#' @param start.col The column where the first allele is found. Default is 2.
#' @param header A boolean variable indicating whether the file has a header row (default is TRUE)
#' @return gty A matrix with all of the genotypes in 0,1,2 format with -9 as missing values
#' @example Gs<-lapply(gty.files,function(file){
#'outname<-paste("colony/",file,sep="")
#'G<-cervus2colonyG(file)
#'write.table(G,outname,row.names = FALSE,col.names = FALSE,sep='\t',quote=FALSE)
#'
#'P<-cervus2colonyP(G,sexes=c(Males="MAL",Females="FEM",Offspring="OFF"))
#'p.outname<-paste("colony/",gsub("genotypes","phenotypes",file),sep="")
#'write.table(P,p.outname,row.names = FALSE,col.names = TRUE,sep='\t',quote=FALSE)
#'return(G)
#'})
#' @export
cervus2coancestryG<-function(gt.name,sep='\t',start.col=4,header = TRUE){
  gt<-read.delim(gt.name,sep=sep,header=header)
  gty<-gt[,c(1,start.col:ncol(gt))]
  colnames(gty)[1]<-"id"
  return(gty)
}

#' Convert a data.frame formatted for CERVUS to input for coancestry
#' @param ids A vector or data.frame with all of the individual IDs. If data.frame is provided, then all of the columns are retained in the final dataset, and the IDs need to be in the first column.
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @return A data.frame with columns ID, offspring, and sex
#' @export
cervus2coancestryP<-function(ids,sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
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

#' Convert a CERVUS file into PLINK v1.7 format
#' @param file A data.frame of genotypes in CERVUS format
#' @param out.dir The directory where the output files belong (default is current director)
#' @param first.allele The location of the first genotype column
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @export
cervus2ped<-function(file,out.dir="./",first.allele=4,sexes=c(Males="MAL",Females="FEM",Offspring="OFF")){
  if(substr(out.dir,nchar(out.dir),nchar(out.dir))!="/"){
    out.dir<-paste(out.dir,"/",sep="")
  }
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
  map<-data.frame(Chr=0,snp=snp.names[seq(1,length(snp.names),2)],
                  distance=0,bp=0)
  write.table(ped,ped.name,row.names = FALSE,col.names=FALSE,quote=FALSE,sep='\t')
  write.table(map,map.name,row.names = FALSE,col.names=FALSE,quote=FALSE,sep='\t')
  return(ped)
}

#' Convert a CERVUS file into PLINK v1.9 format
#' @param file A data.frame of genotypes in CERVUS format
#' @param out.dir The directory where the output files belong (default is current director)
#' @param first.allele The location of the first genotype column
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males, females and offspring
#' @export
cervus2tped<-function(file,out.dir="./",first.allele=4,sexes=c(Males="MAL",Females="FEM",Offspring="OFF"),
                      known.mom=NA,known.dad=NA){
  if(substr(out.dir,nchar(out.dir),nchar(out.dir))!="/"){
    out.dir<-paste(out.dir,"/",sep="")
  }
  gty<-read.delim(file)
  end.snps<-ncol(gty)
  snps<-colnames(gty)[first.allele:end.snps]
  gty$FamID<-0 #no family names
  if(!is.na(known.mom)){ #then family names will come from moms
    gty$FamID[grep(sexes["Females"],gty[,known.mom])]<-paste("Fam",gsub(sexes["Females"],"",gty[grep(sexes["Females"],gty[,known.mom]),known.mom]),sep="")
    gty$FamID[grep(sexes["Females"],gty$ID)]<-paste("Fam",gsub(sexes["Females"],"",gty$ID[grep(sexes["Females"],gty$ID)]),sep="")
  }
  if(!is.na(known.dad)){ #then family names will come from moms
    gty$FamID[grep(sexes["Males"],gty[,known.dad])]<-paste("Fam",gsub(sexes["Males"],"",gty[grep(sexes["Males"],gty[,known.dad]),known.dad]),sep="")
    gty$FamID[grep(sexes["Males"],gty$ID)]<-paste("Fam",gsub(sexes["Males"],"",gty$ID[grep(sexes["Males"],gty$ID)]),sep="")
  }
  gty$sex<--9
  gty$sex[grep(sexes["Males"],gty$ID)]<-1
  gty$sex[grep(sexes["Females"],gty$ID)]<-2
  tped.name<-paste(out.dir,gsub("_genotypes.txt",".tped",file),sep="")
  tfam.name<-paste(out.dir,gsub("_genotypes.txt",".tfam",file),sep="")
  gty[,snps][gty[,snps]==1]<-"A"
  gty[,snps][gty[,snps]==2]<-"T"
  tped<-data.frame(cbind(Chr=1,snp=snp.names[seq(1,length(snps),2)],
                   distance=seq(1,length(snps),2)/length(snps),bp=seq(1,length(snps),2),
                   matrix(as.character(NA),nrow=length(snps)/2,ncol=2*nrow(gty))),
                   stringsAsFactors = FALSE)
  ## Tried to vectorize and failed
  # m<-mapply(function(i,cnt,gty,tped){
  #   mapply(function(ii,ind.cnt,gtyi,tpedcnt){
  #     tpedcnt[,c(ind.cnt,ind.cnt+1)]<-gtyi[ii,]
  #   },1:nrow(gty),seq(5,ncol(tped),2),
  #   MoreArgs = list(gtyi=gty[,c(i,i+1)],tpedcnt=tped[cnt,]))
  # },seq(first.allele,end.snps,2),seq(1:length(snps)/2),
  #        MoreArgs=list(gty=gty,tped=tped))
  cnt<-1
  for(i in seq(first.allele,end.snps,2)){#loop through snps
    ind.cnt<-5 #where the first genotype goes
    for(ii in 1:nrow(gty)){ #loop through individuals
      tped[cnt,c(ind.cnt,ind.cnt+1)]<-gty[ii,c(i,i+1)]
      ind.cnt<-ind.cnt+2
    }
    cnt<-cnt+1
  }
  tfam<-data.frame(gty$FamID,gty$ID,Dad=as.character(gty$Dad),Mom=as.character(gty$Mom),gty$sex,
                   stringsAsFactors = FALSE)
  tfam$Dad[is.na(tfam$Dad)]<--9
  tfam$Mom[is.na(tfam$Mom)]<--9
  write.table(tped,tped.name,row.names = FALSE,col.names=FALSE,quote=FALSE,sep=" ")
  write.table(tfam,tfam.name,row.names = FALSE,col.names=FALSE,quote=FALSE,sep=" ")
  return(tped)
}

#' Use a CERVUS file (or similar) to create add-on files for PRIMUS (e.g., age_file, sex_file)
#' @param gty A data.frame of genotypes in CERVUS/Plink format, with the first few columns including other info
#' @param filename A string that is the prefix for all output files
#' @param age Boolean, indicating whether you want to create an age file
#' @param sex Boolean, indicating whether you want to create a sex file
#' @param sex.col Number indicating which column contains the sex information
#' @param sexes A labeled list with Males, Females, and Offspring elements all designating strings found in ID names that identify individuals as males and females
#' @param off.string The string to use to identify offspring based on the ID. If age=TRUE, either this or age.col must be specified
#' @param age.col The column number containing age information
#' @param affection Boolean, indicating whether you want to create an affection status file
#' @param affection.col The column number containing affection status information
#' @param trait Boolean, indicating whether you want to create a binary trait information file
#' @param trait.col The column number containing trait information
#' @export
makePRIMUSextras<-function(gty,filename,out.dir="./",id.col=1,fid.col=NA,known.mom=NA,known.dad=NA,
                           age=TRUE,sex=TRUE,sex.col=NA,
                           sexes=c(Males="MAL",Females="FEM"),
                           off.string="OFF",age.col=NA,affection=FALSE,affection.col=NA,
                           trait=FALSE,trait.col=NA){
  if(is.na(fid.col)){
    fid<-0
    if(!is.na(known.mom)){ #then family names will come from moms
      fid[grep(sexes["Females"],gty[,known.mom])]<-paste("Fam",gsub(sexes["Females"],"",gty[grep(sexes["Females"],gty[,known.mom]),known.mom]),sep="")
      fid[grep(sexes["Females"],gty$ID)]<-paste("Fam",gsub(sexes["Females"],"",gty$ID[grep(sexes["Females"],gty$ID)]),sep="")
    }
    if(!is.na(known.dad)){ #then family names will come from moms
      fid[grep(sexes["Males"],gty[,known.dad])]<-paste("Fam",gsub(sexes["Males"],"",gty[grep(sexes["Males"],gty[,known.dad]),known.dad]),sep="")
      fid[grep(sexes["Males"],gty$ID)]<-paste("Fam",gsub(sexes["Males"],"",gty$ID[grep(sexes["Males"],gty$ID)]),sep="")
    }
  } else {
    fid<-gty[,fid.col]
  }
  # Do we make a sex file?
  if(isTRUE(sex)){ #then we make a sex file!
    if(is.na(sex.col)){
      sex.out<-as.character(gty[,id.col])
      sex.out[grep(sexes["Males"],sex.out)]<-1
      sex.out[grep(sexes["Females"],sex.out)]<-2
      sex.out[!sex.out %in% c(1,2)]<-0
    } else {
      if(length(grep("^1$",gty[,sex.col]))==0 | 
         length(grep("^2$",gty[,sex.col]))==0 | 
         length(grep("^0$",gty[,sex.col]))==0){ #make sure it's in the correct format
        sex.out<-as.character(gty[,sex.col])
        sex.out[grep(sexes["Males"],sex.out)]<-1
        sex.out[grep(sexes["Females"],sex.out)]<-2
        sex.out[!sex.out %in% c(1,2)]<-0
      } else{
        sex.out<-gty[,sex.col]
      }
    }
    sex.out<-data.frame(FID=fid,IID=gty[,id.col],SEX=sex.out)
    write.table(sex.out,paste(out.dir,filename,".sex",sep=""),sep=" ",
                col.names = FALSE,row.names=FALSE,quote=FALSE)
  } #end of crating sex file
  # Do we make an age file?
  if(isTRUE(age)){
    if(is.na(age.col)){ #then we use IDs to find offspring
      ages<-as.character(gty[,id.col])
      ages[grep(off.string,ages)]<-0
      ages[ages != 0]<-1
    } else{
      ages<-gty[,age.col]
    }
    age.out<-data.frame(FID=fid,IID=gty[,id.col],SEX=ages)
    write.table(age.out,paste(out.dir,filename,".age",sep=""),sep=" ",
                col.names = FALSE,row.names=FALSE,quote=FALSE)
  }
  # Do we make affection status file?
  if(isTRUE(affection)){
    if(class(affection.col)=="character"){ #sanity check
      affection.col<-grep(affection.col,colnames(gty))
    }
    write.table(data.frame(FID=fid,IID=gty[,id.col],AFFECTION_STATUS=gty[,affection.col]),
                paste(out.dir,filename,".affection",sep=""),sep=" ",
                col.names = FALSE,row.names=FALSE,quote=FALSE)
  }
  # Do we make a trait file?
  if(isTRUE(trait)){
    if(class(trait.col)=="character"){ #sanity check
      trait.col<-grep(trait.col,colnames(gty))
    }
    write.table(data.frame(FID=fid,IID=gty[,id.col],TRAIT=gty[,trait.col]),
                paste(out.dir,filename,".btrait",sep=""),sep=" ",
                col.names = FALSE,row.names=FALSE,quote=FALSE)
  }
}

#' Generate Clapper option files
#' @param out.name
#' @param File
#' @param refpopFile=NA
#' @param computeLike=1
#' @param ageFile=NA
#' @param errorRate=0.01
#' @param maxGen=5
#' @param maxSampleDepth=NA
#' @param condLD=1
#' @param back=0.04
#' @param startTemp=100
#' @param tempFact=1.01
#' @param iterPerTemp=40000
#' @param maxIter=30000000
#' @param conv=0.0001
#' @param poissonMean=NA
#' @param beta=30
#' @param numRun=3
#' @param numThreads=2
#' @export
getClapperOptions<-function(out.name,File,refpopFile=NA,computeLike=1,
                            ageFile=NA,errorRate=0.01,maxGen=5,
                            maxSampleDepth=NA,condLD=1,back=0.04,
                            startTemp=100,tempFact=1.01,iterPerTemp=40000,maxIter=as.character("30000000"),
                            conv=as.character("0.0001"),poissonMean=NA,beta=30,
                            numRun=3,numThreads=2){
  if(is.na(refpopFile)) refpopFile<-File
  if(is.na(ageFile)) ageFile<-"N/A"
  if(is.na(maxSampleDepth)) maxSampleDepth<-maxGen
  options<-data.frame(c(File,refpopFile,computeLike,ageFile,errorRate,
                        maxGen,maxSampleDepth,condLD,back,startTemp,
                        tempFact,as.character(iterPerTemp),as.character(maxIter),as.character(conv),
                        poissonMean,beta,numRun,numThreads),
                      c("#fileName","#refPopFileName","#computeLikelihood",
                        "#ageFileName","#errorRate","#maxGen","#maxSampleDepth",
                        "#conditiononld","#back","#startTemp","#tempFact",
                        "#iterPerTemp","#maxIter","#conv","#poissonMean","#beta",
                        "#numRun","#numThreads"))
  out.opt<-options[!is.na(options[,1]),]
  write.table(out.opt,out.name,sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
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
  gty$sex[gty$sex!=0]<-as.numeric(as.factor(gty$sex[gty$sex!=0]))
  gty$output<-1
  gty$ID<-as.numeric(as.factor(gty$ID))
  new.gty<-cbind(gty[,c("ID","age","sex","output")],gty[,first.allele:last.col])
  write.table(new.gty,filename,row.names = FALSE,col.names=FALSE,
              sep=",",quote=FALSE)
}


#' Convert a CERVUS file into input for FaMoz
#' @param gty A data.frame of genotypes in CERVUS format
#' @param filename The name of the file to export the FaMoz file to.
#' @param first.allele The first column with an allele (locus 1, alele 1). Default is 4.
#' @param miss The character specifying missing data. Default is 0. 
#' @param type Default is "snps". For microsatellites it should be "M" followed by the length of repeat-unit.
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
  if(type=="snps"){
    bps<-sample(c("A","C","G","T"),2)
    gty[gty==1]<-bps[1]
    gty[gty==2]<-bps[2]
  }else{
    nalleles<-max(gty[,first.allele:ncol(gty)])
    n<-as.numeric(gsub("M(\\d+)","\\1",type))
    bps<-sample(seq(130,400,n),nalleles,replace = FALSE)
    for(i in 1:nalleles){ gty[gty==i]<-bps[i] }
  }
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


#' Convert a CERVUS file into input for SNPPIT
#' @param gty A data.frame of genotypes in CERVUS format
#' @param filename The name of the file to export the SNPPIT file to.
#' @param first.allele The first column with an allele (locus 1, alele 1). Default is 4.
#' @param error.rates Either a static error rate or a data.frame with all of 
#' @param miss The character specifying missing data. Default is 0.  
#' @param offspring The string that specifies offspring
#' @param sexes If included, it should be a list of strings that specify the two sexes, and those get included in the output file (e.g., c("MAL","FEM"))
#' @param known.mom The column with known moms. Use NA to not include.
#' @param known.dad The column with known dads. Use NA to not include.
#' @export
cervus2snppit<-function(gty,filename,error.rates=0.01,first.allele=4,miss=0,
                        offspring="OFF",sexes=NA,known.mom=2,known.dad=NA){
  # starting info
  nloci<-data.frame(info=c(((ncol(gty)-(first.allele-1))/2),miss))
  rownames(nloci)<-c("NUMLOCI","MISSING_ALLELE")
  suppressWarnings(write.table(nloci,filename,quote=FALSE,col.names = FALSE,row.names = TRUE,sep=" "))
  if(!is.na(sexes)){
    suppressWarnings(write.table("POPCOLUMN_SEX",filename,quote=FALSE,col.names = FALSE,row.names = TRUE,sep=" ",append=TRUE))  
  }
  # error rates
  if(class(error.rates)=="data.frame"){
    suppressWarnings(write.table(error.rates,filename,quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t",append=TRUE))
  } else{
    er<-data.frame(locus.names=colnames(gty)[seq(first.allele,ncol(gty),2)],
                   rates=error.rates)
    suppressWarnings(write.table(er,filename,quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t",append=TRUE))
  }
  # The parents as the POP (I've only got one pop)
  parents<-gty[grep(offspring,gty$ID,invert = TRUE),]
  progeny<-gty[grep(offspring,gty$ID),]
  if(!is.na(sexes)){
    sex<-substr(gty$ID,start = 1,stop=3)
    sex[!sex %in% sexes]<-"?"
    parents[,2]<-sex
    parents<-parents[,-3]
  } else{
    parents<-parents[,-c(2:3)]
  }
  suppressWarnings(write.table("POP Parents",filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
  suppressWarnings(write.table(parents,filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
  # offspring
  if(!is.na(known.mom) & !is.na(known.dad)){
    all.moms<-unique(gty[,known.mom])[!is.na(unique(gty[,known.mom]))]
    all.dads<-unique(gty[,known.dad])[!is.na(unique(gty[,known.dad]))]
    counter<-1
    for(i in 1:length(all.moms)){
      this.mom<-progeny[progeny[,known.mom] %in% all.moms[i],]
      for(ii in 1:length(all.dads)){
        this.group<-progeny[progeny[,known.dad] %in% all.dads[ii],]
        if(!is.na(sexes)){
          output<-this.group[,-3]
          output[,2]<-"?"
        } else{
          output<-this.group[,-c(2:3)]
        }
        cnt.name<-paste("Offspring",counter,sep="")
        suppressWarnings(write.table(paste("OFFSPRING",cnt.name,"Parents"),filename,
                                     append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
        suppressWarnings(write.table(output,filename,append = TRUE,
                                     row.names=FALSE,col.names = FALSE,quote=FALSE))
        counter<-counter+1
      }
    }
    
  } else if(!is.na(known.mom) & is.na(known.dad)){
    all.moms<-unique(gty[,known.mom])[!is.na(unique(gty[,known.mom]))]
    counter<-1
    for(i in 1:length(all.moms)){
      this.mom<-progeny[progeny[,known.mom] %in% all.moms[i],]
      if(!is.na(sexes)){
        output<-this.mom[,-3]
        output[,2]<-"?"
      } else{
        output<-this.mom[,-c(2:3)]
      }
      cnt.name<-paste("Offspring",counter,sep="")
      suppressWarnings(write.table(paste("OFFSPRING",cnt.name,"Parents"),filename,
                                   append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
      suppressWarnings(write.table(output,filename,append = TRUE,
                                   row.names=FALSE,col.names = FALSE,quote=FALSE))
      counter<-counter+1
    }
    
  } else if(is.na(known.mom) & !is.na(known.dad)){
    all.dads<-unique(gty[,known.dad])[!is.na(unique(gty[,known.dad]))]
    counter<-1
    for(i in 1:length(all.dads)){
      this.dad<-progeny[progeny[,known.dad] %in% all.dads[i],]
      if(!is.na(sexes)){
        output<-this.dad[,-3]
        output[,2]<-"?"
      } else{
        output<-this.dad[,-c(2:3)]
      }
      cnt.name<-paste("Offspring",counter,sep="")
      suppressWarnings(write.table(paste("OFFSPRING",cnt.name,"Parents"),filename,
                                   append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
      suppressWarnings(write.table(output,filename,append = TRUE,
                                   row.names=FALSE,col.names = FALSE,quote=FALSE))
      counter<-counter+1
    }
  } else if(is.na(known.mom) & is.na(known.dad)){
    # then they all go into one
    if(!is.na(sexes)){
      output<-progeny[,-3]
      output[,2]<-"?"
    } else{
      output<-progeny[,-c(2:3)]
    }
    suppressWarnings(write.table("OFFSPRING All Parents",filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
    suppressWarnings(write.table(output,filename,append = TRUE,row.names=FALSE,col.names = FALSE,quote=FALSE))
  }
  
  
}
