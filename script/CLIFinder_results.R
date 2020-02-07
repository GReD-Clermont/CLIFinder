#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
missing_args = TRUE

if (length(args)==0 && exists("out1") && exists("out2") && exists("out3") && exists("nfastq") && exists("line_only") && exists("refseq"))
{
  missing_args = FALSE
} else if(length(args)==6)
{
  missing_args = FALSE
  out1 = args[1]
  out2 = args[2]
  out3 = args[3]
  nfastq = as.numeric(args[4])
  line_only = args[5]
  refseq = args[6]
}

if(missing_args)
{
  stop("6 arguments must be supplied or the variables must be set before calling the script.", call.=FALSE)
} else {
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(plyr))

  chim<-read.delim(out1)
  
  chim<-chim[order(chim[,nfastq+7],decreasing=F),]
  chim<-chim[order(chim[,2],decreasing=F),]
  chr<-sub("chr","",as.character(chim[,1]))
  suppressWarnings(chim<-chim[order(as.numeric(chr)),])
  
  grchim <- GRanges(seqnames = chim[,1],
                    IRanges(start = chim[,2],
                            end = chim[,3]),strand=chim[,4])
  
  grchim$ID<-paste("Id_",1:length(chim[,1]),sep="")
  
  mcols(grchim)<-cbind(mcols(grchim),chim[,5:(nfastq+9)])
  
  grfusR<- union(grchim,grchim)
  
  suppressWarnings(position<-as.data.frame(findOverlaps(grfusR,grchim)))
  
  grfusR$dup<- as.vector(table(position[,1]))
  
  suppressWarnings(position2<-as.data.frame(findOverlaps(grfusR[grfusR$dup>1],grchim)))
  
  grfusR2<-grfusR
  strand(grfusR2)<- "+"
  gr3<-union(grfusR2,grfusR2)
  suppressWarnings(position3<- as.data.frame(findOverlaps(gr3,grfusR2)))
  gr3$dup<-as.vector(table(position3[,1]))
  
  suppressWarnings(position3<-as.data.frame(findOverlaps(gr3[gr3$dup>1],grfusR2)))
  grfusR$info<-"no"
  grfusR$info [position3[,2]]<-"overlap sens opposé"
  
  grfusR$ID<-"Id"
  grfusR$ID[position[!duplicated(position[,1]),1]]<-grchim$ID[position[!duplicated(position[,1]),2]]
  
  if(nrow(position2)!=0)
  {
    result <- aggregate(position2[,2] ~ position2[,1], data = position2, paste, collapse = "_")
    grfusR$ID[grfusR$dup>1]<-paste("ID",result[,2],sep="_")
  }
  
  mcols(grfusR)<-cbind(mcols(grfusR), mcols(grchim[position[!duplicated(position[,1]),2]]))
  
  min<-ddply(as.data.frame(grchim), .(seqnames,end,strand), function(x)x[x$Chimera.start==min(x$Chimera.start),])
  min<-ddply(as.data.frame(min), .(seqnames,start,strand), function(x)x[x$Chimera.start==min(x$Chimera.start),])
  max<-ddply(as.data.frame(grchim), .(seqnames,end,strand), function(x)x[x$Chimera.end==max(x$Chimera.end),])
  max<-ddply(as.data.frame(max), .(seqnames,start,strand), function(x)x[x$Chimera.end==max(x$Chimera.end),])
  
  grfusR<-as.data.frame(grfusR)
  grfusR<-grfusR[order(grfusR[,1],grfusR[,2],grfusR[,3],grfusR[,4],decreasing=F),]
  
  grfusR$Chimera.start<- min$Chimera.start
  grfusR$Chimera.end<-max$Chimera.end
  
  datax<-as.data.frame(grfusR)
  colnames(datax)[1:3]<-colnames(chim)[1:3]
  colnames(datax)[5]<-colnames(chim)[4]
  
  grchim2 <- GRanges(seqnames = datax[,nfastq+11],
                     IRanges(start = datax[,nfastq+12],
                             end = datax[,nfastq+13]),strand=datax[,nfastq+14])
  
  mcols(grchim2)<-datax[,-c(4,nfastq+11:nfastq+14)]
  
  grfus<- union(grchim2,grchim2)
  
  suppressWarnings(position<-as.data.frame(findOverlaps(grfus,grchim2)))
  
  grfus$dup<- as.vector(table(position[,1]))
  
  suppressWarnings(position2<-as.data.frame(findOverlaps(grfus[grfus$dup>1],grchim2)))
  
  grfus$ID_final<-"Id"
  grfus$ID_final[position[!duplicated(position[,1]),1]]<-grchim2$ID[position[!duplicated(position[,1]),2]]
  
  if(nrow(position2)!=0)
  {
    result <- aggregate(position2[,2] ~ position2[,1], data = position2, paste, collapse = "_")
    grfus$ID_final[grfus$dup>1]<-paste("Id",result[,2],sep="_")
  }
  
  mcols(grfus)<-cbind(mcols(grfus), mcols(grchim2[position[!duplicated(position[,1]),2]]))
  
  for (i in 0:nfastq)
  {
    mcols(grfus)[grfus$dup>1,11+i] <- mcols(grfus)[grfus$dup>1,11+i] +  mcols(grchim2)[position[duplicated(position[,1]),2],9+i]
  }
  
  grfus2<-grfus
  strand(grfus2)<-"+"
  gr3<-union(grfus2,grfus2)
  
  suppressWarnings(position3<- as.data.frame(findOverlaps(gr3,grfus2)))
  gr3$dup<-as.vector(table(position3[,1]))
  
  suppressWarnings(position3<-as.data.frame(findOverlaps(gr3[gr3$dup>1],grfus2)))
  grfus$info [position3[,2]]<-"overlap sens opposé"
  
  min<-ddply(as.data.frame(grchim2), .(seqnames,end,strand), function(x)x[x$L1.start==min(x$L1.start),])
  min<-ddply(data.frame(min), .(seqnames,start,strand), function(x)x[x$L1.start==min(x$L1.start),])
  max<-ddply(as.data.frame(grchim2), .(seqnames,end,strand), function(x)x[x$L1.end==max(x$L1.end),])
  max<-ddply(data.frame(max), .(seqnames,start,strand), function(x)x[x$L1.end==max(x$L1.end),])
  
  grfus1<-as.data.frame(grfus)
  grfus1<-grfus1[order(grfus1[,1],grfus1[,2],grfus1[,3],grfus1[,4],decreasing=F),]
  
  grfus1$L1.start<- min$L1.start
  grfus1$L1.end<-max$L1.end
  
  dataf<-as.data.frame(grfus1)
  
  result<-( data.frame("Chimera.Chr" = dataf$L1.chromosome,
		       "Chimera.Start" = apply(data.frame(dataf$start,
							  dataf$end,
							  dataf$L1.start,
							  dataf$L1.end),1,min),
		       "Chimera.End" = apply(data.frame(dataf$start,dataf$end,dataf$L1.start,dataf$L1.end), 1, max),
		       "Chimera.Strand" = dataf$L1.strand,
		       "L1.Chr" = dataf$L1.chromosome,
		       "L1.Start" = dataf$L1.start,
		       "L1.End" = dataf$L1.end,
		       "L1.Strand" = dataf$L1.strand,
		       "Unique.Chr" = dataf$seqnames,
		       "Unique.Start" = dataf$start,
		       "Unique.End" = dataf$end,
		       "Unique.Strand" = dataf$strand,
		       "ID_final" = dataf$ID_final,
		       "info" = dataf$info,
		       dataf[,16:(nfastq+16)]))
  
  result<-result[order(result[,2],decreasing=F),]
  chr<-sub("chr","",as.character(result[,1]))
  suppressWarnings(result<-result[order(as.numeric(chr)),])
  options(scipen=10)
  write.table(result,out2,sep="\t",row.names = F,quote = F)
  grchim <- GRanges(seqnames = result$L1.Chr,
                    IRanges(start = result$L1.Start,
                    end = result$L1.End),strand=result$L1.Strand)
  mcols(grchim)<-result
  
  Rep<-read.delim(line_only,skip=1)
  grLINE <- GRanges(seqnames = Rep$genoName,
                    IRanges(start = Rep$genoStart,
                            end = Rep$genoEnd),
                    repStrand = as.character(Rep$strand),
                    repName = as.character(Rep$repName))
  
  Gene<-read.delim(refseq)
  grGene <- GRanges(seqnames = Gene$chrom,
                    IRanges(start = Gene$txStart,
                            end = Gene$txEnd),
                    geneStrand = as.character(Gene$strand),
                    geneName = as.character(Gene$name2))
  
  suppressWarnings(position<-as.data.frame(findOverlaps(grchim,grLINE)))
  suppressWarnings(position2<-as.data.frame(findOverlaps(grchim,grGene)))
  
  grchim$GeneName<-"no_gene"
  grchim$GeneName[position2[,1]]<- grGene$geneName[position2[,2]]
  
  grchim$GeneStrand<-"*"
  grchim$GeneStrand[position2[,1]]<- grLINE$repStrand[position2[,2]]
  
  grchim$repName<-"no"
  grchim$repName[position[,1]]<- grLINE$repName[position[,2]]
  
  grchim$repStart<-0
  grchim$repStart[position[,1]]<-start(grLINE[position[,2]])
  
  grchim$repEnd<-0
  grchim$repEnd[position[,1]]<-end(grLINE[position[,2]])
  
  grchim$repWidth<-0
  grchim$repWidth[position[,1]]<-width(grLINE[position[,2]])
  
  grchim$repStrand<-"*"
  grchim$repStrand[position[,1]]<- grLINE$repStrand[position[,2]]
  
  dup<-position[duplicated(position[,1]),1]
  if(length(dup != 0))
  {
    for (i in 1:length(dup))
    {
      grchim$repName[dup[i]] <-paste(grLINE$repName[position[position[,1]==dup[i],2]],collapse="/")
      grchim$repStart[dup[i]] <-paste(start(grLINE[position[position[,1]==dup[i],2]]),collapse="/")
      grchim$repEnd[dup[i]] <-paste(end(grLINE[position[position[,1]==dup[i],2]]),collapse="/")
      grchim$repWidth[dup[i]] <-paste(width(grLINE[position[position[,1]==dup[i],2]]),collapse="/")
      grchim$repStrand[dup[i]] <-paste(grLINE$repStrand[position[position[,1]==dup[i],2]],collapse="/")
    }
  }
  
  final_result<-as.data.frame(grchim)
  options(scipen=10)
  write.table(final_result[,-c(1:5)],out3,sep="\t",row.names = F,quote = F)
  print("Executed R script")
}

