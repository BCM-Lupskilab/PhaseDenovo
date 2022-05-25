library(data.table)
library(tidyverse)
library(GenomicRanges)
library(regioneR)

PhaseLongRead <- function(path_to_gVCF,min_snp=20,ref_genome,param){
  require(data.table)
  require(dplyr)
  require(VariantAnnotation)
  require(pbmcapply) # for parallel run
  print("**Step1:read vcf file**")
  vcf<- VariantAnnotation::readVcf(file = path_to_gVCF,genome = ref_genome,param = param)
  vcf.gr <- vcf@rowRanges
  rowRanges <- data.table::as.data.table(vcf@rowRanges)
  GT <- VariantAnnotation::geno(vcf)$GT
  PS <- vcf@assays@data$PS # phased block id
  GQ <- vcf@assays@data$GQ # genotype quality
  PR_ID=colnames(GT)[1]
  P1_ID=colnames(GT)[2]
  P2_ID=colnames(GT)[3]
  G1=c('0/0',"0|0")
  G2=c('1/1',"1|1")
  G3=c('0/1',"0|1","1|0")
  GT <- as.data.table(GT)
  setnames(GT,colnames(GT),c("index","P1","P2"))
  GT.anno <- GT %>% mutate(InhFrom=ifelse(index%in%G3&P1%in%G1&P2%in%c(G2,G3),P2_ID,
                                          ifelse(index%in%G3&P1%in%c(G2,G3)&P2%in%G1,P1_ID,"Notphased")))
  PS.GQ <- as.data.table(cbind(PS,GQ))
  setnames(PS.GQ,colnames(PS.GQ),c("index_PS","P1_PS","P2_PS","index_GQ","P1_GQ","P2_GQ"))
  rowRanges_psgq <- cbind(rowRanges,PS.GQ,GT.anno)
  
  ## select the long read block with at least 20 SNP phased
  print("**Step2:assign haplotype based on majority vote**")
  block.ID <-rowRanges_psgq%>%dplyr::filter(index_GQ>20,P1_GQ>20,P2_GQ>20,!is.na(index_PS))%>%
    group_by(index_PS)%>%
    summarise(haplen=max(start)-min(start),snp.count=n()) %>%
    filter(snp.count>min_snp)%>%dplyr::select(index_PS)%>%unlist()
  block_table <- rowRanges_psgq%>%
    dplyr::filter(index_GQ>20,P1_GQ>20,P2_GQ>20,index_PS%in%block.ID,InhFrom!="Notphased")%>%
    dplyr::select(index,index_PS,InhFrom)
  block.ID.p <- pbmcapply::pbmclapply(1:length(block.ID), function(i)
  {
    print(i)
    conf.matrix <- block_table%>%filter(index_PS==block.ID[i])%>%dplyr::select(index,InhFrom)
    conf.matrix <- table(conf.matrix)
    conf.matrix[1,1]
    p <- tryCatch({(conf.matrix[1,1]+conf.matrix[2,2])/sum(conf.matrix)},error=function(e){return(NA)})
    if(is.na(p)){return(NULL)}
    ## use major-vote strategy to determine the haplotype
    if(which.max(c(p,1-p))==1){hap=paste0(rownames(conf.matrix)[1],":",colnames(conf.matrix)[1],";",
                                          rownames(conf.matrix)[2],":",colnames(conf.matrix)[2])}
    if(which.max(c(p,1-p))==2){hap=paste0(rownames(conf.matrix)[2],":",colnames(conf.matrix)[1],";",
                                          rownames(conf.matrix)[1],":",colnames(conf.matrix)[2])}
    return(data.table(PS=block.ID[i],TPratio=max(p,1-p),n.phasedsnp=sum(conf.matrix),haplotype=hap))
  },mc.cores = 1)
  #block.ID.p <- data.table::rbindlist(block.ID.p)
  return(block.ID.p)
}
PhaseDenovo <- function(path_to_gVCF,ref_genome,denovo.gr){
  window.size =  1e6
  param <- GenomicRanges::flank(denovo.gr,both = T,width = window.size)
  phased_block <- PhaseLongRead(path_to_gVCF,min_snp=20,ref_genome,param=param)%>%
    data.table::rbindlist()
  vcf.obj <- VariantAnnotation::readVcf(file = path_to_gVCF,genome = ref_genome,denovo.gr)
  denovo.gr$ps <- vcf.obj@assays@data$PS
  denovo.gr$gt <- geno(vcf.obj)$GT
  hap <- phased_block%>%
    filter(PS==denovo.gr$ps[1,1])%>%
    dplyr::select(haplotype)%>%unlist()%>%base::strsplit(.,";")%>%
    unlist()
  denovo.gr$hap <- hap[pmatch(denovo.gr$gt[1,1],hap)]
  return(denovo.gr)
}


