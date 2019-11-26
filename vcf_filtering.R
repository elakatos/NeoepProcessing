normInd = 10
setwd('~/RNAseq/Neoepitopes/')

fillNormal <- function(x, normInd, altnormInd){
  if (!is.na(altnormInd)) {
    for (ind in altnormInd){
      x[x[,normInd]=='.', normInd] <- x[x[,normInd]=='.', ind]
    }
    x <- x[, -(altnormInd)] #Delete alternative normal columns
  }
  
  return(x)
}

filterSomatic <- function(vcf, normInd){
  normalAD <- sapply(vcf[,normInd], function(x) unlist(strsplit(x, ':'))[2]) #get AD field of normal
  splitAD <- strsplit(normalAD, ',') #Split it to ref and alt
  normAD = as.numeric(map(splitAD,1))
  altAD = as.numeric(map(splitAD,2))
  
  vcf.somatic <- vcf[ ((normAD+altAD)>9 & (normAD+altAD)<30 & altAD<1 )| ((normAD+altAD)>29 & altAD<2 ) | ((normAD+altAD)>99  & altAD<3 ),  ]
  return(vcf.somatic)
}

filterSomatic2 <- function(vcf, normInd){
  normalNV <- sapply(vcf[,normInd], function(x) as.numeric(tail(unlist(strsplit(x, ':')),1)))
  normalNR <- sapply(vcf[,normInd], function(x) as.numeric(tail(unlist(strsplit(x, ':')),2)[1]))
  
  vcf.somatic <- vcf[(normalNV<1 & normalNR>9 & normalNR<30) | (normalNV<2 & normalNR>29) | (normalNV<3 & normalNR>99), ]
  return(vcf.somatic)
}

#fileList <- read.table('IBD_vcf_list.txt', header=T,stringsAsFactors = F, sep='\t')



for (i in 2:nrow(meta.df)){
  print(meta.df$FileName[i])
  print(meta.df$Patient[i])
vcf <- read.table(paste0('~/CRCdata/TCGA_UCEC/VCF_raw/',meta.df$FileName[i]), header=F, sep='\t', stringsAsFactors = F)
names(vcf)<- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', 'NORMAL','TUMOR')


#vcf.pass <- vcf[vcf$FILTER %in% c('PASS', 'alleleBias', 'HapScore', 'SC', 'badReads', 'SC;alleleBias', 'HapScore;alleleBias', 'HapScore;SC'),] #Only with PASS flag
vcf.pass <- subset(vcf, FILTER == 'PASS')
vcf.singlealt <- vcf.pass[!grepl(',',vcf.pass$ALT),] #Only with single ALT allele
#vcf.filled <- fillNormal(vcf.singlealt, normInd, altnormInd) #Unite normal information in one column
vcf.somatic <- filterSomatic(vcf.singlealt, 'NORMAL')
print(dim(vcf.somatic))

cat(c('##fileformat=VCFv4.2\n#'), file=paste0('~/CRCdata/TCGA_UCEC/VCF/',meta.df$Patient[i],'.vcf'))
write.table(vcf.somatic, file=paste0('~/CRCdata/TCGA_UCEC/VCF/',meta.df$Patient[i],'.vcf'),
            sep='\t', row.names=F, quote=F, append=T)
}
