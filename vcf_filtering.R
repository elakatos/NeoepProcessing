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
  normalAD <- sapply(vcf[,normInd], function(x) tail(unlist(strsplit(x, ':')),1)) #get AD field of normal
  splitAD <- strsplit(normalAD, ',') #Split it to ref and alt
  
  vcf.somatic <- vcf[((map(splitAD,1)>99) | (map(splitAD,2)==0)),]
  return(vcf.somatic)
}

filterSomatic2 <- function(vcf, normInd){
  normalNV <- sapply(vcf[,normInd], function(x) as.numeric(tail(unlist(strsplit(x, ':')),1)))
  normalNR <- sapply(vcf[,normInd], function(x) as.numeric(tail(unlist(strsplit(x, ':')),2)[1]))
  
  vcf.somatic <- vcf[(normalNV<1 & normalNR>9 & normalNR<30) | (normalNV<2 & normalNR>29) | (normalNV<3 & normalNR>99), ]
  return(vcf.somatic)
}

fileList <- read.table('IBD_vcf_list.txt', header=T,stringsAsFactors = F, sep='\t')



for (fileName in fileList){
  print(fileName)
vcf <- read.table(fileName, header=F, sep='\t', stringsAsFactors = F)
names(vcf)[1:9] <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')


vcf.pass <- vcf[vcf$FILTER %in% c('PASS', 'alleleBias', 'HapScore', 'SC', 'badReads', 'SC;alleleBias', 'HapScore;alleleBias', 'HapScore;SC'),] #Only with PASS flag
vcf.singlealt <- vcf.pass[!grepl(',',vcf.pass$ALT),] #Only with single ALT allele
#vcf.filled <- fillNormal(vcf.singlealt, normInd, altnormInd) #Unite normal information in one column
vcf.somatic <- filterSomatic2(vcf.singlealt, 'NORMAL')
print(dim(vcf.somatic))

cat(c('#'), file=fileName)
write.table(vcf.somatic, file=fileName, sep='\t', row.names=F, quote=F, append=T)
}
