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
  
  vcf.somatic <- vcf[(floor(normalNV*0.01) >= normalNR), ]
  return(vcf.somatic)
}

fileList <- read.table('IBD_vcf_list.txt', header=T,stringsAsFactors = F, sep='\t')


for (i in 1:(nrow(fileList))){
altnormInd <- as.numeric(unlist(strsplit(fileList[i,,drop=F]$altnormInd,',')))
fileName <- fileList[i,,drop=F]$Sample
print(fileName)

vcf <- read.table(fileName, header=T, sep='\t', stringsAsFactors = F)

vcf.pass <- vcf[vcf$FILTER=='PASS',] #Only with PASS flag
vcf.singlealt <- vcf.pass[!grepl(',',vcf.pass$ALT),] #Only with single ALT allele
vcf.filled <- fillNormal(vcf.singlealt, normInd, altnormInd) #Unite normal information in one column
vcf.somatic <- filterSomatic2(vcf.singlealt, normInd)

write.table(vcf.somatic, file=paste0(fileName,'.somatic'), sep='\t', row.names=F, quote=F)
}
