library(purrr)

readAvinput <- function(sampleFile){
  avinput <- read.table(sampleFile, header=F, stringsAsFactors = F)
  regionNames <- getRegionNames(ncol(avinput)-18, T)
  names(avinput)[c(1,2,3,4,5, 18:ncol(avinput))] <- c('Chrom', 'Start', 'End', 'RefAll', 'AltAll', regionNames)
  return(avinput)
}

readExonicFile <- function(sampleFile){
  exonic <- read.table(sampleFile, header=F, sep='\t', stringsAsFactors = F)
  regionNames <- getRegionNames(ncol(exonic)-21, T)
  names(exonic)[c(1:8, 21:ncol(exonic))] <- c('LineID', 'MutType', 'MutInfo', 'Chrom', 'Start', 'End',
                                                      'RefAll', 'AltAll', regionNames)
  return(exonic)
}

readExonicFileTCGA <- function(sampleFile){
  exonic <- read.table(sampleFile, header=F, sep='\t', stringsAsFactors = F)
  names(exonic)[c(1:8)] <- c('LineID', 'MutType', 'MutInfo', 'Chrom', 'Start', 'End',
                                              'RefAll', 'AltAll')
  return(exonic)
}

getRegionNames <- function(n, normal=F)
{
  if (n<=0){
    return(NULL)
  }
  regionNames <- c(paste0(rep('Region_',n), 0:(n-1)))
  if (normal) {regionNames <- c('Normal', regionNames)}
  return(regionNames)
}

subsetEpTable <- function(epTable, sample, uniqueMutations=F){
  epitopes <- epTable[epTable$Sample==sample,]
  epitopes <- epitopes[, colSums(epitopes!=(-1))>0]
  cat('Total number of epitopes in sample: ', nrow(epitopes))
  if (uniqueMutations){
    epitopes <- epitopes[!duplicated(epitopes$LineID),]
    cat('\n Number of unique mutations giving rise to epitopes: ', nrow(epitopes))
  }
  return(epitopes)
}

computeVaf <- function(readData, colInd){
  readInfo <- strsplit(readData[,colInd], ':')
  vafs <- as.numeric(map(readInfo, 2))/as.numeric(map(readInfo, 1))
  return(vafs)
}

computeVafAD <- function(readData, colInd){
  readVec <- readData[,colInd]
  readVec <- readVec[readVec!='.']
  readInfo <- strsplit(readVec, ':')
  AD <- map(readInfo, 11)
  ADvec <- sapply(AD, function(i) strsplit(i, ','))
  vafs <- as.numeric(map(ADvec, 2))/(as.numeric(map(ADvec, 1))+as.numeric(map(ADvec,2)))
  return(vafs)
}


epTable <- read.table(paste0(dir, '/Neopred_results/Output.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'ID', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]
#barplot(table(epNon$Sample)/table(epTable$Sample)*100, las=2)
epTableStrong <- epTable[epTable$BindLevel=='SB',]

sample = 'Oxford_IBD1.mutectCalls..somatic'
#sampleFile <- paste0(dir, '/avready/',sample,'.avinput')
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
#avinput <- readAvinput(sampleFile)
exonic <- readExonicFile(sampleFileEx)
eps <- subsetEpTable(epTable, sample, unique=T)

tumorColumns <- grep('Region*', names(exonic))
isEpMutation <- (exonic$LineID %in% eps$LineID)
uB <- 0.3
lB <- 0.05
pdf(paste0(dir,':',sample,'_s.pdf'),height=5,width=8)
par(mfrow=c(1,2))
for (i in tumorColumns){
  allVafs <- computeVafAD(exonic, i)
  epVafs <- computeVafAD(exonic[isEpMutation,], i)
  
  allVafsF <- allVafs[(allVafs>lB) & (allVafs< uB)]
  epVafsF <- epVafs[(epVafs>lB) & (epVafs < uB)]
  
  print(length(allVafsF))
  qqplot(allVafsF, epVafsF, pch=19, xlab='All mutations', ylab='Neoepitope mutations', main='QQplot')
  plot.ecdf(allVafsF,col='black', ylab='CDF', xlab='VAF')
  plot.ecdf(epVafsF,col='grey50', add=T)
  print(ks.test(allVafsF, epVafsF))
}
dev.off()


# Epitope distribution ----------------------------------------------------

tumorColumns = grep('Region*', names(eps))
hist(eps$Rank, breaks=20)
hist(eps[rowSums(eps[, tumorColumns])==4,]$Rank, breaks=20  )

epRankClonal <- eps[rowSums(eps[, tumorColumns])==4,]$Rank
epRankNotClonal <- eps[rowSums(eps[, tumorColumns])<4,]$Rank

#setwd('~/RNAseq/Neoepitopes/')
random.data <- read.table('random_proteome_all.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)
names(random.data)[3] <- 'peptide'
random.data.real <- subset(random.data, !((nchar(random.data$peptide)==9) &  (random.data$peptide_pos) %in% c(1,11)) )
random.data.non <- subset(random.data.real, Novelty==0)
random.data.real <- subset(random.data.real, Novelty==1)

random.summary <- data.frame(matrix(vector(), nrow=length(unique(random.data$PatIndex))))
row.names(random.summary) <- unique(random.data$PatIndex)
random.summary$AllPeptides <- sapply(row.names(random.summary), function(x) sum(random.data.real$PatIndex==as.numeric(x)))
random.summary$AllMuts <- sapply(row.names(random.summary),
                                 function(x) length(unique(random.data.real[random.data.real$PatIndex==as.numeric(x),]$Identity)))
random.data.filtered <- subset(random.data.real, BindLevel!='N')

random.summary$Epitopes <- sapply(row.names(random.summary), function(x) sum(random.data.filtered$PatIndex==as.numeric(x)))
random.summary$EpMuts <- sapply(row.names(random.summary),
                                 function(x) length(unique(random.data.filtered[random.data.filtered$PatIndex==as.numeric(x),]$Identity)))

random.summary$SB <- sapply(row.names(random.summary),
                                  function(x) sum((random.data.filtered$PatIndex==as.numeric(x)) * (random.data.filtered$BindLevel=='SB')))


barplot(random.summary$EpMuts/random.summary$AllMuts)
barplot(random.summary$Epitopes/random.summary$EpMuts)
barplot(random.summary$SB/random.summary$Epitopes)

hist(random.data.filtered[random.data.filtered$PatIndex==14,]$Rank, breaks=20)




# HLAs --------------------------------------------------------------------

cat(hlasA, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt')
cat('\n', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat(hlasB, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat('\n', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat(hlasC, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
