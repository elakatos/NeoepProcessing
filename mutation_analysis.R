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

getRegionNames <- function(n, normal=F)
{
  regionNames <- c(paste0(rep('Region',n), 1:n))
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

getPeptideFasta <- function(x, outFileName){
  y <- sapply(1:nrow(x), function(z) cat('>',x[z,'LineID'],'-',x[z,'peptide'],'\n',
                                    x[z,'peptide'],'\n', file=outFileName, append=T, sep=''))
}


dir <- '~/CRCdata/CRCmseq_Set'
epTable <- read.table(paste0(dir, '/Neopred_results/CRCmseq.neoantigens.txt'), header=F, stringsAsFactors = F)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-22), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'ID', 'Score', 'Rank', 'Cand', 'BindLevel')

sample = 'Set.01.snv'
sampleFile <- paste0(dir, '/avready/',sample,'.avinput')
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
avinput <- readAvinput(sampleFile)
exonic <- readExonicFile(sampleFileEx)
eps <- subsetEpTable(epTable, sample)

tumorColumns <- grep('Region*', names(exonic))
isEpMutation <- (exonic$LineID %in% eps$LineID)
par(mfrow=c(1,2))
for (i in tumorColumns){
  allVafs <- computeVaf(exonic, i)
  epVafs <- allVafs[isEpMutation]
  hist(allVafs[allVafs>0], xlab='VAF', ylab='All mutations', main=names(exonic)[i])
  hist(epVafs[epVafs>0], xlab='VAF', ylab='Neoepitope mutations', main=names(exonic)[i])
  print(ks.test(allVafs[allVafs>0], epVafs[epVafs>0]))
}
dev.off()

getPeptideFasta(eps, paste0(dir,'/tmp/',sample,'.all_eps.fasta'))

# Epitope distribution ----------------------------------------------------

hist(-log(eps$Rank), breaks=50  )
