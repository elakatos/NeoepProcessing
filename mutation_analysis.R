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
random.data.real <- subset(random.data, !((nchar(random.data$peptide)==9) &  (random.data$peptide_pos) %in% c(1,11)) )
random.data.real <- subset(random.data.real, Novelty==1)
random.data.real$Sample <- random.data.real$PatIndex
random.data.filtered <- subset(random.data.real, BindLevel!='N')

# random.dataBA <- read.table('random_proteome_all_BA.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)
# random.data.realBA <- subset(random.dataBA, !((nchar(random.dataBA$peptide)==9) &  (random.dataBA$peptide_pos) %in% c(1,11)) )
# random.data.realBA <- subset(random.data.realBA, Novelty==1)
# random.data.realBA$Sample <- random.data.realBA$PatIndex
# random.data.filteredBA <- subset(random.data.realBA, BindLevel!='N')
# 
# names(random.data.filtered)[11] <- 'ID'; names(random.data.filteredBA)[11] <- 'ID'
# random.data.filtered <- getSharedEps(random.data.filtered, random.data.filteredBA)
# names(random.data.filtered)[11] <- 'Identity'

random.summary <- data.frame(matrix(vector(), nrow=length(unique(random.data$PatIndex))))
row.names(random.summary) <- unique(random.data$PatIndex)
random.summary$AllPeptides <- sapply(row.names(random.summary), function(x) sum(random.data.real$PatIndex==as.numeric(x)))
random.summary$AllMuts <- sapply(row.names(random.summary),
                                 function(x) length(unique(random.data.real[random.data.real$PatIndex==as.numeric(x),]$Identity)))

random.summary$Epitopes <- sapply(row.names(random.summary), function(x) sum(random.data.filtered$PatIndex==as.numeric(x)))
random.summary$EpMuts <- sapply(row.names(random.summary),
                                 function(x) length(unique(random.data.filtered[random.data.filtered$PatIndex==as.numeric(x),]$Identity)))

random.summary$SB <- sapply(row.names(random.summary),
                                  function(x) sum((random.data.filtered$PatIndex==as.numeric(x)) * (random.data.filtered$BindLevel=='SB')))


lines(density(random.summary$EpMuts/random.summary$AllMuts), col='red')
barplot(random.summary$Epitopes/random.summary$EpMuts)
barplot(random.summary$SB/random.summary$Epitopes)

hist(random.data.filtered[random.data.filtered$PatIndex==14,]$Rank, breaks=20)



# Compare to normal peptide -----------------------------------------------

filterAllWT <- function(WT, Mut){
  nonWT<- row.names(subset(WT, BindLevel=='N'))
  return(nonWT)
}


filterByWTBinding <- function(dir, epTable){
  WTTable <- read.table(paste0(dir, '/Neopred_results/WT.neoantigens.txt'), header=T,
                        sep = '\t',stringsAsFactors = F, fill=T)
  WTTable <- subset(WTTable, !((nchar(WTTable$peptide)==9) &  (WTTable$peptide_pos) %in% c(1,11)) )
  
  WTTable[nchar(WTTable$peptide)==9,]$peptide_pos <- WTTable[nchar(WTTable$peptide)==9,]$peptide_pos-1
  
  WTTable$Ident <- sapply(1:nrow(WTTable), function(i) paste0(WTTable[i,'Identity'], WTTable[i,'peptide_pos'], WTTable[i, 'hla'], nchar(WTTable[i, 'peptide']), substr(WTTable[i, 'PatIndex'], 1, nchar(WTTable[i,'PatIndex'])-7)))
  row.names(WTTable) <- WTTable$Ident
  epTable$Ident <- sapply(1:nrow(epTable), function(i) paste0(epTable[i,'ID'], epTable[i,'pos'], epTable[i, 'hla'], nchar(epTable[i, 'peptide']),epTable[i, 'Sample']))
  row.names(epTable) <- epTable$Ident
  epTable$MutID <- sapply(1:nrow(epTable), function(i) paste0(epTable[i, 'LineID'], epTable[i, 'Sample']))
  
  WTTable.matched <- WTTable[epTable$Ident,]
  
  nonWT <- filterAllWT(WTTable.matched, epTable)
  epTablemut <- epTable[nonWT,]
  
  p1 = ggplot(epTable, aes(x=Score)) + geom_density(aes(fill='All eps'), alpha=0.4) + geom_density(data=epTablemut, aes(fill='Filtered eps'), alpha=0.4) + theme_minimal()
  
  
  epWTvsMut <- data.frame((table(epTablemut$Sample)/table(epTable$Sample)))
  epDedup <- epTable[!duplicated(epTable$MutID),]
  epMutDedup <- epTablemut[!duplicated(epTablemut$MutID),]
  epWTvsMut <- rbind(epWTvsMut, data.frame((table(epMutDedup$Sample)/table(epDedup$Sample))))
  epWTvsMut$Type <- c(rep('Epitopes', length(unique(epTablemut$Sample))), rep('Mutations', length(unique(epTablemut$Sample))))
    
  p2 = ggplot(epWTvsMut, aes(x=Var1, y=Freq, fill=Type)) + geom_bar(stat='identity', color = 'black',position=position_dodge()) +
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Reds")
  
  print(p1); print(p2);
  
  return(epTablemut)
}

# HLAs --------------------------------------------------------------------

cat(hlasA, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt')
cat('\n', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat(hlasB, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat('\n', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat(hlasC, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
