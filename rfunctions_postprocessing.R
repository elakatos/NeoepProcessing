library(purrr)
library(ggpubr)


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


getStats <- function(eps, sample, table){
  tumorColumns = grep('Region*', names(eps))
  table[sample, 'Total'] = nrow(eps)
  table[sample, 'Total_WB'] = sum(eps$BindLevel=='WB')
  table[sample, 'Total_SB'] = sum(eps$BindLevel=='SB')
  for (region in (names(eps)[tumorColumns])){
    table[sample, paste0('Total_',region)] = sum(eps[,region]==1)
    table[sample, paste0('Total_WB_',region)] = sum((eps[,region]==1) * (eps$BindLevel=='WB'))
    table[sample, paste0('Total_SB_',region)] = sum((eps[,region]==1) * (eps$BindLevel=='SB'))
  }
  clonaleps = eps[rowSums(eps[,tumorColumns,drop=F])==length(tumorColumns),]
  nonclonaleps = eps[rowSums(eps[,tumorColumns,drop=F])!=length(tumorColumns),]
  table[sample, 'Clonal'] = nrow(clonaleps)
  table[sample, 'Clonal_WB'] = sum(clonaleps$BindLevel=='WB')
  table[sample, 'Clonal_SB'] = sum(clonaleps$BindLevel=='SB')
  sharedeps = nonclonaleps[rowSums(nonclonaleps[,tumorColumns,drop=F])>1,]
  subclonaleps = nonclonaleps[rowSums(nonclonaleps[,tumorColumns,drop=F])==1,]
  table[sample, 'Shared'] = nrow(sharedeps)
  table[sample, 'Shared_WB'] = sum(sharedeps$BindLevel=='WB')
  table[sample, 'Shared_SB'] = sum(sharedeps$BindLevel=='SB')
  table[sample, 'Subclonal'] = nrow(subclonaleps)
  table[sample, 'Subclonal_WB'] = sum(subclonaleps$BindLevel=='WB')
  table[sample, 'Subclonal_SB'] = sum(subclonaleps$BindLevel=='SB')
  
  return(table)
}

getStatsTotal <- function(dir, sample, table){
  exonic <- getClonalityTotal(dir, sample)
  table[sample, 'Total'] = nrow(exonic)
  for (region in names(exonic)){
    table[sample, paste0('Total_',region)] = sum(exonic[,region])
  }
  clonal = exonic[rowSums(exonic)==ncol(exonic),]
  nonclonal = exonic[rowSums(exonic)!=ncol(exonic),]
  table[sample, 'Clonal'] = nrow(clonal)
  shared = nonclonal[rowSums(nonclonal)>1,]
  subclonal = nonclonal[rowSums(nonclonal)==1,]
  table[sample, 'Shared'] = nrow(shared)
  table[sample, 'Subclonal'] = nrow(subclonal)
  return(table)
}

getTotalMut <- function(dir, sample, region=''){
  sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
  if (length(grep('TCGA', dir))!=0){
    exonic <- readExonicFileTCGA(sampleFileEx)
  }
  else{
    exonic <- readExonicFile(sampleFileEx)}
  if (region!=''){
    regionMut <- as.numeric(map(strsplit(exonic[,region], ':'),2))>0
    exonic <- exonic[regionMut,]
  }
  #return(nrow(exonic))
  return(sum(exonic$MutType=='nonsynonymous SNV')) #get only nonsynonymous exonic mutations
}

getTotalMutFromFasta <- function(dir, sample){
  #returns only mutations that were inputted into netMHC
  fFile <- paste0(dir, '/fastaFiles/',sample,'.tmp.10.fasta')
  fData <- scan(file=fFile, what='string')
  return(length(fData)/2)
}

recalculateSummaryTable <- function(epTable, summaryTable, mutations=T){
  summaryTableMut <- summaryTable
  for (sample in row.names(summaryTableMut)){
    print(sample)
    eps <- subsetEpTable(epTable, sample, uniqueMutations = mutations)
    summaryTableMut <- getStats(eps, sample, summaryTableMut)
  }
  
  return(summaryTableMut)
}


getMutationTableTotal <- function(dir, sampleNames){
  summaryTableTotal <- data.frame(matrix(vector(), nrow=length(sampleNames)))
  row.names(summaryTableTotal) <- sampleNames
  for (sample in sampleNames){
    summaryTableTotal <- getStatsTotal(dir, sample, summaryTableTotal)
  }
  return(summaryTableTotal)
}


getClonalityTotal <- function(dir, sample){
  sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
  exonic <- readExonicFile(sampleFileEx)
  exonic <- exonic[,grep('Region*',names(exonic))]
  exonic[] <- lapply(1:(ncol(exonic)), function(i) as.numeric(map(strsplit(exonic[,i], ':'),2)))
  exonic <- exonic>0
  return(exonic)
}


getClonalityEp <- function(epTable, sample){
  eps <- subsetEpTable(epTable, sample, uniqueMutations = T)
  tumorcol <- grep('Region*', names(eps))
  eps <- eps[,tumorcol]
  return(eps)
}

getMutPresenceTable <- function(dir, sample, regions, gtInfo){
  fFile <- paste0(dir, '/fastaFiles/',sample,'.tmp.10.fasta')
  fData <- scan(file=fFile, what='string')
  fData <- fData[seq(1,length(fData),by=2)]
  mutLine <- map(strsplit(fData, ';'),1)
  mutIDs <- sapply(mutLine, function(i) substr(i, 2,nchar(i)))
  
  mutCloneTable <- data.frame(matrix(vector(), nrow=length(mutIDs), ncol=length(regions)))
  row.names(mutCloneTable) <- mutIDs
  names(mutCloneTable) <- regions
  
  sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
  exData <- readExonicFile(sampleFileEx)
  row.names(exData) <- exData$LineID
  exData <- exData[,regions,drop=F]
  for (line in mutIDs){
    if (gtInfo=='NV'){
      mutCloneTable[line,] <- sapply(exData[line,], function(i) 1*(as.numeric(unlist(strsplit(i, ':'))[2])>0))
    }
    if (gtInfo=='AD'){
      mutCloneTable[line,] <- sapply(exData[line,], function(i) 1*(as.numeric(strsplit(strsplit(i, ':')[[1]][length(strsplit(i, ':')[[1]])],',')[[1]][2])>0) )
    }
  }
  mutCloneTable[is.na(mutCloneTable)] <- 0
  return(mutCloneTable)
}


getMutationRatios <- function(dir, epTable, ratioTable, gtInfo='NV'){
  
  for (sample in unique(epTable$Sample)){
    eps <- subsetEpTable(epTable, sample, uniqueMutations = T)
    regions <- grep('Region*', names(eps), value = T)
    mutTable <- getMutPresenceTable(dir, sample, regions, gtInfo)
    mutC <- row.names(mutTable)[rowSums(mutTable)==ncol(mutTable)]
    mutP <- row.names(mutTable)[rowSums(mutTable)==1]
    mutS <- row.names(mutTable)[(rowSums(mutTable)>1) & (rowSums(mutTable)<(ncol(mutTable)-1))]
    ratioTable[sample, 'Clonal_All'] <- length(mutC)
    ratioTable[sample, 'Private_All'] <- length(mutP)
    ratioTable[sample, 'Shared_All'] <- length(mutS)
    ratioTable[sample, 'Clonal_Ep'] <- sum(mutC %in% eps$LineID)
    ratioTable[sample, 'Private_Ep'] <- sum(mutP %in% eps$LineID)
    ratioTable[sample, 'Shared_Ep'] <- sum(mutS %in% eps$LineID)
  }
  return(ratioTable)
}

filterByBAPrediction <- function(dir, epTableEL){
  epTableBA <- read.table(paste0(dir, '/Neopred_results/Output_BA.neoantigens.txt'), header=F,
                          sep = '\t',stringsAsFactors = F, fill=T)
  names(epTableBA) <- c('Sample', getRegionNames(ncol(epTableBA)-24), 'LineID', 'Chrom', 'Start',
                        'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                        'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
  epTableBA <- epTableBA[epTableBA$Novelty!=0,]
  epTableBA <- subset(epTableBA, Sample %in% unique(epTableEL$Sample))
  
  epTable <- data.frame(matrix(vector(), ncol=ncol(epTableEL)))
  names(epTable) <- names(epTableEL)
  for (sample in unique(epTableEL$Sample)){
    epsub <- subset(epTableEL, Sample==sample)
    idsEL <- epsub$Identity
    idsBA <- subset(epTableBA, Sample==sample)$Identity
    epsub <- subset(epsub, Identity %in% intersect(idsEL, idsBA))
    epTable <- rbind(epTable, epsub)
  }
  return(epTable)
}


filterAllWT <- function(WT, Mut){
  nonWT<- row.names(subset(WT, BindLevel=='N'))
  return(nonWT)
}


filterByWTBinding <- function(dir, epTable, randomsample = F){
  WTname <- paste0(dir, '/Neopred_results/WT.neoantigens.txt') }
  WTTable <- read.table(WTname, header=T,
                        sep = '\t',stringsAsFactors = F, fill=T)
  WTTable <- subset(WTTable, !((nchar(WTTable$peptide)==9) &  (WTTable$peptide_pos) %in% c(1,11)) )
  
    WTTable[nchar(WTTable$peptide)==9,]$peptide_pos <- WTTable[nchar(WTTable$peptide)==9,]$peptide_pos-1
    WTTable$Ident <- sapply(1:nrow(WTTable), function(i) paste0(WTTable[i,'Identity'], WTTable[i,'peptide_pos'], WTTable[i, 'hla'], nchar(WTTable[i, 'peptide']), substr(WTTable[i, 'PatIndex'], 1, nchar(WTTable[i,'PatIndex'])-7)))
  
  epTable$Ident <- sapply(1:nrow(epTable), function(i) paste0(epTable[i,'Identity'], epTable[i,'pos'], epTable[i, 'hla'], nchar(epTable[i, 'peptide']),epTable[i, 'Sample']))
  epTable$MutID <- sapply(1:nrow(epTable), function(i) paste0(epTable[i, 'LineID'], epTable[i, 'Sample']))
  
  row.names(WTTable) <- WTTable$Ident
  row.names(epTable) <- epTable$Ident
  
  WTTable.matched <- WTTable[epTable$Ident,]
  
  nonWT <- filterAllWT(WTTable.matched, epTable)
  epTablemut <- epTable[nonWT,]
  
  p1 = ggplot(epTable, aes(x=Score)) + geom_density(aes(fill='All eps'), alpha=0.4) + geom_density(data=epTablemut, aes(fill='Filtered eps'), alpha=0.4)+
    theme_minimal() + labs(fill = 'Epitopes')
  
  epWTvsMut <- data.frame((table(epTablemut$Sample)/table(epTable$Sample)))
  epDedup <- epTable[!duplicated(epTable$MutID),]
  epMutDedup <- epTablemut[!duplicated(epTablemut$MutID),]
  epWTvsMut <- rbind(epWTvsMut, data.frame((table(epMutDedup$Sample)/table(epDedup$Sample))))
  epWTvsMut$Type <- c(rep('Epitopes', length(unique(epTablemut$Sample))), rep('Mutations', length(unique(epTablemut$Sample))))
  
  p2 = ggplot(epWTvsMut, aes(x=Var1, y=Freq, fill=Type)) + geom_bar(stat='identity',position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Reds") +
    ggtitle("Percentage of epitopes retained after filtering based on corresponding WT binding") + labs(x="Tumour")
  
  print(p1); print(p2);
  
  return(epTablemut)
}


