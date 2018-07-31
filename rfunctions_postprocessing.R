library(purrr)
library(ggpubr)
library(scales)
library(reshape2)

colReds = c('#fee0d2','#fc9272','#de2d26')
colBlues = c('#deebf7','#9ecae1','#3182bd')


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


getTotalMutFromFasta <- function(dir, sample){
  #returns only mutations that were inputted into netMHC
  fFile <- paste0(dir, '/fastaFiles/',sample,'.tmp.10.fasta')
  fData <- scan(file=fFile, what='string')
  return(length(fData)/2)
}

recalculateSummaryTable <- function(epTable, summaryTable, mutations=T){
  summaryTableMut <- summaryTable
  for (sample in row.names(summaryTableMut)){
    eps <- subsetEpTable(epTable, sample, uniqueMutations = mutations)
    summaryTableMut <- getStats(eps, sample, summaryTableMut)
  }
  summaryTableMut$Neoep <- summaryTable[match(row.names(summaryTableMut), row.names(summaryTable)),'Total']
  
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

plotFilterStatistics <- function(epTable, epTableFiltered){
  
  p1 = ggplot(epTable, aes(x=Score)) + geom_density(aes(fill='All eps'), alpha=0.4) + geom_density(data=epTableFiltered, aes(fill='Filtered eps'), alpha=0.4)+
    labs(fill = 'Epitopes')
  
  #For each sample, compute the number of epitopes before and after filtering
  epFilteredvsAll <- data.frame((table(epTableFiltered$Sample)/table(epTable$Sample)))
  #Generate epTables that only contain one line per mutation
  epDedup <- epTable[!duplicated(epTable$MutID),]
  epDedupFiltered <- epTableFiltered[!duplicated(epTableFiltered$MutID),]
  epFilteredvsAll <- rbind(epFilteredvsAll, data.frame((table(epDedupFiltered$Sample)/table(epDedup$Sample))))
  epFilteredvsAll$Type <- c(rep('Epitopes', length(unique(epTableFiltered$Sample))), rep('Mutations', length(unique(epTableFiltered$Sample))))
  
  p2 = ggplot(epFilteredvsAll, aes(x=Var1, y=Freq, fill=Type)) + geom_bar(stat='identity',position=position_dodge(), colour='black') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Reds") +
    ggtitle("Percentage of neo-epitopes and neo-epitope associated mutations retained after filtering") + labs(x="Tumour")
  
  print(p1); print(p2);
}


filterAllWT <- function(WT, Mut){
  nonWT<- row.names(subset(WT, BindLevel=='N'))
  return(nonWT)
}

filterHigherWT <- function(WT, Mut){
  nonWT <- row.names(Mut[WT$Score < Mut$Score*0.3,])
  return(nonWT)
}


filterByWTBinding <- function(dir, epTable, filt){
  WTname <- paste0(dir, '/Neopred_results/WT.neoantigens.txt') 
  WTTable <- read.table(WTname, header=T,
                        sep = '\t',stringsAsFactors = F, fill=T)
  WTTable <- subset(WTTable, !((nchar(WTTable$peptide)==9) &  (WTTable$peptide_pos) %in% c(1,11)) )
  
  WTTable[nchar(WTTable$peptide)==9,]$peptide_pos <- WTTable[nchar(WTTable$peptide)==9,]$peptide_pos-1
  #Generate unique identifiers to match WT and mutated table entries
  WTTable$Ident <- sapply(1:nrow(WTTable), function(i) paste0(WTTable[i,'Identity'], WTTable[i,'peptide_pos'], WTTable[i, 'hla'], nchar(WTTable[i, 'peptide']), substr(WTTable[i, 'PatIndex'], 1, nchar(WTTable[i,'PatIndex'])-7)))
  epTable$Ident <- sapply(1:nrow(epTable), function(i) paste0(epTable[i,'Identity'], epTable[i,'pos'], epTable[i, 'hla'], nchar(epTable[i, 'peptide']),epTable[i, 'Sample']))
  epTable <- epTable[!duplicated(epTable$Ident),]
  row.names(WTTable) <- WTTable$Ident
  row.names(epTable) <- epTable$Ident
  sharedIdent <- intersect(WTTable$Ident, epTable$Ident) #To ensure that no NAs are introduced, in case no WT information available for some reason
  if (length(sharedIdent)<length(epTable$Ident)){print('Warning: some neoepitopes might miss matching wild-type information!\n')}
  WTTable.matched <- WTTable[sharedIdent,]; epTable <- epTable[sharedIdent,]
  
  #Generate unique mutation identifier for mutated table entries
  epTable$MutID <- sapply(1:nrow(epTable), function(i) paste0(epTable[i, 'LineID'], epTable[i, 'Sample']))

  #Filter
  if (filt %in% c('all', 'a')){
    nonWT <- filterAllWT(WTTable.matched, epTable)
  }
  else if (filt %in% c('higher', 'h')){
      nonWT <- filterHigherWT(WTTable.matched, epTable)
  }
  else {
    print('Filtering mode is not specified, all entries are retained.\n')
    nonWT <- sharedIdent
    }
  epTablemut <- epTable[nonWT,]
  
  plotFilterStatistics(epTable, epTablemut)
  
  return(epTablemut)
}

filterRandomByWTBinding <- function(random.data, filt){
  WTTable <- read.table('random_wt_proteome_all.txt', header=T,sep = '\t',stringsAsFactors = F, fill=T)
  WTTable <- subset(WTTable, !((nchar(WTTable$peptide)==9) &  (WTTable$peptide_pos) %in% c(1,11)) )
  #Match WT and mutated table entries (utilising the 1-1 correspondence)
  WTTable.matched <- WTTable[row.names(random.data),]
  
  #Generate unique mutation identifier for mutated table entries
  random.data$MutID <- sapply(1:nrow(random.data), function(i) paste0(random.data[i, 'Identity'],'**', random.data[i, 'PatIndex']))
  #Filter
  
  if (filt %in% c('all', 'a')){
    nonWT <- filterAllWT(WTTable.matched, random.data)
  }
  else if (filt %in% c('higher', 'h')){
    nonWT <- filterHigherWT(WTTable.matched, random.data)
  }
  else {
    print('Filtering mode is not specified, all entries are retained.\n')
    nonWT <- row.names(random.data)
  }
  random.data.filtered <- random.data[nonWT,]
  
  plotFilterStatistics(random.data, random.data.filtered)
  
  return(random.data.filtered)
  
}

processSummaryOfSampleSet <- function(dir, epTable, prefix){
  summaryTable <- read.table(paste0(dir,'/Neopred_results/',prefix,'.neoantigens.summarytable.txt'), header=T, row.names=1)
  summaryTable <- summaryTable[unique(epTable$Sample),]
  
  #Adjust summary table in case neo-epitopes have been filtered
  summaryTable <- recalculateSummaryTable(epTable, summaryTable, mutations = F)
  #Compute statistics from deduplicated information (1 entry per mutation in epTable)
  summaryTableMut <- recalculateSummaryTable(epTable, summaryTable)
  
  #Plot the percentage ratios of clonality information of neo-eps and neo-ep mutations
  clonality <- rbind(melt(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')])),
                     melt(as.matrix(summaryTableMut[,c('Clonal', 'Shared','Subclonal')])))
  clonality$Type <- c(rep('P',3*nrow(summaryTable)), rep('M', 3*nrow(summaryTable)))
  pc <- ggplot(data=clonality, aes(x=Type, fill=Var2)) + geom_bar(aes(y=value),stat="identity", position = position_fill(reverse=T), colour='black') +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values=rev(colReds)) + labs(x = "Tumour", y = "", fill="Clonality") + ggtitle("Clonality of neoepitopes") +
    facet_grid(. ~ Var1)
  #print(pc)
  
  #Plot the average number of epitopes (peptides) produced by mutations that lead to epitopes
  eps <- melt(as.matrix(summaryTable[, 'Total', drop=F]/summaryTableMut[,'Total',drop=F]))
  pe <- ggplot(data=eps, aes(x=Var1, y=value)) + geom_bar(stat='identity', fill=colBlues[3], color='black') +
    labs(x = "Tumour", y = "") + ggtitle("Average neo-epitopes per neo-ep mutation") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  #print(pe)
  
  #Plot the number of mutations giving rise to neo-epitopes over all missense (exonic, aa-changing) mutations
  summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMutFromFasta(dir, x))
  pem <- ggplot(data=summaryTableMut, aes(x=row.names(summaryTableMut), y=Total/Total_MUT)) + geom_bar(stat='identity', fill=colBlues[3], color='black') +
    scale_y_continuous(labels = percent_format()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x = "Tumour", y = "") + ggtitle("Percentage of missense mutations producing at least one neo-epitope")
  #print(pem)
  
  #Is there a connection between the percentage of clonal mutations and percentage of neo-ep mutations?
  pemc <- ggplot(data=summaryTableMut, aes(x=Clonal/Total, y=Total/Total_MUT)) + geom_point(size=3) +
    geom_smooth(method='lm', color=colBlues[3]) + labs(x="Percentage of mutations clonal", y="Percentage of mutation associated with neo-eps")
  #print(pemc)
  
  return(summaryTableMut)
}

processSummaryOfRandomSet <- function(random.data.real, random.data.filtered){
  #Create data frame to hold data
  random.summary <- data.frame(matrix(vector(), nrow=length(unique(random.data.real$PatIndex))))
  row.names(random.summary) <- unique(random.data.real$PatIndex)
  
  #Calculate all mutations and all novel peptides produced by them
  random.summary$AllPeptides <- sapply(row.names(random.summary), function(x) sum(random.data.real$PatIndex==as.numeric(x)))
  random.summary$AllMuts <- sapply(row.names(random.summary),
                                   function(x) length(unique(random.data.real[random.data.real$PatIndex==as.numeric(x),]$Identity)))
  #Calculate amount of neo-peptides that are neo-epitopes and all mutations producing associated with neo-epitopes
  random.summary$Epitopes <- sapply(row.names(random.summary), function(x) sum(random.data.filtered$PatIndex==as.numeric(x)))
  random.summary$EpMuts <- sapply(row.names(random.summary),
                                  function(x) length(unique(random.data.filtered[random.data.filtered$PatIndex==as.numeric(x),]$Identity)))
  random.summary$SB <- sapply(row.names(random.summary),
                              function(x) sum((random.data.filtered$PatIndex==as.numeric(x)) * (random.data.filtered$BindLevel=='SB')))
  
  return(random.summary)
}

