library(vioplot)
library(ggpubr)

#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Generate summary table for MUTATIONS ------------------------------------

getSharedEps <- function(epTableEL, epTableBA){
  epTable <- data.frame(matrix(vector(), ncol=ncol(epTableEL)))
  names(epTable) <- names(epTableEL)
  for (sample in unique(epTableEL$Sample)){
    epsub <- subset(epTableEL, Sample==sample)
    idsEL <- epsub$ID
    idsBA <- subset(epTableBA, Sample==sample)$ID
    epsub <- subset(epsub, ID %in% intersect(idsEL, idsBA))
    epTable <- rbind(epTable, epsub)
  }
  return(epTable)
}

getStats <- function(eps, sample, table){
  tumorColumns = grep('Region*', names(eps))
  eps = eps[!duplicated(eps$LineID),]
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

getMutationTable <- function(epTable, summaryTable){
summaryTableMut <- summaryTable
for (sample in row.names(summaryTableMut)){
  print(sample)
  eps <- subsetEpTable(epTable, sample)
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
  mutS <- row.names(mutTable)[(rowSums(mutTable)>1) & (rowSums(mutTable)<ncol(mutTable))]
  ratioTable[sample, 'Clonal_All'] <- length(mutC)
  ratioTable[sample, 'Private_All'] <- length(mutP)
  ratioTable[sample, 'Shared_All'] <- length(mutS)
  ratioTable[sample, 'Clonal_Ep'] <- sum(mutC %in% eps$LineID)
  ratioTable[sample, 'Private_Ep'] <- sum(mutP %in% eps$LineID)
  ratioTable[sample, 'Shared_Ep'] <- sum(mutS %in% eps$LineID)
}
return(ratioTable)
}


# Get statistics for Polyp and Set ---------------------------------------------------

setwd('~/CRCdata')

mutRatiosBatch = list()
mutRatioTable = data.frame(matrix(vector(), ncol=6))
names(mutRatioTable) <- c('Clonal_All', 'Clonal_Ep','Private_All','Private_Ep','Shared_All', 'Shared_Ep')

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Set')
#dirList <- c('IBD')

for (dir in dirList){
pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/Output.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'ID', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epTable <- epTable[epTable$Novelty!=0,]
epTable <- epTable[epTable$Sample!='Set.10.snv',] #to disregard the replication of Set10

# epTableBA <- read.table(paste0(dir, '/Neopred_results/Output_BA.neoantigens.txt'), header=F,
#                       sep = '\t',stringsAsFactors = F, fill=T)
# names(epTableBA) <- c('Sample', getRegionNames(ncol(epTableBA)-24), 'LineID', 'Chrom', 'Start',
#                     'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
#                     'Gl', 'Ip', 'Il', 'Icore', 'ID', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
# epTableBA <- epTableBA[epTableBA$Novelty!=0,]
# epTableBA <- epTableBA[epTableBA$Sample!='Set.10.snv',]
# 
# epTable <- getSharedEps(epTable, epTableBA)

summaryTable <- read.table(paste0(dir,'/Neopred_results/Output.neoantigens.summarytable.txt'), header=T, row.names=1)
summaryTable <- summaryTable[unique(epTable$Sample),]

barcolors = c('firebrick', 'darkorange3', 'goldenrod1')
barplot(t(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')]/summaryTable$Total)), col=barcolors, legend=c('Clonal', 'Shared','Subclonal'), las=2)
barplot(t(as.matrix(summaryTable[,c('Total_WB', 'Total_SB')]/summaryTable$Total)), col=barcolors, legend=c('Weak binders', 'Strong binders'), las=2)

summaryTableMut <- getMutationTable(epTable, summaryTable)

barplot(t(as.matrix(summaryTableMut[,c('Clonal', 'Shared','Subclonal')]/summaryTableMut$Total)), col=barcolors, legend=c('Clonal', 'Shared','Subclonal'), las=2)
barplot(t(as.matrix(summaryTableMut[,c('Total_WB', 'Total_SB')]/summaryTableMut$Total)), col=barcolors, legend=c('Weak binders', 'Strong binders'), las=2)

barplot(summaryTable$Total/summaryTableMut$Total, col=barcolors[1], main='Average neoepitopes per neoep mutation')
summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMutFromFasta(dir, x))
barplot(summaryTableMut$Total/summaryTableMut$Total_MUT, col=barcolors[1], main='Neoepitope mutations/ all mutations')

mutRatioTable <- getMutationRatios(dir, epTable, mutRatioTable)
mutRatiosBatch[dir] <- list(summaryTableMut$Total/summaryTableMut$Total_MUT)
dev.off()
}


mutRatiosBatch['Random_proteome'] <- list(random.summary$EpMuts/random.summary$AllMuts)


# Analyse Polyp, Set and Random -------------------------------------------

ggqqplot(mutRatiosBatch[[2]])


plot(log(mutRatioTable$Clonal_All), log(mutRatioTable$Clonal_Ep), pch=19, col=c(rep(2,5),rep(3,11)))
segments(2,2, 5.5, 0.85*5.5, col='grey50')

var.test(mutRatiosBatch[[1]], mutRatiosBatch[[2]])
t.test(mutRatiosBatch[[1]], mutRatiosBatch[[2]])

pdf('Neoepitope_ratio.pdf', height=5, width=8)
plot.ecdf(mutRatiosBatch[[2]], col='darkred', xlim=c(0.5, 0.95))
plot.ecdf(mutRatiosBatch[[1]], col='steelblue4',add=T)
#plot.ecdf(mutRatiosBatch[[3]], col='darkgreen', add=T)
plot.ecdf(mutRatiosBatch[[3]], col='grey75', add=T)

plot(density(mutRatiosBatch[[2]]), col='darkred', xlim=c(0.5, 0.95), ylim=c(0, 8))
lines(density(mutRatiosBatch[[1]]), col='steelblue4')
lines(density(mutRatiosBatch[[3]]), col='grey75')
dev.off()
t.test(mutRatiosBatch[[1]], mutRatiosBatch[[2]])


plot(density((mutRatioTable$Clonal_Ep/mutRatioTable$Clonal_All)[1:5]), ylim=c(0, 7.5))
lines(density((mutRatioTable$Private_Ep/mutRatioTable$Private_All)[1:5]), col='darkred')



# TCGA sample -------------------------------------------------------------

dir = 'TCGA_COAD'
pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/TCGA_COAD.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'ID', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]

summaryTable <- read.table(paste0(dir,'/Neopred_results/TCGA_COAD.neoantigens.summarytable.txt'), header=T, row.names=1)

barcolors = c('firebrick', 'darkorange3', 'goldenrod1')
#barplot(t(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')]/summaryTable$Total)), col=barcolors, legend=c('Clonal', 'Shared','Subclonal'), las=2)
barplot(t(as.matrix(summaryTable[,c('Total_WB', 'Total_SB')]/summaryTable$Total)), col=barcolors, legend=c('Weak binders', 'Strong binders'), las=2)

summaryTableMut <- getMutationTable(epTable, summaryTable)
summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMut(dir, x))
dev.off()


# General mutations stats -------------------------------------------------

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Set')
DNDAs = list()

for (dir in dirList){
summaryTable <- read.table(paste0(dir,'/Neopred_results/CRCmseq.neoantigens.summarytable.txt'), header=T, row.names=1)
dnda = vector()
for (sample in row.names(summaryTable)){
  sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
  exonic <- readExonicFile(sampleFileEx)
  dnda = c(dnda, (sum(exonic$MutType=='nonsynonymous SNV')/sum(exonic$MutType=='synonymous SNV')))
}
DNDAs[dir] = list(dnda)
}

