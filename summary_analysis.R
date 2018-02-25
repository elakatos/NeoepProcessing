library(vioplot)

#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Generate summary table for MUTATIONS ------------------------------------

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
  clonaleps = eps[rowSums(eps[,tumorColumns])==length(tumorColumns),]
  nonclonaleps = eps[rowSums(eps[,tumorColumns])!=length(tumorColumns),]
  table[sample, 'Clonal'] = nrow(clonaleps)
  table[sample, 'Clonal_WB'] = sum(clonaleps$BindLevel=='WB')
  table[sample, 'Clonal_SB'] = sum(clonaleps$BindLevel=='SB')
  sharedeps = nonclonaleps[rowSums(nonclonaleps[,tumorColumns])>1,]
  subclonaleps = nonclonaleps[rowSums(nonclonaleps[,tumorColumns])==1,]
  table[sample, 'Shared'] = nrow(sharedeps)
  table[sample, 'Shared_WB'] = sum(sharedeps$BindLevel=='WB')
  table[sample, 'Shared_SB'] = sum(sharedeps$BindLevel=='SB')
  table[sample, 'Subclonal'] = nrow(subclonaleps)
  table[sample, 'Subclonal_WB'] = sum(subclonaleps$BindLevel=='WB')
  table[sample, 'Subclonal_SB'] = sum(subclonaleps$BindLevel=='SB')
  
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
  return(nrow(exonic))
}

getMutationTable <- function(epTable, summaryTable){
summaryTableMut <- summaryTable
for (sample in row.names(summaryTableMut)){
  eps <- subsetEpTable(epTable, sample)
  summaryTableMut <- getStats(eps, sample, summaryTableMut)
}

return(summaryTableMut)
}


# #What's the ratio between produced neoepitopes and mutations
# barplot(summaryTable$Total/summaryTableMut$Total)
# 
# #What's the distribution of neoepitopes per mutation
# summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMut(dir, x))
# barplot(summaryTableMut$Total/summaryTableMut$Total_MUT)
# barplot(summaryTableMut$Total_SB/summaryTableMut$Total_MUT)


getMutationRatios <- function(dir, summaryTableMut){
tumorColumns <- grep('Total_Region*', names(summaryTableMut))
epMut <- vector()
allMut <- vector()
for (sample in row.names(summaryTableMut)){
  tc <- tumorColumns[summaryTableMut[sample, tumorColumns]!=0]
  for (region in (names(summaryTableMut)[tc])){
    epMut <- c(epMut, summaryTableMut[sample, region])
    allMut <- c(allMut, getTotalMut(dir, sample, substr(region, 7, nchar(region))))
  }
}
plot(epMut, allMut, pch=19, xlab='Neoepitope mutations', ylab='All mutations')
plot(epMut/allMut, pch=19, ylab='Neoepitope/all mutations')
return(epMut/allMut)
}


# Compare Polyp and Set ---------------------------------------------------

setwd('~/CRCdata')

mutRatios = list()
mutRatiosBatch = list()

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Set')

for (dir in dirList){
pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/CRCmseq.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'ID', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]

summaryTable <- read.table(paste0(dir,'/Neopred_results/CRCmseq.neoantigens.summarytable.txt'), header=T, row.names=1)

barcolors = c('firebrick', 'darkorange3', 'goldenrod1')
barplot(t(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')]/summaryTable$Total)), col=barcolors, legend=c('Clonal', 'Shared','Subclonal'), las=2)
barplot(t(as.matrix(summaryTable[,c('Total_WB', 'Total_SB')]/summaryTable$Total)), col=barcolors, legend=c('Weak binders', 'Strong binders'), las=2)

summaryTableMut <- getMutationTable(epTable, summaryTable)

barplot(summaryTable$Total/summaryTableMut$Total, col=barcolors[1], main='Average neoepitopes per neoep mutation')
summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMut(dir, x))

mutRatios[dir] <- list(getMutationRatios(dir, summaryTableMut))
mutRatiosBatch[dir] <- list(summaryTableMut$Total/summaryTableMut$Total_MUT)
dev.off()
}


pdf('CRCmseq_comparison_summary.pdf', height=5, width=8)
qqplot(mutRatios[[2]], mutRatios[[1]], pch=19, xlab='Neoepitope/all mutations in Carcinoma', ylab='Neoepitope/all mutations in Adenoma', main='QQplot')
plot.ecdf(mutRatios[[2]], col='firebrick3', xlim=c(0.35, 0.65), ylab='CDF',
          xlab='Neoepitope/all mutations', main=paste0('KS test p-value: ', ks.test(mutRatios[[1]], mutRatios[[2]])$p.value),
          legend=c('Carcinoma', 'Adenoma'))
plot.ecdf(mutRatios[[1]], col='skyblue3',add=T)
plot.ecdf(mutRatiosBatch[[2]], col='darkred', add=T)
plot.ecdf(mutRatiosBatch[[1]], col='steelblue4',add=T)

vioplot(mutRatios[[2]], mutRatios[[1]], col='wheat3', names = c('Carcinoma', 'Adenoma'))
title('Neoeptiope/all mutation ratio')
dev.off()





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

