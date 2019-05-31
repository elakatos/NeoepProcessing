library(purrr)
library(ggpubr)
library(scales)
library(reshape2)


# Neoantigen summary functions --------------------------------------------

colReds = c('#fee0d2','#fc9272','#de2d26')
colBlues = c('#deebf7','#9ecae1','#3182bd')

# Read in .avinput file that was the base of Annovar annotation
readAvinput <- function(sampleFile){
  avinput <- read.table(sampleFile, header=F, stringsAsFactors = F)
  regionNames <- getRegionNames(ncol(avinput)-18, T)
  names(avinput)[c(1,2,3,4,5, 18:ncol(avinput))] <- c('Chrom', 'Start', 'End', 'RefAll', 'AltAll', regionNames)
  return(avinput)
}

# Read in multi-region exonic annotated file
readExonicFile <- function(sampleFile){
  exonic <- read.table(sampleFile, header=F, sep='\t', stringsAsFactors = F)
  regionNames <- getRegionNames(ncol(exonic)-21, T)
  names(exonic)[c(1:8, 21:ncol(exonic))] <- c('LineID', 'MutType', 'MutInfo', 'Chrom', 'Start', 'End',
                                              'RefAll', 'AltAll', regionNames)
  return(exonic)
}

# Read in exonic annotated file produced from (single region) TCGA
readExonicFileTCGA <- function(sampleFile){
  exonic <- read.table(sampleFile, header=F, sep='\t', stringsAsFactors = F)
  names(exonic)[c(1:8)] <- c('LineID', 'MutType', 'MutInfo', 'Chrom', 'Start', 'End',
                             'RefAll', 'AltAll')
  return(exonic)
}

# Generate region names of multi-region antigen prediction file
getRegionNames <- function(n, normal=F)
{
  if (n<=0){
    return(NULL)
  }
  regionNames <- c(paste0(rep('Region_',n), 0:(n-1)))
  if (normal) {regionNames <- c('Normal', regionNames)}
  return(regionNames)
}

# Take a single sample out of antigen prediction table
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

# Compute VAF values from exonic files with information NR:NV
computeVaf <- function(readData, colInd, nrInd, nvInd){
  readInfo <- strsplit(readData[,colInd], ':')
  vafs <- as.numeric(map(readInfo, nvInd))/as.numeric(map(readInfo, nrInd))
  return(vafs)
}

# Compute VAF values from exonic files with allelic depth information: ref,alt
computeVafAD <- function(readData, colInd){
  readVec <- readData[,colInd]
  readVec <- readVec[readVec!='.']
  readInfo <- strsplit(readVec, ':')
  AD <- map(readInfo, 2)
  ADvec <- sapply(AD, function(i) strsplit(i, ','))
  vafs <- as.numeric(map(ADvec, 2))/(as.numeric(map(ADvec, 1))+as.numeric(map(ADvec,2)))
  return(vafs)
}

# Count number of antigenic mutations per regions and categories
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


# Get the number of mutations that were evaluated for antigenicity
getTotalMutFromFasta <- function(dir, sample){
  #returns only mutations that were inputted into netMHC
  fFile <- paste0(dir, '/fastaFiles/',sample,'.tmp.10.fasta')
  fData <- scan(file=fFile, what='string')
  return(sum(grepl('>', fData)))
}

# Compute table of antigenic mutations
recalculateSummaryTable <- function(epTable, summaryTable, mutations=T){
  summaryTableMut <- summaryTable
  for (sample in row.names(summaryTableMut)){
    eps <- subsetEpTable(epTable, sample, uniqueMutations = mutations)
    epsTot <- nrow(subsetEpTable(epTable, sample, F))
    summaryTableMut <- getStats(eps, sample, summaryTableMut)
    summaryTableMut[sample, 'Neoep'] <- epsTot
  }
  return(summaryTableMut)
}


processSummaryOfSampleSet <- function(dir, epTable, prefix){
  summaryTable <- read.table(paste0(dir,'/Neopred_results/',prefix,'.neoantigens.summarytable.txt'), header=T, row.names=1)
  summaryTable <- summaryTable[unique(epTable$Sample),]
  
  #Compute statistics from deduplicated information (1 entry per mutation in epTable)
  summaryTableMut <- recalculateSummaryTable(epTable, summaryTable)
  summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMutFromFasta(dir, x))
  
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


# HLA processing functions ------------------------------------------------

# Create a label for filtering depending on output of CN measurement
labelHLAAI <- function(line){
  if (is.na(line$HLA_type1copyNum_withBAFBin_upper)){return(NaN)}
  lab <- 'none'
  if (line$PVal_unique<0.01){lab <- 'AI'}
  if ((line$HLA_type1copyNum_withBAFBin_lower>1.5) & (line$HLA_type2copyNum_withBAFBin_lower > 1.5 ))
  { lab <- 'gain'}
  if ( (lab=='AI') & (line$HLA_type1copyNum_withBAFBin_upper< 0.7) & (line$HLA_type1copyNum_withBAFBin < 0.5 ))
  { lab <- 'LOH'}
  if ((lab=='AI') & (line$HLA_type2copyNum_withBAFBin < 0.5) & (line$HLA_type2copyNum_withBAFBin_upper < 0.7 ))
  { lab <- 'LOH'}
  return(lab)
}

# Create a label for filtering depending on the reliability
labelHLArel <- function(line, mincov){
  if (is.na(line$HLA_type1copyNum_withBAFBin_upper)){return(F)}
  lab <- T
  confint1 <- line$HLA_type1copyNum_withBAFBin_upper - line$HLA_type1copyNum_withBAFBin_lower
  confint2 <- line$HLA_type2copyNum_withBAFBin_upper - line$HLA_type2copyNum_withBAFBin_lower
  if ( ((confint1>1.5) & (line$HLA_type1copyNum_withBAFBin<1.5) & (line$HLA_type1copyNum_withBAFBin>0)) |  ((confint2>1.5) & (line$HLA_type2copyNum_withBAFBin<1.5) & (line$HLA_type2copyNum_withBAFBin>0)) )
  { lab <- F}
  if ( (confint1>6) | (confint2>6))
  {lab <- F}
  if ( ((abs(line$HLA_type1copyNum_withBAFBin)<0.4) & (line$HLA_type1copyNum_withBAFBin_upper>0.8)) |   ((abs(line$HLA_type2copyNum_withBAFBin)<0.4) & (line$HLA_type2copyNum_withBAFBin_upper>0.8)))
  {lab <- F}
  if ( line$numMisMatchSitesCov < mincov )
  {lab <- F}
  return(lab)
}

# Evaluate and filter results of LOHHLA output
analyseLohhla <- function(lohhla.master, clin.df){
  lohhla.signif <- lohhla.master[!is.na(lohhla.master$PVal_unique),]
  lohhla.signif <- lohhla.signif[lohhla.signif$PVal_unique<0.01,]
  
  lohhla.patients <- data.frame(row.names = unique(lohhla.master$region))
  lohhla.patients$ID <- sapply(row.names(lohhla.patients), function(x) substr(x,1,12))
  lohhla.patients$AI <- sapply(row.names(lohhla.patients), function(x) x %in% lohhla.signif$region)
  lohhla.patients$CN <- sapply(row.names(lohhla.patients),
                               function(x) lohhla.master[lohhla.master$region==x, 'Label'][1]!='NaN')
  lohhla.patients$LOH <- sapply(row.names(lohhla.patients),
                                function(x) ('LOH' %in% lohhla.master[lohhla.master$region==x, 'Label']))
  lohhla.patients$LOSS <- sapply(row.names(lohhla.patients),
                                 function(x) ('loss' %in% lohhla.master[lohhla.master$region==x, 'Label']))
  lohhla.patients$NORM <- sapply(row.names(lohhla.patients),
                                 function(x) sum(lohhla.master[lohhla.master$region==x, 'Label']!='none')==0)
  lohhla.patients$HIGH <- sapply(row.names(lohhla.patients),
                                 function(x) ('gain' %in% lohhla.master[lohhla.master$region==x, 'Label']))
  
  lohhla.patients$MSI <- clin.df[match(lohhla.patients$ID, clin.df$Patient), 'MSI']=='MSI-H'
  lohhla.patients$HYP <- clin.df[match(lohhla.patients$ID, clin.df$Patient), 'Hypermut']==1
  return(list(sig=lohhla.signif, pat=lohhla.patients))
}
