
# Helper functions --------------------------------------------------------

getGeneScores <- function(epTable, geneTable){
  altNameTable <- data.frame(matrix(vector(), ncol=1)); names(altNameTable) <- 'Alt.name'
  
  exprTable <- data.frame(matrix(0, ncol=length(unique(epTable$Sample)), nrow=length(unique(geneTable$Gene.name))))
  names(exprTable) <- unique(epTable$Sample)
  row.names(exprTable) <- unique(geneTable$Gene.name)
  
  for (i in row.names(epTable)){
    samp <- epTable[i, 'Sample']
    x <- epTable[i,]
    if (x$Gene.name %in% row.names(exprTable)){
      #We can use Gene.name as row-identifier
      gN <- x$Gene.name
    }else{
      if (x$Gene.name %in% row.names(altNameTable)){
        #we already came across this and use the recorded row-identifier
        gN <- altNameTable[x$Gene.name,'Alt.name']
      }else{
        #we identify the change by chromosome and position
        gene <- subset(geneTable, (Chrom == x$Chrom) & (Start < x$Start) & (End > x$Start) )
        if (nrow(gene)>=1){
          gN <- as.character(gene$Gene.name)[1]
          altNameTable[x$Gene.name,] <- gN
        }else{cat("No match found for gene ", x$Gene.name,", omitting from analysis.\n", sep='')
          next}
      }
    }
    exprTable[gN,samp] <- exprTable[gN,samp]+ epTable[i, 'ImmuneScore'] 
  }
  return(exprTable)
}



# Analysis of CRCmseq -----------------------------------------------------


geneTable <- read.table('~/Dropbox/Code/mart_export_hg19.txt', sep='\t', header=T)
names(geneTable) <- c('Gene.ID', 'Start', 'End', 'Chrom', 'Gene.name')

dir <- '~/CRCdata/CRCmseq_Set'
epTable <- read.table(paste0(dir, '/Neopred_results/Output_BA.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-24), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]

#Plot score distribution for a selected sample
sample <- 'Set.04.snv'
eps <- subsetEpTable(epTable, sample)
ggplot(eps, aes(x = Score)) + geom_histogram()


#Create immunescore
epTable$ImmuneScore <- epTable$Score

#Assign score to genes
epTable$Gene.name <- sapply(epTable$Gene, function(x) unlist(strsplit(x, ':'))[1]  )
exprTable <- getGeneScores(epTable, geneTable)
exprTableFiltered <- exprTable[rowSums(exprTable > 0)>2,]



# TCGA analysis -----------------------------------------------------------
geneTable <- read.table('~/Dropbox/Code/mart_export.txt', sep='\t', header=T)
names(geneTable) <- c('Gene.ID', 'Start', 'End', 'Chrom', 'Gene.name')


dir <- 'TCGA_CRC'
epTable <- read.table(paste0(dir, '/Neopred_results/Output_BA.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-24), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]

#Plot score distribution for a selected sample
sample <- "TCGA-AG-3881.MSS"
eps <- subsetEpTable(epTable, sample)
ggplot(eps, aes(x = Affinity)) + geom_histogram()
ggplot(eps[eps$Affinity<500,], aes(x = Score)) + geom_histogram()

#excluded ones: hypermutated patients
excl <- c('TCGA-AA-3977.MSS', 'TCGA-AG-A002')
epTable <- subset(epTable, !(Sample %in% excl))

#Create immunescore
epTable$ImmuneScore <- epTable$Affinity

#Assign score to genes
epTable$Gene.name <- sapply(epTable$Gene, function(x) unlist(strsplit(x, ':'))[1]  )
exprTable <- getGeneScores(epTable, geneTable)
exprTableFiltered <- exprTable[rowSums(exprTable > 0)>2,]

