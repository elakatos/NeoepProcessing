#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Get statistics for Polyp and Set ---------------------------------------------------

setwd('~/CRCdata')

analysisPostfix <- 'WTallfilteredStrong'

mutRatiosBatch = list()
mutRatioTable = data.frame(matrix(vector(), ncol=6))
names(mutRatioTable) <- c('Clonal_All', 'Clonal_Ep','Private_All','Private_Ep','Shared_All', 'Shared_Ep')

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Set')
#dirList <- c('IBD')
prefix <- 'Output_BA'

for (dir in dirList){
  
  pdf(paste0(dir, '_summary_',analysisPostfix,'.pdf'), height = 5, width=8)
  epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                        sep = '\t',stringsAsFactors = F, fill=T)
  names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-24), 'LineID', 'Chrom', 'Start',
                      'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                      'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
  epTable <- epTable[epTable$Novelty!=0,]
  #epTable <- epTable[epTable$Sample!='Set.10.recalled.snv',] #to disregard the replication of Set10
  
  #Filter epTable according to WT or alternative binding prediction
  #epTable <- filterByBAPrediction(dir, epTable)
  #epTable <- filterByWTBinding(dir, epTable, 'all')
  #epTable <- epTable[epTable$BindLevel=='SB',]
  
  summaryTableMut <- processSummaryOfSampleSet(dir, epTable, prefix)
  mutRatioTable <- getMutationRatios(dir, epTable, mutRatioTable)
  mutRatiosBatch[dir] <- list(summaryTableMut$Total/summaryTableMut$Total_MUT)
  dev.off()
}


# Read in random data

random.data <- read.table('random_proteome_all.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)
random.data <- subset(random.data, !((nchar(random.data$peptide)==9) &  (random.data$peptide_pos) %in% c(1,11)) )
random.data <- subset(random.data, Novelty==1)
random.data$Sample <- random.data$PatIndex
random.data.filtered <- subset(random.data, BindLevel!='N')

#Filter according to WT or other binding information
random.data.filtered <- filterRandomByWTBinding(random.data.filtered, 'a')
random.data.filtered <- random.data.filtered[random.data.filtered$BindLevel=='SB',]

random.summary <- processSummaryOfRandomSet(random.data, random.data.filtered)

mutRatiosBatch['Random'] <- list(random.summary$EpMuts/random.summary$AllMuts)

#Compute percentage of neo-ep associated mutations for mutations of specific clonality
mutRatioTable$Clonal_Ratio <- mutRatioTable$Clonal_Ep/mutRatioTable$Clonal_All
mutRatioTable$Private_Ratio <- mutRatioTable$Private_Ep/mutRatioTable$Private_All
mutRatioTable$Subclonal_Ratio <- (mutRatioTable$Private_Ep+mutRatioTable$Shared_Ep)/(mutRatioTable$Private_All+mutRatioTable$Shared_All)
mutRatioTable[mutRatioTable$Private_Ep<10,'Private_Ratio'] <- NA

# Analyse Polyp, Set and Random -------------------------------------------

#are the mutation ratios normally distributed?
ggqqplot(mutRatiosBatch[[2]])

mycolors = c(colReds, colBlues,"#999999")
nA <- 5
nC <- 11

#Plot percentage of neo-ep mutations
mRB <- data.frame('ratio' = unlist(mutRatiosBatch))
mRB$set <- as.factor(c(rep('Adenoma', nA), rep('Carcinoma', nC), rep('Random', 50)))
mycomp <- list( c("Adenoma", "Carcinoma"), c("Adenoma", "Random"), c('Carcinoma', 'Random') )

p1 <- ggplot(mRB, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill='black') +
  scale_fill_manual(values=mycolors[c(3,6,7)]) +
  stat_compare_means(label = "p.signif",comparisons=mycomp, hide.ns = T)

mRClonal <- data.frame('ratio' = c(mutRatioTable$Clonal_Ratio, mutRatioTable$Subclonal_Ratio, unlist(mutRatiosBatch[['Random']])))
mRClonal$set <- as.factor(c(rep('Ad-Clonal', nA), rep('Car-Clonal', nC), rep('Ad-Subclonal', nA), rep('Car-Subclonal', nC), rep('Random', 50)))
mycomp = list(c('Ad-Clonal', 'Ad-Subclonal'),c('Car-Clonal', 'Car-Subclonal'),c('Ad-Clonal', 'Car-Clonal'), c('Car-Clonal', 'Random'), c('Ad-Clonal', 'Random') )

p2 <- ggplot(mRClonal, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='black') +
  scale_fill_manual(values=mycolors[c(6,5,3,2,7)]) +
  stat_compare_means(label = "p.signif",comparisons=mycomp, hide.ns = T) #label.y = c(0.9, 0.97, 1.01, 1.01, 1.07)

carcClonal <- subset(mRClonal, set %in% c('Car-Clonal', 'Car-Subclonal', 'Random'))
#p3 <- ggplot(carcClonal, aes(x=set, y=ratio, fill=set)) + geom_violin() +
#  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='black') +
#  scale_fill_manual(values=mycolors[c(3,2,7)])


p4 = ggpaired(carcClonal[carcClonal$set!='Random',], x='set', y='ratio', fill='set',line.color='gray40',
              palette=mycolors[c(2,5)], line.size=0.4, point.size=2, ggtheme=theme_gray(), alpha=0.4) +
  stat_compare_means(paired=T)

pdf(paste0('Neoepitope_mutations_', analysisPostfix, '.pdf'), width=8, height=5)
print(p1);print(p2);print(p4)
dev.off()



# TCGA sample -------------------------------------------------------------

dir = 'TCGA_CRC'
prefix='Output_BA'
#pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-24), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]
epTable <- epTable[epTable$Affinity<500,]

summaryTableMut <- processSummaryOfSampleSet(dir, epTable, prefix)
summaryTableMut$MutRatio <- summaryTableMut$Total/summaryTableMut$Total_MUT
summaryTableMut$MSI <- sapply(row.names(summaryTableMut), function(x) substr(x, 14,16))

ggplot(summaryTableMut, aes(x=Total_MUT, y=MutRatio, colour=MSI)) + geom_point()

#Check if there are shared epitopes/mutations
epTable$Gene.name <- sapply(epTable$Gene, function(x) unlist(strsplit(x, ':'))[1]  )
epTable$mutation <- apply(epTable, 1,function(x) paste0(x['Chrom'], ':',x['Start'] ) )
epTable.mut.dedup <- epTable[!duplicated(epTable[, c('Sample', 'LineID')]),]
epTable.pep.dedup <- epTable[!duplicated(epTable[, c('Sample', 'peptide')]),]

mutation.table <- table(epTable.mut.dedup$mutation)
mutation.table.sorted <- mutation.table[order(-mutation.table)]

epitope.table <- table(epTable.pep.dedup$peptide)
epitope.table.sorted <- epitope.table[order(-epitope.table)]

# HLA types ---------------------------------------------------------------

hlaConvert <- function(x){
  if (is.na(x)){
    y <- NA
  }
  else{
    y <- paste0('HLA-', toupper(substr(x,5,5)),substr(x,7,8),':',substr(x,10,11))
  }
  return(y)
}

nnGrep <- function(x){
  i <- gregexpr("HLA", x) #find the two HLAs mentioned
  origHLA <- substr(x, i[[1]][1], i[[1]][1]+9)
  nn <- substr(x, i[[1]][2], i[[1]][2]+9)
  j <- regexpr("[0-9]\\.[0-9][0-9][0-9]", x) #match the format of distance
  dist <- substr(x, j, j+4)
  return(list(HLA=origHLA, NN=nn, Distance=dist))
}

#Read in, filter out possibly incorrect predictions and convert to netMHCpan format
#dir = 'TCGA_CRC'
#hlas <- read.table(paste0('~/RNAseq/Neoepitopes/',dir,'/hlatypes.txt'), sep='\t', header=T, row.names=1, stringsAsFactors = F)
hlasOrig <- read.table('~/Dropbox/Code/TCGA/hlatypes_tcga_1.txt', sep='\t', header=T, row.names=1, stringsAsFactors = F)
hlasCorrect <- subset(hlasOrig, !( (HLA.A_1=='hla_a_01_01_01_01') & (is.na(HLA.A_2)) & (HLA.B_1== 'hla_b_07_02_01') & (is.na(HLA.B_2)) & (HLA.C_1=='hla_c_01_02_01') & (is.na(HLA.C_2))  ))
hlas <- as.data.frame(apply(hlasCorrect, 2,function(r) sapply(r, function(x) hlaConvert(x))) )
hlaList <- as.vector(as.matrix(hlas))

#Build nearest neighbour table from encountered HLA types
hlaNN <- data.frame(matrix(vector(), ncol=3))
names(hlaNN) <- c('HLA', 'NN', 'Distance')
hlaMaps <- scan(file='~/Dropbox/Code/TCGA/hla_mappings.txt', what=character(), sep='\n') #Collection of lines stating nearest neighbours
for (i in 1:length(hlaMaps)){ hlaNN[i,] <- nnGrep(hlaMaps[i]) }
hlaNN <- hlaNN[!duplicated(hlaNN$HLA),]

hlaList.mapped <- mapvalues(hlaList, from=hlaNN$HLA, to=hlaNN$NN)
hlas.mapped <- as.data.frame(sapply(hlas, function(x) mapvalues(x, from=hlaNN$HLA, to=hlaNN$NN)))

#todo <- setdiff(hlaList, hlaNN$HLA)
#todo <- sapply(todo, function(x) paste0(substr(x, 1, 7),substr(x, 9, 10)))

#Get patients who have a particular HLAtype and NA all other hla
allele<-'HLA-A02:01'
nonAllele <- hlas.mapped!=allele

hlasOut <- hlasCorrect; hlasOut[nonAllele] <- NA
hlasOut <- hlasOut[rowSums(is.na(hlasOut))<6,]
#Further filtering: exclude known MSI and hyper-mutated ones
nonHyper <- subset(clin.df, (MSI %in% c('MSS', NA)) &(Hypermut %in% c(0, NA)) )
hlasOut <- hlasOut[row.names(hlasOut) %in% nonHyper$Patient,]
write.table(hlasOut, paste0('~/Dropbox/Code/TCGA/hlatypes_',sub(':','', allele),'.txt'), sep='\t', quote=F)
