#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Get statistics for Polyp and Set ---------------------------------------------------

setwd('~/CRCdata')

analysisPostfix <- 'WTallfiltered'

mutRatiosBatch = list()
mutRatioTable = data.frame(matrix(vector(), ncol=6))
names(mutRatioTable) <- c('Clonal_All', 'Clonal_Ep','Private_All','Private_Ep','Shared_All', 'Shared_Ep')

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Set')
#dirList <- c('IBD')

for (dir in dirList){
  
  pdf(paste0(dir, '_summary_',analysisPostfix,'.pdf'), height = 5, width=8)
  epTable <- read.table(paste0(dir, '/Neopred_results/Output.neoantigens.txt'), header=F,
                        sep = '\t',stringsAsFactors = F, fill=T)
  names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                      'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                      'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
  epTable <- epTable[epTable$Novelty!=0,]
  epTable <- epTable[epTable$Sample!='Set.10.recalled.snv',] #to disregard the replication of Set10
  
  #Filter epTable according to WT or alternative binding prediction
  #epTable <- filterByBAPrediction(dir, epTable)
  epTable <- filterByWTBinding(dir, epTable)
  
  summaryTableMut <- processSummaryOfSampleSet(dir, epTable)
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
random.data.filtered <- filterRandomByWTBinding(random.data.filtered)

random.summary <- processSummaryOfRandomSet(random.data, random.data.filtered)

mutRatiosBatch['Random_proteome'] <- list(random.summary$EpMuts/random.summary$AllMuts)

#Compute percentage of neo-ep associated mutations for mutations of specific clonality
mutRatioTable$Clonal_Ratio <- mutRatioTable$Clonal_Ep/mutRatioTable$Clonal_All
mutRatioTable$Private_Ratio <- mutRatioTable$Private_Ep/mutRatioTable$Private_All
mutRatioTable$Subclonal_Ratio <- (mutRatioTable$Private_Ep+mutRatioTable$Shared_Ep)/(mutRatioTable$Private_All+mutRatioTable$Shared_All)
mutRatioTable[mutRatioTable$Private_Ep<10,'Private_Ratio'] <- NA

# Analyse Polyp, Set and Random -------------------------------------------

#are the mutation ratios normally distributed?
ggqqplot(mutRatiosBatch[[2]])

mycolors = c(colReds, colBlues,"#999999")

#Plot percentage of neo-ep mutations
mRB <- data.frame('ratio' = unlist(mutRatiosBatch))
mRB$set <- as.factor(c(rep('Adenoma', 5), rep('Carcinoma', 11), rep('Random', 50)))
mycomp <- list( c("Adenoma", "Carcinoma"), c("Adenoma", "Random"), c('Carcinoma', 'Random') )

p1 <- ggplot(mRB, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill='black') +
  scale_fill_manual(values=mycolors[c(3,6,7)]) +
  stat_compare_means(label = "p.signif",comparisons=mycomp, hide.ns = T)

mRClonal <- data.frame('ratio' = c(mutRatioTable$Clonal_Ratio, mutRatioTable$Subclonal_Ratio, unlist(mutRatiosBatch[[3]])))
mRClonal$set <- as.factor(c(rep('Ad-Clonal', 5), rep('Car-Clonal', 11), rep('Ad-Subclonal', 5), rep('Car-Subclonal', 11), rep('Random', 50)))
mycomp = list(c('Ad-Clonal', 'Ad-Subclonal'),c('Car-Clonal', 'Car-Subclonal'),c('Ad-Clonal', 'Car-Clonal'), c('Car-Clonal', 'Random'), c('Ad-Clonal', 'Random') )

p2 <- ggplot(mRClonal, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='black') +
  scale_fill_manual(values=mycolors[c(6,5,3,2,7)]) +
  stat_compare_means(label = "p.signif",comparisons=mycomp, hide.ns = T, label.y = c(0.9, 0.97, 1.01, 1.01, 1.07))

carcClonal <- data.frame('ratio' = c(mutRatioTable$Clonal_Ratio[6:16], mutRatioTable$Subclonal_Ratio[6:16], unlist(mutRatiosBatch[[3]])))
carcClonal$set <- as.factor(c(rep('Clonal', 11), rep('Subclonal', 11), rep('Random', 50)))

p3 <- ggplot(carcClonal, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  scale_x_discrete(limits=c("Clonal", "Subclonal", "Random")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='black') +
  scale_fill_manual(values=mycolors[c(2,3,5)])



p4 = ggpaired(carcClonal[1:22,], x='set', y='ratio', fill='set',line.color='gray40',
              palette=mycolors[c(2,5)], line.size=0.4, point.size=2, ggtheme=theme_gray(), alpha=0.4) +
  stat_compare_means(paired=T)

pdf(paste0('Neoepitope_mutations_', analysisPostfix, '.pdf'), width=8, height=5)
print(p1);print(p2);print(p4)
dev.off()

# TCGA sample -------------------------------------------------------------

dir = 'TCGA_COAD'
pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/TCGA_COAD.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
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

