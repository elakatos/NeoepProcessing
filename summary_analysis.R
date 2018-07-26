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
prefix='Total'
#pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel')
#epNon <- epTable[epTable$Novelty==0,]
#epTable <- epTable[epTable$Novelty!=0,]
#epTable <- epTable[epTable$Affinity<500,]
epTable <- epTable[epTable$BindLevel=='SB',]

summaryTableMut <- processSummaryOfSampleSet(dir, epTable, prefix)
summaryTableMut$MutRatio <- summaryTableMut$Total/summaryTableMut$Total_MUT
summaryTableMut <- summaryTableMut[order(summaryTableMut$MutRatio),]

summaryTableMut$MSI <- clin.df[match(row.names(summaryTableMut), clin.df$Patient), 'MSI']
summaryTableMut$MSI[summaryTableMut$MSI=='MSS'] <- 'MSI-L'
summaryTableMut$Hyp <- clin.df[match(row.names(summaryTableMut), clin.df$Patient), 'Hypermut']
summaryTableMut$B2M <- clin.df[match(row.names(summaryTableMut), clin.df$Patient), 'B2M']>0
summaryTableMut$VS <- clin.df[match(row.names(summaryTableMut), clin.df$Patient), 'VS']
summaryTableMut$Age <- clin.df[match(row.names(summaryTableMut), clin.df$Patient), 'Age']
summaryTableMut$Race <- clin.df[match(row.names(summaryTableMut), clin.df$Patient), 'Race']
#summaryTableMut$MSI <- clin.df.verified[match(row.names(summaryTableMut), clin.df.verified$Patient), 'MSI']
#summaryTableMut$Hyp <- clin.df.verified[match(row.names(summaryTableMut), clin.df.verified$Patient), 'Hypermut']

# Compare neoep-mutation ratio in MSI and MSS
summaryTableMut.sub <- subset(summaryTableMut, MSI != 'POLE')
summaryTableMut.sub <- subset(summaryTableMut, Total_MUT > 20)
ggplot(summaryTableMut.sub, aes(x=Total_MUT, y=MutRatio, colour=Hyp)) + geom_point() + scale_x_continuous(trans='log1p')
ggplot(summaryTableMut.sub, aes(x=MSI, y=MutRatio, fill=as.factor(MSI))) + geom_violin() +
stat_compare_means(comparisons = list(c('MSI-H', 'MSI-L')))
ggplot(summaryTableMut.sub, aes(x=Hyp, y=MutRatio, fill=as.factor(Hyp))) + geom_violin() +
  stat_compare_means(comparisons = list(c(1,0)))
ggplot(summaryTableMut.sub, aes(x=B2M, y=MutRatio, fill=B2M)) + geom_violin() +
  stat_compare_means(comparisons = list(c(FALSE,TRUE)))
ggplot(summaryTableMut.sub, aes(x=VS, y=MutRatio, fill=as.factor(VS))) + geom_violin() +
  stat_compare_means(comparisons = list(c('alive','dead')))

summaryTableMut.sub$Age <- cut(summaryTableMut.sub$Age, c(min(summaryTableMut.sub$Age, na.rm=T),
                                                          median(summaryTableMut.sub$Age, na.rm=T),
                                                          max(summaryTableMut.sub$Age, na.rm=T)))
summaryTableMut.sub$Age <- cut(summaryTableMut.sub$Age, quantile(summaryTableMut.sub$Age, na.rm=T))
ggplot(summaryTableMut.sub, aes(x=Age, y=MutRatio, fill=as.factor(Age))) + geom_violin() +
  stat_compare_means(comparisons = list(c('(31.2,57.7]','(75.6,90.1]')))
summaryTableMut.sub$Race[summaryTableMut.sub$Race=='not reported'] <- NA
ggplot(summaryTableMut.sub, aes(x=Race, y=MutRatio, fill=Race)) + geom_violin() +
  stat_compare_means(comparisons = list(c('black or african american','asian')))

#Lowest/highest mutation ratios?



#Check if there are shared epitopes/mutations
epTable$Gene.name <- sapply(epTable$Gene, function(x) unlist(strsplit(x, ':'))[1]  )
epTable$mutation <- apply(epTable, 1,function(x) paste0(x['Chrom'], ':',x['Start'] ) )
epTable.mut.dedup <- epTable[!duplicated(epTable[, c('Sample', 'LineID')]),]
epTable.pep.dedup <- epTable[!duplicated(epTable[, c('Sample', 'peptide')]),]

mutation.table <- table(epTable.mut.dedup$mutation)
mutation.table.sorted <- mutation.table[order(-mutation.table)]

epitope.table <- table(epTable.pep.dedup$peptide)
epitope.table.sorted <- epitope.table[order(-epitope.table)]
