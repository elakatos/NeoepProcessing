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
  
  summaryTable <- processSummaryOfSampleSet(dir, epTable, prefix)
  mutRatioTable <- getMutationRatios(dir, epTable, mutRatioTable)
  mutRatiosBatch[dir] <- list(summaryTable$Total/summaryTable$Total_MUT)
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

dir = '../TCGA_CRC'
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
#epTable <- epTable[epTable$BindLevel=='SB',]

summaryTable <- processSummaryOfSampleSet(dir, epTable, prefix)

summaryTable$MutRatio <- summaryTable$Total/summaryTable$Total_MUT
summaryTable <- summaryTable[order(summaryTable$MutRatio),]

summaryTable$MSI <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'MSI']
summaryTable$Hyp <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Hypermut']
summaryTable$B2M <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'B2M']>0
summaryTable$VS <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'VS']
summaryTable$Age <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Age']/365
summaryTable$Race <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Race']
summaryTable$Stage <- gsub('[abc]', '', clin.df[match(row.names(summaryTable), clin.df$Patient), 'Stage'])
summaryTable[summaryTable=='not reported'] <- NA

summaryTable$LOH <- lohhla.patients[match(row.names(summaryTable), lohhla.patients$ID), 'LOH']
summaryTable$MUT <- lohhla.patients[match(row.names(summaryTable), lohhla.patients$ID), 'MUT']
#summaryTable$MSI <- clin.df.verified[match(row.names(summaryTable), clin.df.verified$Patient), 'MSI']
#summaryTable$Hyp <- clin.df.verified[match(row.names(summaryTable), clin.df.verified$Patient), 'Hypermut']

# Compare neoep-mutation ratio in MSI and MSS
summaryTable.sub <- subset(summaryTable, MSI != 'POLE')
summaryTable.sub <- subset(summaryTable.sub, Total_MUT > 30)
ylb <- 'RATIO of neo-epitope associated mutations'

summaryTable.sub$MutRatio <- summaryTable.sub$Total #Instead look at the total number of neo-epitope mutations
ylb <- 'TOTAL neo-epitope associated mutations'

p0 <- ggplot(summaryTable.sub, aes(x=Total_MUT, y=MutRatio, colour=as.factor(Hyp))) + geom_point() + scale_x_continuous(trans='log10') +scale_y_continuous(trans='log10') +
  stat_cor(label.x.npc = 'centre') + stat_cor(mapping=aes(x=Total_MUT,y=MutRatio,colour=Hyp),label.x.npc='centre', label.y.npc='bottom') +
  scale_color_manual(values=c(rgb(0.2,0.35,0.45), rgb(0.35,0.55,0.8))) +
  labs(colour='Hypermutated', x='Total somatic missense mutations', y=ylb)
p1 <- ggplot(summaryTable.sub, aes(x=MSI, y=MutRatio, fill=as.factor(MSI))) + geom_violin() +
stat_compare_means(comparisons = list(c('MSI-H', 'MSI-L'))) +
  guides(fill=F) + labs(y=ylb, title='All non-POLE tumours')
p2 <- ggplot(summaryTable.sub, aes(x=Hyp, y=MutRatio, fill=as.factor(Hyp))) + geom_violin() +
  stat_compare_means(comparisons = list(c(1,0))) +
  guides(fill=F) + labs(y=ylb, x='Hypermutated', title='All non-POLE tumours')
p3 <- ggplot(summaryTable.sub, aes(x=B2M, y=MutRatio, fill=B2M)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in B2M', title='All non-POLE tumours')
p4 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$MUT),], aes(x=MUT, y=MutRatio, fill=MUT)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in HLA', title='All non-POLE tumours')
p5 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$LOH),], aes(x=LOH, y=MutRatio, fill=as.factor(LOH))) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='LOH in HLA', title='All non-POLE tumours')
p6 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$VS),], aes(x=VS, y=MutRatio, fill=as.factor(VS))) + geom_violin() +
  stat_compare_means(comparisons = list(c('alive','dead'))) +
  guides(fill=F) + labs(y=ylb, x='Vital status', title='All non-POLE tumours')
p7 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Stage),], aes(x=Stage, y=MutRatio, fill=as.factor(Stage))) + geom_violin() +
  #stat_compare_means(comparisons = list(c('stge i', 'stge ii'), c('stge ii', 'stge iii'), c('stge iii', 'stge iv'),
  #                                      c('stge i', 'stge iii'), c('stge ii', 'stge iv'), c('stge i', 'stge iv')))
  guides(fill=F) + labs(y=ylb, title='All non-POLE tumours')
summaryTable.sub$Age <- cut(summaryTable.sub$Age, quantile(summaryTable.sub$Age, na.rm=T))
p8 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Age),], aes(x=Age, y=MutRatio, fill=as.factor(Age))) + geom_violin() +
  guides(fill=F) + labs(y=ylb,title='All non-POLE tumours')
p9 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Race),], aes(x=Race, y=MutRatio, fill=Race)) + geom_violin() +
  #stat_compare_means(comparisons = list(c('black or african american','white'), c('black or african american','asian'), c('asian','white'))) +
  guides(fill=F) + labs(y=ylb,title='All non-POLE tumours')

summaryTable.sub2 <- subset(summaryTable.sub, !LOH & (MSI=='MSI-H'))
p42 <- ggplot(summaryTable.sub2[!is.na(summaryTable.sub2$MUT),], aes(x=MUT, y=MutRatio, fill=MUT)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in HLA', title='MSI tumours without LOH')

summaryTable.sub2 <- subset(summaryTable.sub, (MSI == 'MSI-L'))
p52 <- ggplot(summaryTable.sub2[!is.na(summaryTable.sub2$LOH),], aes(x=LOH, y=MutRatio, fill=LOH)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='LOH in HLA', title='MSS tumours')


pdf('~/Dropbox/Code/TCGA/Figures/Mutneoep_filtered_stats.pdf', width=6.5, height=5)
#print(p0); print(p1); print(p2); print(p3); print(p4); print(p5);print(p6); print(p7); print(p8); print(p9); print(p42); print(p52)

print(p0);print(p1+scale_y_log10());  print(p2+scale_y_log10());  print(p3+scale_y_log10());  print(p4+scale_y_log10());  print(p5+scale_y_log10());  print(p6+scale_y_log10());  print(p7+scale_y_log10());  print(p8+scale_y_log10());  print(p9+scale_y_log10());  print(p42+scale_y_log10());  print(p52+scale_y_log10());
dev.off()


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
