#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Get statistics for Polyp and Set ---------------------------------------------------

setwd('~/CRCdata')


colReds = c('#fee0d2','#fc9272','#de2d26')
colBlues = c('#deebf7','#9ecae1','#3182bd')

mutRatiosBatch = list()
mutRatioTable = data.frame(matrix(vector(), ncol=6))
names(mutRatioTable) <- c('Clonal_All', 'Clonal_Ep','Private_All','Private_Ep','Shared_All', 'Shared_Ep')

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Set')
#dirList <- c('IBD')

for (dir in dirList){
  
pdf(paste0(dir, '_summary_allbinding.pdf'), height = 5, width=8)
epTable <- read.table(paste0(dir, '/Neopred_results/Output.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epTable <- epTable[epTable$Novelty!=0,]
epTable <- epTable[epTable$Sample!='Set.10.snv',] #to disregard the replication of Set10

# epTable <- filterByBAPrediction(dir, epTable)

#epTable <- filterByWTBinding(dir, epTable)

summaryTable <- read.table(paste0(dir,'/Neopred_results/Output.neoantigens.summarytable.txt'), header=T, row.names=1)
summaryTable <- summaryTable[unique(epTable$Sample),]

#Adjust summary table in case neo-epitopes have been filtered
summaryTable <- recalculateSummaryTable(epTable, summaryTable, mutations = F)
summaryTableMut <- recalculateSummaryTable(epTable, summaryTable)

clonality <- rbind(melt(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')])),
                   melt(as.matrix(summaryTableMut[,c('Clonal', 'Shared','Subclonal')])))
clonality$Type <- c(rep('P',3*nrow(summaryTable)), rep('M', 3*nrow(summaryTable)))
pc <- ggplot(data=clonality, aes(x=Type, fill=Var2)) + geom_bar(aes(y=value),stat="identity", position = position_fill(reverse=T)) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=rev(colReds)) + labs(x = "Tumour", y = "", fill="Clonality") + ggtitle("Clonality of neoepitopes") +
  facet_grid(. ~ Var1)
print(pc)

eps <- melt(as.matrix(summaryTable[, 'Total', drop=F]/summaryTableMut[,'Total',drop=F]))
pe <- ggplot(data=eps, aes(x=Var1, y=value)) + geom_bar(stat='identity') + scale_fill_manual(values=colBlues[1]) +
  labs(x = "Tumour", y = "") + ggtitle("Average neo-epitopes per neo-ep mutation") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(pe)


summaryTableMut$Total_MUT <- sapply(row.names(summaryTableMut), function(x) getTotalMutFromFasta(dir, x))
pem <- ggplot(data=summaryTableMut, aes(x=row.names(summaryTableMut), y=Total/Total_MUT)) + geom_bar(stat='identity', fill=colBlues[3]) +
  scale_y_continuous(labels = percent_format()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Tumour", y = "") + ggtitle("Percentage of missense mutations producing at least one neo-epitope")
print(pem)

mutRatioTable <- getMutationRatios(dir, epTable, mutRatioTable)
mutRatiosBatch[dir] <- list(summaryTableMut$Total/summaryTableMut$Total_MUT)
dev.off()
}


mutRatiosBatch['Random_proteome'] <- list(random.summary$EpMuts/random.summary$AllMuts)


mutRatioTable$Clonal_Ratio <- mutRatioTable$Clonal_Ep/mutRatioTable$Clonal_All
mutRatioTable$Private_Ratio <- mutRatioTable$Private_Ep/mutRatioTable$Private_All
mutRatioTable$Subclonal_Ratio <- (mutRatioTable$Private_Ep+mutRatioTable$Shared_Ep)/(mutRatioTable$Private_All+mutRatioTable$Shared_All)
mutRatioTable[mutRatioTable$Private_Ep<12,'Private_Ratio'] <- NA

# Analyse Polyp, Set and Random -------------------------------------------

ggqqplot(mutRatiosBatch[[2]])

t.test(mutRatiosBatch[[1]], mutRatiosBatch[[2]])

#Plot with ggplot

mycolors = c("#1984c0","#b32200", "#999999", "#19a0c0","#e23e00", "#e28a00")

mRB <- data.frame('ratio' = unlist(mutRatiosBatch))
mRB$set <- as.factor(c(rep('Adenoma', 5), rep('Carcinoma', 11), rep('Random', 50)))

p <- ggplot(mRB, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, fill='black') +
  theme_classic() + scale_fill_manual(values=mycolors[1:3])

mRClonal <- data.frame('ratio' = c(mutRatioTable$Clonal_Ratio, mutRatioTable$Subclonal_Ratio, unlist(mutRatiosBatch[[3]])))
mRClonal$set <- as.factor(c(rep('Ad-Clonal', 5), rep('Car-Clonal', 11), rep('Ad-Subclonal', 5), rep('Car-Subclonal', 11), rep('Random', 49)))

p2 <- ggplot(mRClonal, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='black') +
  theme_classic() + scale_fill_manual(values=mycolors[c(1,4,2,5,3)])

carcClonal <- data.frame('ratio' = c(mutRatioTable$Clonal_Ratio[6:16], mutRatioTable$Subclonal_Ratio[6:16], unlist(mutRatiosBatch[[3]])))
carcClonal$set <- as.factor(c(rep('Clonal', 11), rep('Subclonal', 11), rep('Random', 49)))

p3 <- ggplot(carcClonal, aes(x=set, y=ratio, fill=set)) + geom_violin() +
  scale_x_discrete(limits=c("Clonal", "Subclonal", "Random")) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, fill='black') +
  theme_classic() + scale_fill_manual(values=mycolors[c(2,3,5)])



p4 = ggpaired(carcClonal[1:22,], x='set', y='ratio', color='set', line.color='gray', palette=mycolors[c(2,5)], line.size=0.4)

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

