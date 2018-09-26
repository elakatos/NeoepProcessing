#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Get statistics for Polyp and Set ---------------------------------------------------

setwd('~/CRCdata')

analysisPostfix <- 'all'

mutRatiosBatch = list()
mutRatioTable = data.frame(matrix(vector(), ncol=6))
names(mutRatioTable) <- c('Clonal_All', 'Clonal_Ep','Private_All','Private_Ep','Shared_All', 'Shared_Ep')

dirList <- c('CRCmseq_Polyp', 'CRCmseq_Ha','CRCmseq_Set')
#dirList <- c('IBD')
prefix <- 'Output_BA'

for (dir in dirList){
  
  #pdf(paste0(dir, '_summary_',analysisPostfix,'.pdf'), height = 5, width=8)
  epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                        sep = '\t',stringsAsFactors = F, fill=T)
  names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-24), 'LineID', 'Chrom', 'Start',
                      'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                      'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
  epTable <- epTable[epTable$Novelty!=0,]
  epTable <- epTable[!(epTable$Sample %in% c('Set.10.recalled.snv', 'Set.10.snv', 'Set.01.snv')),] #to disregard the replication of Set10
  
  #Filter epTable according to WT or alternative binding prediction
  #epTable <- filterByBAPrediction(dir, epTable)
  #epTable <- filterByWTBinding(dir, epTable, 'all')
  #epTable <- epTable[epTable$BindLevel=='SB',]
  recoTable <- read.table(paste0(dir, '/Neopred_results/PredictedRecognitionPotentials.txt'),
                          sep='\t', stringsAsFactors = F, header=T)
  recoTable$AntigenID <- apply(recoTable, 1,function(x) paste0(x['Sample'], ':',x['Mutation'], ':',x['MutantPeptide'] ) )
  epTable$AntigenID <- apply(epTable, 1,function(x) paste0(x['Sample'], ':',x['Identity'], ':',x['peptide'] ) )
  recoTable.imm <- subset(recoTable, NeoantigenRecognitionPotential>1e-4)
  epTable <- subset(epTable, AntigenID %in% recoTable.imm$AntigenID)
  
  summaryTable <- processSummaryOfSampleSet(dir, epTable, prefix)
  mutRatioTable <- getMutationRatios(dir, epTable, mutRatioTable)
  mutRatiosBatch[dir] <- list(summaryTable$Total/summaryTable$Total_MUT)
  #dev.off()
}


# Read in random data

random.data <- read.table('random_proteome_all.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)
random.data <- subset(random.data, !((nchar(random.data$peptide)==9) &  (random.data$peptide_pos) %in% c(1,11)) )
#random.data <- subset(random.data, Novelty==1)
random.data$Sample <- random.data$PatIndex
random.data.filtered <- subset(random.data, BindLevel!='N')

random.data.filtered$AntigenID <- apply(random.data.filtered, 1,function(x) paste0(x['Sample'], ':',x['Identity'], ':',x['peptide'] ) )

random.recopo <- read.table('random.PredictedRecognitionPotentials.txt', sep='\t', header=T, stringsAsFactors = F)
random.recopo$AntigenID <- apply(random.recopo, 1,function(x) paste0(x['Sample'], ':',x['Mutation'], ':',x['MutantPeptide'] ) )
random.recopo.imm <- subset(random.recopo, NeoantigenRecognitionPotential>1e-4)
random.data.filtered <- subset(random.data.filtered, AntigenID %in% random.recopo.imm$AntigenID)

#Filter according to WT or other binding information
#random.data.filtered <- filterRandomByWTBinding(random.data.filtered, 'a')
#random.data.filtered <- random.data.filtered[random.data.filtered$BindLevel=='SB',]

random.summary <- processSummaryOfRandomSet(random.data, random.data.filtered)

mutRatiosBatch['Random'] <- list(random.summary$EpMuts/random.summary$AllMuts)

#Compute percentage of neo-ep associated mutations for mutations of specific clonality
mutRatioTable$Clonal_Ratio <- mutRatioTable$Clonal_Ep/mutRatioTable$Clonal_All
mutRatioTable$Private_Ratio <- mutRatioTable$Private_Ep/mutRatioTable$Private_All
mutRatioTable$Subclonal_Ratio <- (mutRatioTable$Private_Ep+mutRatioTable$Shared_Ep)/(mutRatioTable$Private_All+mutRatioTable$Shared_All)
mutRatioTable[mutRatioTable$Private_Ep<10,'Private_Ratio'] <- NA

# Analyse Polyp, Set and Random -------------------------------------------

mr.df <- data.frame(x=unlist(mutRatiosBatch), var=c(rep('Adenoma',8), rep('Carcinoma',9),
                                                    rep('Random', 50)))
ggplot(mr.df, aes(x=var, y=x, fill=var)) + geom_violin() +
  geom_boxplot(width=0.05, fill='grey80') +
  stat_compare_means(comparisons=list(c('Adenoma', 'Carcinoma'), c('Adenoma', 'Random'), c('Carcinoma', 'Random')))

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
#epTable <- epTable[epTable$BindLevel=='SB',]

summaryTable <- processSummaryOfSampleSet(dir, epTable, prefix)

summaryTable$MutRatio <- summaryTable$Total/summaryTable$Total_MUT
summaryTable$NeoepRatio <- summaryTable$Neoep/summaryTable$Total

summaryTable$MSI <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'MSI']
summaryTable$Hyp <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Hypermut']
summaryTable$B2M <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'B2M']>0
summaryTable$VS <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'VS']
summaryTable$Age <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Age']/365
summaryTable$Race <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Race']
summaryTable$Cancer <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Cancer']
summaryTable$Stage <- gsub('[abc]', '', clin.df[match(row.names(summaryTable), clin.df$Patient), 'Stage'])
summaryTable[summaryTable=='not reported'] <- NA

summaryTable$LOH <- lohhla.patients[match(row.names(summaryTable), lohhla.patients$ID), 'LOH']
summaryTable$MUT <- sapply(row.names(summaryTable), function(x) x %in% muthla.signif$individual)
summaryTable$PDL <- lohhla.patients[match(row.names(summaryTable), lohhla.patients$ID), 'PDL']
summaryTable$CYT <- lohhla.patients[match(row.names(summaryTable), lohhla.patients$ID), 'CYT']

immTable <- read.table('~/Dropbox/Code/TCGA/data_CRC_fullimmune.txt',stringsAsFactors = F, header=T, sep='\t')
summaryTable$IP <- immTable[match(row.names(summaryTable), immTable$sampleID), 'Immune.phenotype']
#summaryTable$MSI <- clin.df.verified[match(row.names(summaryTable), clin.df.verified$Patient), 'MSI']
#summaryTable$Hyp <- clin.df.verified[match(row.names(summaryTable), clin.df.verified$Patient), 'Hypermut']

# Compare neoep-mutation ratio in MSI and MSS
summaryTable.sub <- subset(summaryTable, MSI != 'POLE')
summaryTable.sub <- subset(summaryTable.sub, Total_MUT > 30)
mlab <- 'All non-POLE tumours'
ylb <- 'RATIO of neo-epitope associated mutations'

summaryTable.sub$MutRatio <- summaryTable.sub$NeoepRatio
ylb <- 'AVG NUMBER of neo-epitopes from one mutation'

summaryTable.sub$MutRatio <- summaryTable.sub$Total #Instead look at the total number of neo-epitope mutations
ylb <- 'TOTAL neo-epitope associated mutations'

summaryTable.sub$MutRatio <- summaryTable.sub$Neoep
ylb <- 'TOTAL neo-epitopes'

summaryTable.sub$MutRatio <- summaryTable.sub$Total_MUT
ylb <- 'TOTAL somatic missense mutations'

p0 <- ggplot(summaryTable.sub, aes(x=Total_MUT, y=MutRatio, colour=as.factor(Hyp))) + geom_point(size=2) + scale_x_continuous(trans='log10') +
  stat_cor(label.x.npc = 'centre', size=5) + stat_cor(mapping=aes(x=Total_MUT,y=MutRatio,colour=''),label.x.npc='centre', label.y.npc='bottom', size=5) +
  scale_color_manual(values=c('black','darkseagreen4', 'darkorange3')) +
  labs(colour='Hypermutated', x='Total somatic missense mutations', y=ylb) + theme_bw() + guides(color=F)
p1 <- ggplot(summaryTable.sub, aes(x=MSI, y=MutRatio, fill=as.factor(MSI))) + geom_violin() +
  stat_compare_means(comparisons = list(c('MSI-H', 'MSI-L')), aes(label = paste0("p = ", ..p.format..))) +
  guides(fill=F) + labs(y=ylb, title=mlab) + theme_bw() + scale_fill_manual(values=c('darkorange3', 'darkseagreen4'))
p2 <- ggplot(summaryTable.sub, aes(x=Hyp, y=MutRatio, fill=as.factor(Hyp))) + geom_violin() +
  stat_compare_means(comparisons = list(c(1,0))) +
  guides(fill=F) + labs(y=ylb, x='Hypermutated', title=mlab)
p3 <- ggplot(summaryTable.sub, aes(x=B2M, y=MutRatio, fill=B2M)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in B2M', title=mlab)
p4 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$MUT),], aes(x=MUT, y=MutRatio, fill=MUT)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in HLA', title=mlab)
p5 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$LOH),], aes(x=LOH, y=MutRatio, fill=as.factor(LOH))) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='LOH in HLA', title=mlab)
p6 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$VS),], aes(x=VS, y=MutRatio, fill=as.factor(VS))) + geom_violin() +
  stat_compare_means(comparisons = list(c('alive','dead'))) +
  guides(fill=F) + labs(y=ylb, x='Vital status', title=mlab)
p7 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Stage),], aes(x=Stage, y=MutRatio, fill=as.factor(Stage))) + geom_violin() +
  #stat_compare_means(comparisons = list(c('stge i', 'stge ii'), c('stge ii', 'stge iii'), c('stge iii', 'stge iv'),
  #                                      c('stge i', 'stge iii'), c('stge ii', 'stge iv'), c('stge i', 'stge iv'))) +
  guides(fill=F) + labs(y=ylb, title=mlab)
#summaryTable.sub$Age <- cut(summaryTable.sub$Age, quantile(summaryTable.sub$Age, na.rm=T))
#p8 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Age),], aes(x=Age, y=MutRatio, fill=as.factor(Age))) + geom_violin() +
#  guides(fill=F) + labs(y=ylb,title=mlab)
p8 <- ggplot(summaryTable.sub, aes(x=Age, y=MutRatio, colour=MSI)) + geom_point(size=2) + stat_cor(mapping=aes(x=Age,y=MutRatio,colour=Age),size=5) +
  guides(colour=F) + labs(y=ylb,title=mlab) + theme_bw() + scale_color_manual(values=c('darkseagreen4', 'darkorange3'))
p9 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Race),], aes(x=Race, y=MutRatio, fill=Race)) + geom_violin() +
  stat_compare_means(comparisons = list(c('black or african american','white'), c('black or african american','asian'), c('asian','white'))) +
  guides(fill=F) + labs(y=ylb,title=mlab)
p10 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$IP),], aes(x=as.factor(IP), y=MutRatio, fill=as.factor(IP))) + geom_violin() +
  #stat_compare_means(comparisons = list(c('4','6'), c('1','4'), c('3','4'))) +
  guides(fill=F) + labs(y=ylb,title=mlab)
p11 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Cancer),], aes(x=Cancer, y=MutRatio, fill=Cancer)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TCGA-COAD','TCGA-READ'))) +
  guides(fill=F) + labs(y=ylb,title=mlab)
p12 <- ggplot(summaryTable.sub, aes(x=PDL, y=MutRatio, colour=-PDL)) + geom_point() + stat_cor() +
  guides(colour=F) + labs(y=ylb,title=mlab)
p13 <- ggplot(summaryTable.sub, aes(x=CYT, y=MutRatio, colour=-CYT)) + geom_point() + stat_cor() +
  guides(colour=F) + labs(y=ylb,title=mlab)

summaryTable.sub2 <- subset(summaryTable.sub, !LOH & (MSI=='MSI-H'))
p42 <- ggplot(summaryTable.sub2[!is.na(summaryTable.sub2$MUT),], aes(x=MUT, y=MutRatio, fill=MUT)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in HLA', title='MSI tumours without LOH')

# summaryTable.sub2 <- subset(summaryTable.sub, !MUT & (MSI=='MSI-H'))
# p43 <- ggplot(summaryTable.sub2[!is.na(summaryTable.sub2$LOH),], aes(x=LOH, y=MutRatio, fill=LOH)) + geom_violin() +
#   stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
#   guides(fill=F) + labs(y=ylb, x='LOH in HLA', title='MSI tumours')

summaryTable.sub2 <- subset(summaryTable.sub, (MSI == 'MSI-L'))
p52 <- ggplot(summaryTable.sub2[!is.na(summaryTable.sub2$LOH),], aes(x=LOH, y=MutRatio, fill=LOH)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='LOH in HLA', title='MSS tumours')


pdf('~/Dropbox/Code/TCGA/Figures/Mutratio_recopo_stats_STAD.pdf', width=6.5, height=5)
print(p0); print(p1); print(p2); print(p3); print(p4); print(p5);#print(p12); print(p13);
print(p6); print(p7); print(p9); #print(p10);print(p11);
print(p8); print(p42); print(p52)

#print(p0+scale_y_log10());print(p1+scale_y_log10());  print(p2+scale_y_log10());  print(p3+scale_y_log10());  print(p4+scale_y_log10());  print(p5+scale_y_log10());print(p12+scale_y_log10());print(p13+scale_y_log10());  print(p6+scale_y_log10());  print(p7+scale_y_log10());  print(p9+scale_y_log10()); print(p10+scale_y_log10()); print(p11+scale_y_log10());print(p8+scale_y_log10());  print(p42+scale_y_log10());  print(p52+scale_y_log10());
dev.off()


pcrcrat <- p1 + geom_boxplot(width=0.05, fill='grey80') + labs(y='Relative neoantigen burden', title='',x='') +
  theme(text=element_text(size=16))

pcrcage1 <- ggplot(summaryTable.sub, aes(x=Age, y=MutRatio, colour=MSI)) + geom_point(size=2) + stat_cor(mapping=aes(x=Age,y=MutRatio,colour=Age),size=5) +
  guides(colour=F) + labs(y='Relative neoantigen burden') + theme_bw() + scale_color_manual(values=c('darkorange3','darkseagreen4')) +
  scale_y_continuous(limits=c(0.0,0.25)) + scale_x_continuous(breaks=c(30,60,90)) +
  theme(text=element_text(size=16))

pcrcage2 <- ggplot(summaryTable.sub, aes(x=Age, y=Neoep, colour=MSI)) + geom_point(size=2) + stat_cor(mapping=aes(x=Age,y=Neoep,colour=Age),size=5) +
  guides(colour=F) + labs(y='Total neoantigen burden') + theme_bw() + scale_color_manual(values=c('darkorange3','darkseagreen4')) +
  scale_y_log10(breaks=c(100,1000)) + scale_x_continuous(breaks=c(30,60,90)) +
  theme(text=element_text(size=16))

pcrcage3 <- ggplot(summaryTable.sub, aes(x=Age, y=Total_MUT, colour=MSI)) + geom_point(size=2) + stat_cor(mapping=aes(x=Age,y=Total_MUT,colour=Age),size=5) +
  guides(colour=F) + labs(y='Total missense mutation burden') + theme_bw() + scale_color_manual(values=c('darkorange3','darkseagreen4')) +
  scale_y_log10(breaks=c(100,1000)) + scale_x_continuous(breaks=c(30,60,90)) +
  theme(text=element_text(size=16))

pstadrat <- p1 + geom_boxplot(width=0.05, fill='grey80') + labs(y='Relative neoantigen burden', title='',x='') +
  theme(text=element_text(size=16))

#Lowest/highest mutation ratios?



#Check if there are shared epitopes/mutations
epTable$Gene.name <- sapply(epTable$Gene, function(x) unlist(strsplit(x, ':'))[1]  )
epTable$mutation <- apply(epTable, 1,function(x) paste0(x['Chrom'], ':',x['Start'], ':',x['AltAll'] ) )
epTable.mut.dedup <- epTable[!duplicated(epTable[, c('Sample', 'LineID')]),]
epTable.pep.dedup <- epTable[!duplicated(epTable[, c('Sample', 'peptide')]),]

mutation.table <- table(epTable.mut.dedup$mutation)
mutation.table.sorted <- mutation.table[order(-mutation.table)]

epitope.table <- table(epTable.pep.dedup$peptide)
epitope.table.sorted <- epitope.table[order(-epitope.table)]

gene.table <- table(epTable.mut.dedup$Gene.name)
gene.table.sorted <- gene.table[order(-gene.table)]

#Compare neoep mutations with generic mutations

mutTable <- data.frame(matrix(vector()))

for (samp in unique(epTable$Sample)){
muts <- tryCatch(
  {read.table(paste0(dir,'/VCF/',samp,'.vcf'),sep='\t', stringsAsFactors = F)[,c(1,2,4,5)]},
  error=function(e){return(NA)}
)
if(!is.na(muts)){
  muts$Sample <- samp
  mutTable <- rbind(mutTable, muts)
}
}
names(mutTable)[1:4] <- c('Chrom', 'Start', 'RefAll', 'AltAll')

mutTable$mutation <- apply(mutTable, 1,function(x) paste0(x['Chrom'], ':',x['Start'], ':',x['AltAll'] ) )
all.table <- table(mutTable$mutation); all.table.sorted <- all.table[order(-all.table)]

mutCompare <- data.frame(matrix(vector(), nrow=length(mutation.table)))
row.names(mutCompare) <- names(mutation.table)
mutCompare$Neoep <- mutation.table[match(row.names(mutCompare), names(mutation.table))]
mutCompare$All <- all.table[match(row.names(mutCompare), names(all.table))]

mutCompare.signif <- subset(mutCompare, All>25)


############################################################################
# Recognition Potential ---------------------------------------------------

epTable$AntigenID <- apply(epTable, 1,function(x) paste0(x['Sample'], ':',x['Identity'], ':',x['peptide'] ) )


recoTable <- read.table('~/RNAseq/Neoepitopes/TCGA_CRC/Neopred_results/PredictedRecognitionPotentials.txt',
                        stringsAsFactors = F, header=T)

recoTable$AntigenID <- apply(recoTable, 1,function(x) paste0(x['Sample'], ':',x['Mutation'], ':',x['MutantPeptide'] ) )

recoTable.imm <- subset(recoTable, NeoantigenRecognitionPotential>1e-4)

epTable.imm <- subset(epTable, AntigenID %in% recoTable.imm$AntigenID)

recoTable.sorted <- recoTable.imm[order(-recoTable.imm$NeoantigenRecognitionPotential),]
recoTable.pat1 <- recoTable.sorted[!duplicated(recoTable.sorted$Sample),]

ggplot(recoTable.pat1, aes(x=NeoantigenRecognitionPotential)) +geom_histogram(bins=40)

############################################################################
# Add expression ----------------------------------------------------------

expr <- read.table('~/Dropbox/Code/TCGA/CRC_RNA_expression.tpm', sep='\t', header=T, row.names=1, stringsAsFactors = F)
gene.key <- read.table('~/Dropbox/Code/mart_export.txt', sep='\t', header=T)

gene.med = apply(expr, 1, median)
names(gene.med) <- sapply(names(gene.med), function(z) unlist(strsplit(z, '\\.'))[1])
gene.key$Median <- gene.med[match(gene.key$Gene.stable.ID, names(gene.med))]
gene.key <- gene.key[order(-gene.key$Median),]

epTable$Gene.median <- gene.key[match(epTable$Gene.name, gene.key$Gene.name), 'Median']

epTable$Gene.expr <- apply(epTable, 1, function(x)
  expr[gene.key[match(x['Gene.name'], gene.key$Gene.name),'Gene.stable.ID'],gsub('-', '.',x['Sample'])])

