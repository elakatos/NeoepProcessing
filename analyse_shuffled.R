setwd('~/RNAseq/Neoepitopes/Randomised_profiles/Neopred_results')

getTotalMutFromFasta <- function(dir, sample){
  #returns only mutations that were inputted into netMHC
  fFile <- paste0(dir, '/fastaFiles/',sample,'.tmp.10.fasta')
  fData <- scan(file=fFile, what='string')
  return(sum(grepl('>', fData)))
}


mainSummaryTable <- data.frame(matrix(vector()))
thresh.Reco <- 1e-1

for (i in 1:14){
epTable <- read.table(paste0('Total_',i,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epTable <- subset(epTable, Novelty==1)

st <- data.frame(Sample=unique(epTable$Sample), stringsAsFactors = F)
st$Type <- unlist(strsplit(st$Sample[1],'\\.'))[1]

epTable.dedup <- epTable[!duplicated(paste0(epTable$Sample,epTable$LineID)),]
st$Binding_mutation <- sapply(st$Sample, function(x) sum(epTable.dedup$Sample==x))
epTable.strong <- subset(epTable, BindLevel=='SB')
epTable.strong.dedup <- epTable.strong[!duplicated(paste0(epTable.strong$Sample, epTable.strong$LineID)),]
st$Strong_binding_mutation <- sapply(st$Sample, function(x) sum(epTable.strong.dedup$Sample==x))

recoTable <- read.table(paste0('recopo_',i,'/PredictedRecognitionPotentials.txt'),
                        sep='\t', stringsAsFactors = F, header=T)
recoTable$AntigenID <- paste0(recoTable$Sample,':',recoTable$Mutation,':',recoTable$MutantPeptide)
epTable$AntigenID <- paste0(epTable$Sample, ':',epTable$Identity, ':', epTable$peptide)
recoTable.imm <- subset(recoTable, NeoantigenRecognitionPotential>thresh.Reco)
epTable.imm <- subset(epTable, AntigenID %in% recoTable.imm$AntigenID)
epTable.imm.dedup <- epTable.imm[!duplicated(paste0(epTable.imm$Sample, epTable.imm$LineID)),]

st$Neoantigen_mutation <- sapply(st$Sample, function(x) sum(epTable.imm.dedup$Sample==x))

st$Total_missense_mutation <- sapply(st$Sample, function(x) getTotalMutFromFasta('.', x))

mainSummaryTable <- rbind(mainSummaryTable, st)
}

write.table(mainSummaryTable,file='Total.summary.txt',sep='\t', quote=F, row.names=F)



# Plot outputs ------------------------------------------------------------

st.df <- read.table('~/RNAseq/Neoepitopes/Randomised_profiles/Neopred_results/Total.summary.txt', sep='\t', header=T, stringsAsFactors = F)

st.df$Supertype <- sapply(st.df$Type, function(x) unlist(strsplit(x, '_'))[1])

st.sub <- subset(st.df, Type %in% c('MSI_average','MSS_average'))
ggplot(st.sub, aes(x=Type, fill=Type, y=Neoantigen_mutation/Total_missense_mutation)) +
  geom_violin() + geom_boxplot(fill='grey70',width=0.05) +
  theme_mypub() + guides(fill=F) +
  labs(y='Proportion of antigenic mutations',x='') +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=c('darkorange3','darkseagreen4','grey50')) +
  stat_compare_means(comparisons=list(c('MSI_average','MSS_average')))

st.sub <- subset(st.df, Type %in% c('TcellLow_average','TcellMedium_average','TcellHigh_average'))
st.sub$Type <- factor(st.sub$Type, c('TcellLow_average','TcellMedium_average','TcellHigh_average'))
ggplot(st.sub, aes(x=Type, fill=Type, y=Neoantigen_mutation/Total_missense_mutation)) +
  geom_violin() + geom_boxplot(fill='grey70',width=0.05) +
  theme_mypub() + guides(fill=F) +
  labs(y='Proportion of antigenic mutations',x='T-cell score') +
  scale_x_discrete(labels=c('Low_average', 'Medium_average', 'High_average')) +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  stat_compare_means(comparisons=list(c('TcellMedium_average','TcellLow_average'),
                                      c('TcellMedium_average','TcellHigh_average'),
                                      c('TcellLow_average','TcellHigh_average')))


st.sub <- subset(st.df, !grepl('average',Type))
st.sub$Supertype <- factor(st.sub$Supertype, c('TcellLow','TcellMedium','TcellHigh'))
ggplot(st.sub, aes(x=Supertype, fill=Supertype, y=Neoantigen_mutation/Binding_mutation)) +
  geom_violin() + geom_boxplot(fill='grey70',width=0.05) +
  theme_mypub() + guides(fill=F) +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  stat_compare_means(comparisons=list(c('TcellMedium','TcellLow'),
                                      c('TcellMedium','TcellHigh'),
                                      c('TcellLow','TcellHigh')))



# Check WT vs MT binding --------------------------------------------------

i <- 3
recoTable <- read.table(paste0('recopo_',i,'/PredictedRecognitionPotentials.txt'),
                        sep='\t', stringsAsFactors = F, header=T)
bindingTable <- read.table(paste0('recopo_',i,'/Neoantigens.WTandMTtable.txt'),
                           sep='\t',stringsAsFactors = F, header=T)
recoTable$MT.binding <- bindingTable[match(recoTable$NeoantigenID, bindingTable$ID),'MT.SCORE']
recoTable$WT.binding <- bindingTable[match(recoTable$NeoantigenID, bindingTable$ID),'WT.SCORE']

ggplot(recoTable, aes(x=WT.binding, y=MT.binding)) +
  geom_point(alpha=0.05) + scale_x_log10() + scale_y_log10() +
  geom_abline(intercept = 0, slope=1, alpha=0.7, colour='firebrick') +
  geom_vline(xintercept = 500, alpha=0.7, colour='firebrick') +
  theme_mypub()

table(recoTable$WT.binding>500)
table(recoTable$WT.binding>500, recoTable$NeoantigenRecognitionPotential>1e-1)
table(recoTable$WT.binding>500, recoTable$NeoantigenRecognitionPotential>1)

