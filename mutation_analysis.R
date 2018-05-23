
dir <- 'IBD'
epTable <- read.table(paste0(dir, '/Neopred_results/Output.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]

epTable <- filterByWTBinding(dir, epTable, 'a')

epTableStrong <- epTable[epTable$BindLevel=='SB',]

vafTestDF <- data.frame(matrix(vector(), ncol=4))
names(vafTestDF) <- c('Sample', 'Region', 'UB', 'pValue')


for (sample in unique(epTable$Sample)){
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
exonic <- readExonicFile(sampleFileEx)
eps <- subsetEpTable(epTableStrong, sample, unique=T)

tumorColumns <- grep('Region*', names(exonic))
isEpMutation <- (exonic$LineID %in% eps$LineID)
uB <- 1.0
lB <- 0.0
pdf(paste0(dir,'_',sample,'.pdf'),height=5,width=8)
par(mfrow=c(1,2))
for (i in tumorColumns){
  allVafs <- computeVafAD(exonic, i)
  epVafs <- computeVafAD(exonic[isEpMutation,], i)
  nonepVafs <- computeVafAD(exonic[!isEpMutation,], i)
  
  allVafsF <- allVafs[(allVafs>lB) & (allVafs< uB)]
  epVafsF <- epVafs[(epVafs>lB) & (epVafs < uB)]
  nonepVafsF <- nonepVafs[(nonepVafs>lB) & (nonepVafs < uB)]
  
  qqplot(nonepVafsF, epVafsF, pch=19, xlab='All mutations', ylab='Neoepitope mutations', main='QQplot')
  plot.ecdf(nonepVafsF,col='black', ylab='CDF', xlab='VAF')
  plot.ecdf(epVafsF,col='grey50', add=T)
  print(sample)
  print(ks.test(nonepVafsF, epVafsF, alternative='less'))
}
dev.off()
}


#VAF plotting of selected
sample = 'Oxford_IBD2.mutectCalls..somatic'
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
exonic <- readExonicFile(sampleFileEx)
eps <- subsetEpTable(epTableStrong, sample, unique=T)

isEpMutation <- (exonic$LineID %in% eps$LineID)

uB<-0.7
lB<-0.025

i=22+0
allVafs <- computeVafAD(exonic, i)
epVafs <- computeVafAD(exonic[isEpMutation,], i)
nonepVafs <- computeVafAD(exonic[!isEpMutation,], i)
allVafsF <- allVafs[(allVafs>lB) & (allVafs< uB)]
epVafsF <- epVafs[(epVafs>lB) & (epVafs < uB)]
nonepVafsF <- nonepVafs[(nonepVafs>lB) & (nonepVafs < uB)]
vafDF <- data.frame('vaf'=allVafsF, 'type'='All mutations')
epDF <- data.frame('vaf'=epVafsF, 'type'='Neo-epitope mutations')
nonepDF <- data.frame('vaf'=nonepVafsF, 'type'='Non neo-epitope mutations')
DF <- rbind(nonepDF, epDF)

mycols = c('#d0cc9e','#4165d1', '#e13512')


ggplot(nonepDF, aes(x=vaf, y=..scaled..)) + geom_density(fill='grey30',alpha=0.4, adjust=0.8) +
  geom_density(data=epDF, aes(x=vaf, y=..scaled..), alpha=0.4, fill='red', adjust=1)

# p1 = ggplot(DF, aes(x=vaf, y=..scaled.., fill = type)) + geom_density(alpha=0.5, adjust=1) +
#   scale_fill_manual(values=mycols)
# 
# pdf('~/Dropbox/Neoepitopes/Prelim_vaf_S10:R2.pdf', width=8, height=5)
# pl <- p1+scale_x_continuous(limits=c(0.01, 0.7)) + theme_bw() +
#   theme(text = element_text(size=16 ,family='sans'), legend.position = c(0.8, 0.7), legend.background = element_rect(colour='black')) +
#   labs(x='Variant allele frequency', y='Frequency') + guides(fill=guide_legend(title=NULL))
# print(pl)
# dev.off()


# Example VAF generation + plotting

clonalVAFs <- rbinom(5, 40, 0.45)/40
subclonalVAfs <- vector()
for (i in 1:5){
  genVAFs <- rbinom(5*(2^i), 50, (0.5/(2^i)+rnorm(1,0,0.04)))/50
  subclonalVAfs <- c(subclonalVAfs, genVAFs)
}

vafDF <- data.frame('vaf' = c(clonalVAFs, subclonalVAfs))
vafDF <- vafDF[vafDF$vaf>0.025,,drop=F]
ggplot(vafDF_saved_ep, aes(x=vaf, y=..scaled..)) + geom_density(adjust=1)


p1 = ggplot(vafDF_saved, aes(x=vaf, y=..scaled..)) + geom_density(fill='#889174', adjust=0.8)

pdf('~/Dropbox/Neoepitopes/Example_VAF_i.pdf', width=6, height=5)
pl <- p1+scale_x_continuous(limits=c(0.02, 0.7)) + scale_y_continuous(breaks=c(0.25,0.75)) +
  theme_bw() + theme(text = element_text(size=20 ,family='sans')) +
  labs(x='', y='') + guides(fill=guide_legend(title=NULL))
print(pl)
dev.off()

# Epitope distribution ----------------------------------------------------

tumorColumns = grep('Region*', names(eps))
hist(eps$Rank, breaks=20)
hist(eps[rowSums(eps[, tumorColumns])==4,]$Rank, breaks=20  )

epRankClonal <- eps[rowSums(eps[, tumorColumns])==4,]$Rank
epRankNotClonal <- eps[rowSums(eps[, tumorColumns])<4,]$Rank





# random.dataBA <- read.table('random_proteome_all_BA.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)
# random.data.realBA <- subset(random.dataBA, !((nchar(random.dataBA$peptide)==9) &  (random.dataBA$peptide_pos) %in% c(1,11)) )
# random.data.realBA <- subset(random.data.realBA, Novelty==1)
# random.data.realBA$Sample <- random.data.realBA$PatIndex
# random.data.filteredBA <- subset(random.data.realBA, BindLevel!='N')
# 
# random.data.filtered <- getSharedEps(random.data.filtered, random.data.filteredBA)


lines(density(random.summary$EpMuts/random.summary$AllMuts), col='red')
barplot(random.summary$Epitopes/random.summary$EpMuts)
barplot(random.summary$SB/random.summary$Epitopes)

hist(random.data.filtered[random.data.filtered$PatIndex==14,]$Rank, breaks=20)

# HLAs --------------------------------------------------------------------

cat(hlasA, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt')
cat('\n', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat(hlasB, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat('\n', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
cat(hlasC, sep=' ', file='~/RNAseq/Software/NeoepProcessing/HLA_list.txt', append=T)
