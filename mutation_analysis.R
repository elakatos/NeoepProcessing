
epTable <- read.table(paste0(dir, '/Neopred_results/Output.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epNon <- epTable[epTable$Novelty==0,]
epTable <- epTable[epTable$Novelty!=0,]
#barplot(table(epNon$Sample)/table(epTable$Sample)*100, las=2)
epTableStrong <- epTable[epTable$BindLevel=='SB',]

sample = 'Oxford_IBD1.mutectCalls..somatic'
#sampleFile <- paste0(dir, '/avready/',sample,'.avinput')
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
#avinput <- readAvinput(sampleFile)
exonic <- readExonicFile(sampleFileEx)
eps <- subsetEpTable(epTable, sample, unique=T)

tumorColumns <- grep('Region*', names(exonic))
isEpMutation <- (exonic$LineID %in% eps$LineID)
uB <- 0.3
lB <- 0.05
pdf(paste0(dir,':',sample,'_s.pdf'),height=5,width=8)
par(mfrow=c(1,2))
for (i in tumorColumns){
  allVafs <- computeVafAD(exonic, i)
  epVafs <- computeVafAD(exonic[isEpMutation,], i)
  
  allVafsF <- allVafs[(allVafs>lB) & (allVafs< uB)]
  epVafsF <- epVafs[(epVafs>lB) & (epVafs < uB)]
  
  print(length(allVafsF))
  qqplot(allVafsF, epVafsF, pch=19, xlab='All mutations', ylab='Neoepitope mutations', main='QQplot')
  plot.ecdf(allVafsF,col='black', ylab='CDF', xlab='VAF')
  plot.ecdf(epVafsF,col='grey50', add=T)
  print(ks.test(allVafsF, epVafsF))
}
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
