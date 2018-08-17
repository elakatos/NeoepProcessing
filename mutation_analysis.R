
dir <- '../TCGA_CRC'
prefix <- 'Total'
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel')
#epNon <- epTable[epTable$Novelty==0,]
#epTable <- epTable[epTable$Novelty!=0,]
#epTable <- filterByWTBinding(dir, epTable, 'a')
#epTableStrong <- epTable[epTable$BindLevel=='SB',]

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

# Entire VAF of CRCmseq samples

vcf <- read.table('~/RNAseq/Neoepitopes/CRCmseq_Set/avready/Set.10.snv.avinput',sep='\t', stringsAsFactors = F)
names(vcf)[c(1:5, 18:ncol(vcf))] <- c('chr','start','end','ref','alt',getRegionNames(ncol(vcf)-18,T))

v <- getRegionNames(ncol(vcf)-18)
vaf.data <- data.frame(matrix(vector(), ncol=length(v),nrow=nrow(vcf)))
vaf.data[,1:length(v)] <- sapply(v, function(z) computeVaf(vcf,z))

ggplot(melt(vaf.data), aes(x=value, fill=variable)) + geom_histogram(alpha=0.6, position='dodge')
ggplot(vaf.data[vaf.data$X4>(-0)& vaf.data$X4<0.8 ,], aes(x=X4)) + geom_histogram(bins=40) +
  geom_histogram(data=vaf.data[rowSums(vaf.data==0)==0,], aes(x=X4), fill='red', alpha=0.5,bins=40)

#ggplot(vaf.data[rowSums(vaf.data==0)==0,], aes(x=X13)) + geom_histogram(bins=30)

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


###########################################################################
# TCGA analysis -----------------------------------------------------------

dir <- '../TCGA_CRC'
prefix <- 'Total'
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-23), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel')

epTable$mutID <- apply(epTable, 1,function(x) paste0(x['Sample'], ':',x['LineID'] ) )

#epTable <- subset(epTable, Affinity<200)


allMutVAFs <- data.frame(matrix(vector(), ncol=3)); names(allMutVAFs) <- c('Sample', 'LineID', 'VAF')

for(sample in unique(epTable$Sample)){
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
exonic <- readExonicFile(sampleFileEx)

tmpVAF <- data.frame(Sample=sample, LineID=exonic$LineID, VAF=computeVafAD(exonic, 22))
allMutVAFs <- rbind(allMutVAFs, tmpVAF)
}

allMutVAFs$mutID <- apply(allMutVAFs, 1,function(x) paste0(x['Sample'], ':',x['LineID'] ) )

pdf(paste0(dir, '_VAFs_collection_rp.pdf'),width=8,height=5)

fmax=0.8; fmin=0.1
steps <- seq(fmax,fmin,by=(-1e-2))

for (samp in row.names(subset(summaryTable, (Total>200)))){
#for (samp in row.names(subset(summaryTable, (Total>120) & !(MUT) & !(LOH) & !(B2M) & !(PDL)))){
  
pAll <- ggplot(subset(allMutVAFs, (Sample==samp)), aes(x=VAF)) + geom_histogram(bins=30)
pEp <- ggplot(subset(allMutVAFs, (Sample==samp) & (mutID %in% epTable$mutID)), aes(x=VAF, y=..density..)) + geom_histogram(bins=30, fill='firebrick3', alpha=0.5) +
  labs(title=samp) + scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  geom_histogram(data=subset(allMutVAFs, (Sample==samp) & !(mutID %in% epTable$mutID)), aes(x=VAF,y=..density..), bins=30, fill='skyblue4', alpha=0.5)

vafAll<-subset(allMutVAFs, (Sample==samp))$VAF
vafEp <- subset(allMutVAFs, (Sample==samp) & (mutID %in% epTable$mutID))$VAF

cumvaf <- data.frame(invf = (1/steps), cumvaf=sapply(steps, function(x) sum(vafAll>=x))/length(vafAll),
                     cumvafEp=sapply(steps, function(x) sum(vafEp>=x))/length(vafEp))
pCum <- ggplot(cumvaf, aes(x=invf, y=cumvaf)) + geom_line() + geom_line(data=cumvaf, aes(x=invf, y=cumvafEp), colour='red') +
  labs(x='1/f', y='Cumulative frequency', title=summaryTable[samp,'MSI'])

grid.arrange(pEp, pCum, nrow=1)

}
dev.off()



# VAF vs epitope-strength -------------------------------------------------

samp <- 'TCGA-AD-6889'
d <- epTable.dedup[epTable.dedup$Sample==samp,]
d$VAF <- allMutVAFs[match(d$mutID, allMutVAFs$mutID),'VAF']
ggplot(d, aes(x=Rank, y=VAF)) + geom_point() + scale_y_continuous(limits=c(0, 0.25)) + stat_cor()

