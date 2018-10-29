

# Get Random tables -------------------------------------------------------

random.data <- read.table('random_proteome_all.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)
random.data.wt <- read.table('random_wt_proteome_all.txt', sep='\t',row.names=NULL, header=T,stringsAsFactors = F)

to.exclude <- ((nchar(random.data$peptide)==9) &  (random.data$peptide_pos) %in% c(1,11))
random.data <- random.data[!to.exclude,]
random.data.wt <- random.data.wt[!to.exclude,]

no.binder <- random.data$Affinity > 500
random.data <- random.data[!no.binder,]
random.data.wt <- random.data.wt[!no.binder,]

mutandwt.df <- data.frame(ID=1:nrow(random.data), MUTATION_ID=random.data$Identity,
                          Sample=random.data$PatIndex,
                          WT.PEPTIDE=random.data.wt$peptide, MT.PEPTIDE=random.data$peptide,
                          MT.ALLELE=random.data$hla,
                          WT.SCORE=random.data.wt$Affinity, MT.SCORE=random.data$Affinity,
                          HLA='hla', CHOP_SCORE=1)

write.table(mutandwt.df,file='~/RNAseq/Neoepitopes/random.Neoantigens.WTandMTtable.txt', quote=F, row.names=F, sep='\t')

# VAF ---------------------------------------------------------------------


dir <- '../CRCmseq_Set'
prefix <- 'Set6_WGS'
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample', getRegionNames(ncol(epTable)-24), 'LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity', 'Rank', 'Cand', 'BindLevel','Novelty')
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
sample = 'Set.06.WGS.snv'
sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
exonic <- readExonicFile(sampleFileEx)
eps <- subsetEpTable(epTableStrong, sample, unique=T)

isEpMutation <- (exonic$LineID %in% eps$LineID)

uB<-0.7
lB<-0.025

i=22+0
allVafs <- computeVaf(exonic, i, 5, 6)
epVafs <- computeVaf(exonic[isEpMutation,], i, 5, 6)
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

epvaf <- data.frame(vafs=allVafs, eps = isEpMutation)
ggplot(epvaf[epvaf$vafs>0,], aes(x=vafs, fill=eps)) + geom_histogram(alpha=0.5, position='identity')

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

# Entire VAF of CRCmseq samples -------------------------------------------------

vcf <- read.table('~/CRCdata/CRCmseq_Set/avready/Set.06.WGS.snv.avinput',sep='\t', stringsAsFactors = F)
names(vcf)[c(1:5, 18:ncol(vcf))] <- c('chr','start','end','ref','alt',getRegionNames(ncol(vcf)-18,F),'Normal')

v <- getRegionNames(ncol(vcf)-18)
vaf.data <- data.frame(matrix(vector(), ncol=length(v),nrow=nrow(vcf)))
vaf.data[,1:length(v)] <- sapply(v, function(z) computeVaf(vcf,z, 5, 6)) #compute vafs with indices defined for NR and NV

ggplot(melt(vaf.data), aes(x=value, fill=variable)) + geom_histogram(alpha=0.6, position='dodge',bins=50)
ggplot(vaf.data[vaf.data$X4>(-0)& vaf.data$X4<0.8 ,], aes(x=X4)) + geom_histogram(bins=80) #+
  geom_histogram(data=vaf.data[rowSums(vaf.data==0)==0,], aes(x=X4), fill='red', alpha=0.5,bins=40)

f <- seq(0.4,0.05, by=-1e-3)
cumvaf <- sapply(f, function(x) sum((vaf.data[,4]*1.25)>x))
cumvaf.df <- data.frame(invf=1/f, cumvaf=cumvaf)
ggplot(cumvaf.df, aes(x=invf, y=cumvaf)) + geom_line()
#ggplot(vaf.data[rowSums(vaf.data==0)==0,], aes(x=X13)) + geom_histogram(bins=30)

# Epitope distribution ----------------------------------------------------

tumorColumns = grep('Region*', names(eps))
hist(eps$Rank, breaks=20)
hist(eps[rowSums(eps[, tumorColumns])==4,]$Rank, breaks=20  )

epRankClonal <- eps[rowSums(eps[, tumorColumns])==4,]$Rank
epRankNotClonal <- eps[rowSums(eps[, tumorColumns])<4,]$Rank

#Based on recognition potential

recoTable <- read.table('~/CRCdata/CRCmseq_Set/Neopred_results/PredictedRecognitionPotentials.txt',
                        stringsAsFactors = F, header=T)

recoTable.imm <- recoTable[recoTable$NeoantigenRecognitionPotential>1e-4,]
ggplot(subset(recoTable.imm, Sample=='Set.09.Distal.snv' & A<30), aes(x=NeoantigenRecognitionPotential)) + geom_density()


epSave <- epTable[,c('Score', 'Sample', 'Rank', 'Affinity')]
pp <- ggplot(epSave, aes(x=Score, colour=Sample)) + geom_density(adjust=1.2, size=1.1) +guides(col=F) +theme_bw()
ggplot(subset(recoTable.imm,  A<200), aes(x=NeoantigenRecognitionPotential, colour=Sample)) + geom_density() + scale_x_log10()

eptest <- epTable[epTable$Sample=='Set.09.Proximal.snv',]
regno <- sum(eptest[1,2:14]!=-1)
eptest.cl <- eptest[rowSums(eptest[,2:(1+regno)])==regno,]
eptest.sc <- eptest[rowSums(eptest[,2:(1+regno)])==1,]

ggplot(eptest.sc, aes(x=Score)) + geom_density(col='skyblue4') + geom_density(data=eptest.cl, aes(x=Score), col='darkred') +
  labs(title=paste0('p = ',ks.test(eptest.sc$Score, eptest.cl$Score)$p.value,'\n(One-sided: p = ',ks.test(eptest.sc$Score, eptest.cl$Score,alternative='less')$p.value,')'))


r.cl <- subset(recoTable.imm, Sample==eptest$Sample[1] & (Mutation %in% eptest.cl$Identity))
r.sc <- subset(recoTable.imm, Sample==eptest$Sample[1] & (Mutation %in% eptest.sc$Identity))

ggplot(r.sc, aes(x=A)) + geom_density(col='skyblue4') + geom_density(data=r.cl, aes(x=A), col='darkred') +
  labs(title=paste0('p = ',ks.test(r.sc$A, r.cl$A)$p.value,'\n(One-sided: p = ',ks.test(r.sc$A, r.cl$A,alternative='less')$p.value,')'))


exonic <- readExonicFile(paste0(dir, '/avannotated/','Set.09.Distal.snv','.avannotated.exonic_variant_function'))


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




# Subclonal neo-epitopes and LOH ------------------------------------------

samp <- 'STM003.merged'
LOHregs <- paste0('Region_', c(3,4,2))
noLOHregs <- paste0('Region_',c(1))

x <- subset(epTable, Sample == samp)
noregs <- sum(x[1,2:10]!=-1)
x.sc <- x[rowSums(x[,2:(2+noregs-1)])==1,]
eps.sc <- colSums(x.sc[,2:(2+noregs-1)])

eps.sc.LOH <- mean(eps.sc[LOHregs])
eps.sc.noLOH <- mean(eps.sc[noLOHregs])

loh.subclonalno.df[nrow(loh.subclonalno.df)+1,] <- c(samp, eps.sc.LOH, eps.sc.noLOH,
                                                     paste0(LOHregs,collapse=','),
                                                     paste0(noLOHregs,collapse=','))

x.noLOH <- x[(rowSums(x[,noLOHregs,drop=F])>0) & (rowSums(x[,LOHregs,drop=F])==0),]
x.LOH <- x[(rowSums(x[,noLOHregs,drop=F])==0) & (rowSums(x[,LOHregs,drop=F])>0),]
t.noLOH <- table(x.noLOH$hla)
t.LOH <- table(x.LOH$hla)

hla.lost <- 'HLA-B*08:01'
hla.kept <- 'HLA-B*07:02'

loh.neoeps.df[nrow(loh.neoeps.df)+1,] <- c(samp, substr(hla.lost, 5,5),'LOH',
                                           t.LOH[hla.lost],t.LOH[hla.kept], sum(t.LOH) )
loh.neoeps.df[nrow(loh.neoeps.df)+1,] <- c(samp, substr(hla.lost, 5,5),'noLOH',
                                           t.noLOH[hla.lost],t.noLOH[hla.kept], sum(t.noLOH) )

loh.neoeps.df$LostHLARatio <- loh.neoeps.df$BoundToLost/loh.neoeps.df$BoundToAll
loh.neoeps.df$ID <- apply(loh.neoeps.df,1, function(x) paste0(x['Sample'],':',x['Allele']))


pLines <- ggplot(loh.neoeps.df, aes(x=Subclone, y=LostHLARatio, group=ID, color=Allele)) +
  geom_line(size=1.2) + geom_point(size=2.5) + scale_x_discrete(expand=c(0.1,0.1)) +
  scale_y_continuous(labels=percent) +
  scale_color_manual(values=c('#7577c9','#51a77f','#e9813d')) +
  labs(y='Percentage of subclonal neoantigens bound to\nlost allele/ total neoantigens', x='Subclone') + 
  theme_bw() + theme(text=element_text(size=16))

loh.bound.df <- loh.neoeps.df[,c('ID', 'BoundToLost', 'BoundToKept', 'Subclone')]
loh.bound.df <- loh.bound.df[order(loh.bound.df$Subclone,decreasing=T),]

pLines2 <- ggplot(melt(loh.bound.df), aes(x=variable, y=value, group=paste0(ID,Subclone), colour=Subclone, shape=ID)) +
  geom_line(size=1.2) + geom_point(size=3) + scale_x_discrete(expand=c(0.1,0.1), labels=c('Lost allele', 'Retained allele')) +
  scale_y_log10() + scale_color_manual(values=c('#4180ae', '#cf5f5f')) +
  scale_shape_manual(values=c(16,17,18,15,4,8)) +
  labs(y='Neoantigens bound to HLA allele', x='Allele') + 
  theme_bw() + theme(text=element_text(size=16)) + guides(shape=F)


p1 <- ggpaired(loh.neoeps.df, x='Subclone', y='LostHLARatio', id='Sample', fill='Subclone',
         line.color = "grey35", line.size = 0.4) + stat_compare_means(paired=T, label.x.npc ='centre') +
  scale_y_continuous(labels=percent) + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  labs(y='Percentage of subclonal neo-epitopes\n lost allele/ kept allele', x='') + guides(fill=F) + theme_bw()
p2 <- ggpaired(loh.neoeps.df, x='Subclone', y='LostTotalRatio', id='Sample', fill='Subclone',
         line.color = "grey35", line.size = 0.4) + stat_compare_means(paired=T, label.x.npc ='centre') +
  scale_y_continuous(labels=percent) + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  labs(y='Percentage of subclonal neo-epitopes\n lost allele/ total neo-epitopes', x='') + guides(fill=F) + theme_bw()
p3 <- ggpaired(subset(loh.neoeps.df, Subclone=='LOH'), cond1='BoundToLost', cond2='BoundToKept', fill='condition',
               line.color = "grey35", line.size = 0.4) + stat_compare_means(paired=T, label.x.npc ='centre') +
  scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  labs(y='Number of neo-epitopes bound to allele\n in subclone with LOH', x='') + guides(fill=F) + theme_bw()
pb <- ggplot(loh.neoeps.df, aes(x=Subclone, y=LostHLARatio, fill=Subclone)) + geom_bar(stat='identity') + facet_wrap(~paste0(Sample,'_',Allele)) +
  scale_y_continuous(labels=percent) + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  labs(y='Percentage of subclonal neo-epitopes\n lost allele/ kept allele', x='') + guides(fill=F) + theme_bw()

p0 <- ggpaired(loh.subclonalno.df, cond1='SubclonalEpsLOH', cond2='SubclonalEpsnoLOH', fill='condition',
               line.color = "grey35", line.size = 0.4) + stat_compare_means(paired=T, label.x.npc ='centre') +
  scale_fill_manual(values=c('skyblue4', 'firebrick3')) + scale_y_log10() + scale_x_discrete(labels=c("SubclonalEpsnoLOH" = "without LOH", "SubclonalEpsLOH" = "with LOH")) +
  labs(y='Number of private neo-epitopes in clone', x='') + guides(fill=F) + theme_bw()


pdf('~/CRCdata/HLA_LOH/CRCmseq/Neoeps_subclonal_LOH.pdf', width=8, height=6)
print(p0);print(pb);print(p1);print(p2);print(p3)
dev.off()

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
epTable$AntigenID <- apply(epTable, 1,function(x) paste0(x['Sample'], ':',x['Identity'], ':',x['peptide'] ) )


recoTable <- read.table('../TCGA_CRC/Neopred_results/PredictedRecognitionPotentials.txt',
                        stringsAsFactors = F, header=T)

recoTable$AntigenID <- apply(recoTable, 1,function(x) paste0(x['Sample'], ':',x['Mutation'], ':',x['MutantPeptide'] ) )

recoTable.imm <- subset(recoTable, NeoantigenRecognitionPotential>1e-2)

epTable.imm <- subset(epTable, AntigenID %in% recoTable.imm$AntigenID)

#epTable <- subset(epTable, Affinity<200)

ploidy.df <- read.table('~/Dropbox/Code/TCGA/CRC_ploidy_master.txt', sep='\t', stringsAsFactors = F, header=T, row.names=1)

allTotVAFs <- data.frame(matrix(vector(), ncol=6)); names(allMutVAFs) <- c('Sample', 'Chr', 'Start', 'VAF', 'CN', 'CCF')

# generate CCFs of mutations

# for(sample in unique(epTable$Sample)){
# #sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
# #exonic <- readExonicFile(sampleFileEx)
# vcf <- read.table(paste0(dir, '/VCF/',sample,'.vcf'), sep='\t', stringsAsFactors = F)
# names(vcf)[c(1,2,11)] <- c('Chr', 'Start', 'Region_0')
#   
# #tmpVAF <- data.frame(Sample=sample, Chr=exonic$Chrom, Start=exonic$Start,  LineID=exonic$LineID, VAF=computeVafAD(exonic, 'Region_0'))
# tmpVAF <- data.frame(Sample=sample, Chr=vcf$Chr, Start=vcf$Start,  VAF=computeVafAD(vcf, 'Region_0'))
# 
# 
# cnafile <- cna.df[cna.df$Patient==sample,'FileName']
# if (length(cnafile)==0){next}
# cna <- read.table(paste0('~/Dropbox/Code/TCGA/CRC_CNA/',cnafile), sep='\t', header=T, stringsAsFactors = F)
# cna$Chromosome <- paste0('chr', cna$Chromosome)
# cna$Segment_Mean <- (((2^(cna$Segment_Mean))*2 -2 )/ploidy.df[sample,'tumorPurity'] + 2)
# #cna$Segment_Mean <- (2^(cna$Segment_Mean))
# 
# tmpVAF$CN <- sapply(1:nrow(tmpVAF), function(i) getCNAofMut(cna,tmpVAF[i,]))
# tmpVAF$CCF <- (tmpVAF$VAF*tmpVAF$CN)*(1/ploidy.df[sample,'tumorPurity'])
# allTotVAFs <- rbind(allTotVAFs, tmpVAF)
# }
# 
# getCNAofMut <- function(cna, mutLine){
#   cnstate <- cna[cna$Chromosome==mutLine$Chr & (cna$Start < mutLine$Start) & (cna$End > mutLine$Start),
#       'Segment_Mean']
#   if (length(cnstate)==0){return(NA)}
#   return(cnstate)
# }

allMutVAFs$mutID <- apply(allMutVAFs, 1,function(x) paste0(x['Sample'], ':',x['LineID'] ) )

epnums <- table(epTable$Sample)
escape.df$EpNumberRP <- epnums[ escape.df$Patient]

sampsToCheck <- subset(escape.df, Escape %in% c('AI', 'NONE') & EpNumberRP > 30  )$Patient

pdf(paste0(dir, '_VAFs_collection_rp.pdf'),width=8,height=5)

fmax=0.8; fmin=0.25
steps <- seq(fmax,fmin,by=(-1e-2))

for (samp in sampsToCheck ){
#for (samp in row.names(subset(summaryTable, (Total>120) & !(MUT) & !(LOH) & !(B2M) & !(PDL)))){
  
#pAll <- ggplot(subset(allTotVAFs, (Sample==samp)), aes(x=CCF)) + geom_histogram(bins=30)
pEp <- ggplot(subset(allTotVAFs, (Sample==samp)), aes(x=CCF, y=..density..)) + geom_histogram(bins=40, fill='grey30', alpha=0.5) +
  geom_histogram(data=subset(allMutVAFs, (Sample==samp) & (mutID %in% epTable$mutID)), aes(x=CCF, y=..density..),bins=25, fill='firebrick3', alpha=0.5) +
  labs(title=samp) +
  geom_histogram(data=subset(allMutVAFs, (Sample==samp) & !(mutID %in% epTable$mutID)), aes(x=CCF,y=..density..), bins=25, fill='skyblue4', alpha=0.5) +
  scale_x_continuous(limits=c(0, 2))

vafAll<-na.omit(subset(allMutVAFs, (Sample==samp))$CCF)
vafTot<-na.omit(subset(allTotVAFs, (Sample==samp))$CCF)
vafEp <- na.omit(subset(allMutVAFs, (Sample==samp) & (mutID %in% epTable$mutID))$CCF)

# cumvaf <- data.frame(invf = (1/steps), cumvaf=sapply(steps, function(x) sum(vafAll>=x))/length(vafAll),
#                      cumvafEp=sapply(steps, function(x) sum(vafEp>=x))/length(vafEp),
#                      cumvafTot=sapply(steps, function(x) sum(vafTot>=x))/length(vafTot))
cumvaf <- data.frame(invf = (1/steps),
                     cumvafEp=sapply(steps, function(x) sum(vafEp>=x)))
cumvaf$cumvafEp <- cumvaf$cumvafEp - min(cumvaf$cumvafEp)
cumvaf$cumvafEp <- cumvaf$cumvafEp/max(cumvaf$cumvafEp)

pCum <- ggplot(cumvaf, aes(x=invf, y=cumvaf)) + geom_line(colour='skyblue4') + geom_line(data=cumvaf, aes(x=invf, y=cumvafEp), colour='firebrick3') +
  geom_line(data=cumvaf, aes(x=invf, y=cumvafTot), colour='grey20') +
  labs(x='1/CCF', y='Cumulative frequency', title=paste0(escape.df[escape.df$Patient==samp, c('EpNumberRP', 'Escape', 'CYT')],collapse=';'))

grid.arrange(pEp, pCum, nrow=1)

}
dev.off()

fmax=0.6; fmin=0.2
steps <- seq(fmax,fmin,by=(-1e-2))

vafEp <- na.omit(subset(allMutVAFs, (Sample %in% goodSamples) & (mutID %in% epTable$mutID))$CCF)
vafTot <- na.omit(subset(allTotVAFs, (Sample %in% goodSamples))$CCF)
vafMut <- na.omit(subset(allMutVAFs, (Sample %in% goodSamples))$CCF)
cumvaf <- data.frame(invf = (1/steps),
                     cumvafEp=sapply(steps, function(x) sum(vafEp>=x)),
                     cumvafTot=sapply(steps, function(x) sum(vafTot>=x)),
                     cumvafEx=sapply(steps, function(x) sum(vafMut>=x)))
cumvaf$cumvafEp <- cumvaf$cumvafEp - min(cumvaf$cumvafEp);cumvaf$cumvafEp <- cumvaf$cumvafEp/max(cumvaf$cumvafEp)
cumvaf$cumvafTot <- cumvaf$cumvafTot - min(cumvaf$cumvafTot);cumvaf$cumvafTot <- cumvaf$cumvafTot/max(cumvaf$cumvafTot)
cumvaf$cumvafEx <- cumvaf$cumvafEx - min(cumvaf$cumvafEx);cumvaf$cumvafEx <- cumvaf$cumvafEx/max(cumvaf$cumvafEx)
cumvaf <- cumvaf[,c('invf','cumvafTot', 'cumvafEx', 'cumvafEp')]

ggplot(cumvaf, aes(x=invf*2, y=cumvafTot)) + geom_line(color='grey35',size=1.2) +
  geom_point(data=cumvaf, aes(x=invf*2, y=cumvafEx), color="#8150a3", size=3, shape=15) +
  geom_line(data=cumvaf, aes(x=invf*2, y=cumvafEp), color="#d63027", size=1.2) +
  guides(color=T)

pvafallcrc <- ggplot(melt(cumvaf, id='invf'), aes(x=invf*2, y=value, color=variable, shape=variable)) +
  geom_point(size=3) + geom_line(size=1.6) +
  scale_color_manual(values=c('grey35',"#8150a3","#d63027"), labels=c('All', 'Exonic', 'Neoantigen')) +
  scale_shape_manual(values=c(26, 15, 26)) +
  labs(color='Mutations', x='Inverse allelic frequency 1/f', y='Cumulative frequency distribution') +
  theme_bw() + theme(text=element_text(size=16)) + guides(shape=F)

# VAF/clonality of neoepitopes --------------------------------------------

# Show that in TCGA neoepitopes are clonal

allTotVAFs <- read.table('~/Dropbox/Code/TCGA/CRC_allVAF_master_file.txt',
                         sep='\t', stringsAsFactors = F, header=T)
allMutVAFs <- read.table('~/Dropbox/Code/TCGA/CRC_exonicVAF_master_file.txt',
                         sep='\t', stringsAsFactors = F, header=T)
allMutVAFs$mutID <- apply(allMutVAFs, 1,function(x) paste0(x['Sample'], ':',x['LineID'] ) )

#filter samples out by purity+ploidy
goodSamples <- gsub('.tumour', '', row.names(subset(ploidy.df, tumorPurity>0.4 & tumorPloidy < 3.6)))

allEpVAFs <- subset(allMutVAFs, mutID %in% epTable.imm$mutID)

maxEpVAFs <- sapply(goodSamples, function(x) max(allEpVAFs[allEpVAFs$Sample==x,'CCF'], na.rm=T)  )
maxEpVAFs[maxEpVAFs==-Inf] <- NA
maxEpVAFs <- maxEpVAFs[order(-maxEpVAFs)]; maxEpVAFs[maxEpVAFs>1] <- 1

mepvaf.df <- na.omit(data.frame(value=maxEpVAFs, ind=1:length(maxEpVAFs),
                        tc = tcavg[match(names(maxEpVAFs), gsub('\\.', '-',names(tcavg)))],
                        cyt = cyts[match(names(maxEpVAFs), gsub('\\.', '-',names(cyts)))]))

orderRows <- order(mepvaf.df[mepvaf.df$value==1, 'cyt'])
mepvaf.df.ord <- rbind(mepvaf.df[orderRows,], mepvaf.df[mepvaf.df$value!=1,])
mepvaf.df.ord$ind <- 1:nrow(mepvaf.df.ord)

pmaxvaf.tcga <- ggplot(mepvaf.df.ord, aes(x=ind, y=value, colour=cyt/2.1293)) +
  geom_point(alpha=0.8, size=3) +
  scale_y_log10(labels=percent)  +
  immColor +
  labs(x='TCGA CRC sample', y='Proportion of cells sharing the most common neoantigen', colour='Normalised\nT-cell avg') +
  theme_bw() + theme(text=element_text(size=16), axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Is there any connection between immune escape and number of total/subclonal neoantigens?


escape.df$Total_neoep <- NA; escape.df$Subclonal_neoep <- NA

for (samp in goodSamples){
  vaftmp <- subset(allMutVAFs, Sample == samp)
  escape.df[escape.df$Patient==samp,'Total_neoep'] <- sum((vaftmp$mutID %in% epTable.imm$mutID), na.rm=T)
  escape.df[escape.df$Patient==samp,'Subclonal_neoep'] <- sum((vaftmp$mutID %in% epTable.imm$mutID) & (vaftmp$CCF < 0.6) & (vaftmp$CCF > 0.3), na.rm=T)
}

x <- subset(escape.df, MSI=='MSI-H')

psc <- ggplot(x, aes(x=(Escape %in% c('NONE','AI')), y=Subclonal_neoep, fill=(Escape %in% c('NONE','AI')))) +
  geom_boxplot() + stat_compare_means(label.x.npc = 'centre') + guides(fill=F) +
  labs(y='High subclonal neantigens (0.3 < CCF < 0.6)', x='Immune escape') + scale_x_discrete(labels=c('Yes (n=42)', 'No (n=7)'))

# VAF vs epitope-strength -------------------------------------------------

samp <- 'TCGA-AD-6889'
d <- epTable.dedup[epTable.dedup$Sample==samp,]
d$VAF <- allMutVAFs[match(d$mutID, allMutVAFs$mutID),'VAF']
ggplot(d, aes(x=Rank, y=VAF)) + geom_point() + scale_y_continuous(limits=c(0, 0.25)) + stat_cor()



# Compare to shuffled control ---------------------------------------------

realSummary <- read.table('~/CRCdata/TCGA_CRC/Neopred_results/Total.neoantigens.summarytable.txt', header=T, stringsAsFactors = F)

x <- realSummary$Sample[25]
randomNeoeps <- c()
for (i in 1:28){
randomSummary <- read.table(paste0('~/CRCdata/TCGA_CRC/Neopred_results/Shuffled_',i,'.neoantigens.summarytable.txt'), header=T, stringsAsFactors = F)
if (x %in% randomSummary$Sample){randomNeoeps <- c(randomNeoeps, randomSummary[randomSummary$Sample==x,'Total'])}
}

compare.df <- data.frame(Random=randomSummary$Total, Real= realSummary[match(randomSummary$Sample, realSummary$Sample),'Total']  )

ggpaired(compare.df, cond1='Random', cond2='Real', fill='condition', line.color='grey70') +
  stat_compare_means(paired=T) + scale_y_continuous(trans='log1p')
  

###########################################################################
# Mutations in normal neoepitopes -----------------------------------------

norm.bed <- read.table(paste0(dir, '/a0201_negative_mutations.bed'), sep='\t', stringsAsFactors = F)
names(norm.bed) <- c('transcript', 'epStart', 'epEnd', 'protLength', 'file', 'transcript2', 'mutStart', 'mutEnd',
                     'mutPos', 'mutAll', 'mutCons', 'impact')

patient <- 'TCGA-AZ-4315'
pat.muts <- subset(norm.bed, file=paste0(patient, '.vep.bed'))

pat.vcf <- read.table(paste0(dir, '/VCF/',patient, '.vcf'), sep='\t', stringsAsFactors = F)
names(pat.vcf) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', 'NORMAL', 'TUMOUR')

pat.vcf$VAF <- computeVafAD(pat.vcf, 'TUMOUR')
pat.vcf$MUTID <- apply(pat.vcf, 1, function(z) paste0(z['CHROM'],':',gsub(' ','',z['POS'])))
pat.vcf.neg <- subset(pat.vcf, MUTID %in% pat.muts$mutPos)

f <- seq(0.75, 0.25, by=-0.01)
vaf.all <- sapply(f, function(x) sum(pat.vcf$VAF>=x))
vaf.neg <- sapply(f, function(x) sum(pat.vcf.neg$VAF>=x))
vaf.df <- data.frame(all=vaf.all, neg=vaf.neg, invf=1/f)
pa <- ggplot(vaf.df, aes(x=invf, y=vaf.all)) + geom_line()
pn <- ggplot(vaf.df, aes(x=invf, y=vaf.neg)) + geom_line()
grid.arrange(pa, pn, nrow=1)
