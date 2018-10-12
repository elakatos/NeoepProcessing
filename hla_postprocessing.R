#############################################################################
# HLA types ---------------------------------------------------------------
# Gather HLA types from file list
# Convert to format recognized by netMHCpan and substitute with nearest neighbour to
# obtain the type netMHCpan would use in its run (e.g. HLA-02:04 == HLA-02:01)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

hlaConvert <- function(x){
  if (is.na(x)){
    y <- NA
  }
  else{
    y <- paste0('HLA-', toupper(substr(x,5,5)),substr(x,7,8),':',substr(x,10,11))
  }
  return(y)
}

nnGrep <- function(x){
  i <- gregexpr("HLA", x) #find the two HLAs mentioned
  origHLA <- substr(x, i[[1]][1], i[[1]][1]+9)
  nn <- substr(x, i[[1]][2], i[[1]][2]+9)
  j <- regexpr("[0-9]\\.[0-9][0-9][0-9]", x) #match the format of distance
  dist <- substr(x, j, j+4)
  return(list(HLA=origHLA, NN=nn, Distance=dist))
}

#Read in, filter out possibly incorrect predictions and convert to netMHCpan format
#dir = 'TCGA_CRC'
#hlas <- read.table(paste0('~/RNAseq/Neoepitopes/',dir,'/hlatypes.txt'), sep='\t', header=T, row.names=1, stringsAsFactors = F)
hlaFile <- '~/Dropbox/Code/TCGA/hlatypes_total_CRC.txt'

hlasOrig <- read.table(hlaFile, sep='\t', header=T, row.names=1, stringsAsFactors = F)
hlasCorrect <- subset(hlasOrig, !( (HLA.A_1=='hla_a_01_01_01_01') & (is.na(HLA.A_2)) & (HLA.B_1== 'hla_b_07_02_01') & (is.na(HLA.B_2)) & (HLA.C_1=='hla_c_01_02_01') & (is.na(HLA.C_2))  ))
hlas <- as.data.frame(apply(hlasCorrect, 2,function(r) sapply(r, function(x) hlaConvert(x))) )
hlaList <- as.vector(as.matrix(hlas))

#Build nearest neighbour table from encountered HLA types
hlaNN <- data.frame(matrix(vector(), ncol=3))
names(hlaNN) <- c('HLA', 'NN', 'Distance')
hlaMaps <- scan(file='hla_mappings.txt', what=character(), sep='\n') #Collection of lines stating nearest neighbours
for (i in 1:length(hlaMaps)){ hlaNN[i,] <- nnGrep(hlaMaps[i]) }
hlaNN <- hlaNN[!duplicated(hlaNN$HLA),]

hlaList.mapped <- mapvalues(hlaList, from=hlaNN$HLA, to=hlaNN$NN)
hlas.mapped <- as.data.frame(sapply(hlas, function(x) mapvalues(x, from=hlaNN$HLA, to=hlaNN$NN)))

#todo <- setdiff(hlaList, hlaNN$HLA)
#todo <- sapply(todo, function(x) paste0(substr(x, 1, 7),substr(x, 9, 10)))

# Generate HLA file for a particular type ---------------------------------

allele<-'HLA-C12:03'
nonAllele <- hlas.mapped!=allele
hlasOut <- hlasCorrect; hlasOut[nonAllele] <- NA #NA all other entries in the hla table
hlasOut <- hlasOut[rowSums(is.na(hlasOut))<6,]
#Further filtering: exclude known MSI and hyper-mutated ones
#nonHyper <- subset(clin.df, (MSI %in% c('MSS', NA)) &(Hypermut %in% c(0, NA)) )
#hlasOut <- hlasOut[row.names(hlasOut) %in% nonHyper$Patient,]
outDir <- '~/Dropbox/Code/TCGA/'
write.table(hlasOut, paste0(outDir,'hlatypes_',sub(':','', allele),'.txt'), sep='\t', quote=F)


# Shuffle HLA haplotypes --------------------------------------------------
for (i in 31:80){
hlasShuffled$HLA.A_1 <- sample(hlasShuffled$HLA.A_1)
hlasShuffled$HLA.A_2 <- sample(hlasShuffled$HLA.A_2)
hlasShuffled$HLA.B_1 <- sample(hlasShuffled$HLA.B_1)
hlasShuffled$HLA.B_2 <- sample(hlasShuffled$HLA.B_2)
hlasShuffled$HLA.C_1 <- sample(hlasShuffled$HLA.C_1)
hlasShuffled$HLA.C_2 <- sample(hlasShuffled$HLA.C_2)
hlasOut <- hlasShuffled[sample(1:nrow(hlasShuffled),50),]
write.table(hlasOut, file=paste0('~/Dropbox/Code/TCGA/hlatypes_total_shuffled_',i,'.txt'), sep='\t', quote=F)
}
#####################################################################################
# LOHHLA analysis ---------------------------------------------------------

# Create a label for filtering depending on the reliability and output of CN measurement
labelHLAAI <- function(line){
  if (is.na(line$HLA_type1copyNum_withBAFBin_upper)){return(NaN)}
  lab <- 'none'
  if (line$PVal_unique<0.01){lab <- 'AI'}
  if ((line$HLA_type1copyNum_withBAFBin_lower>1.5) & (line$HLA_type2copyNum_withBAFBin_lower > 1.5 ))
  { lab <- 'gain'}
  if ( (lab=='AI') & (line$HLA_type1copyNum_withBAFBin_upper< 0.7) & (line$HLA_type1copyNum_withBAFBin < 0.5 ))
  { lab <- 'LOH'}
  if ((lab=='AI') & (line$HLA_type2copyNum_withBAFBin < 0.5) & (line$HLA_type2copyNum_withBAFBin_upper < 0.7 ))
  { lab <- 'LOH'}
  #if ((lab %in% c('AI','LOH')) & (line$HLA_type1copyNum_withBAFBin_upper<0.5) & (line$HLA_type2copyNum_withBAFBin_upper <0.5 ))
  #{ lab <- 'loss'}
  return(lab)
}

labelHLArel <- function(line, mincov){
  if (is.na(line$HLA_type1copyNum_withBAFBin_upper)){return(F)}
  lab <- T
  confint1 <- line$HLA_type1copyNum_withBAFBin_upper - line$HLA_type1copyNum_withBAFBin_lower
  confint2 <- line$HLA_type2copyNum_withBAFBin_upper - line$HLA_type2copyNum_withBAFBin_lower
  if ( ((confint1>1.5) & (line$HLA_type1copyNum_withBAFBin<1.5) & (line$HLA_type1copyNum_withBAFBin>0)) |  ((confint2>1.5) & (line$HLA_type2copyNum_withBAFBin<1.5) & (line$HLA_type2copyNum_withBAFBin>0)) )
  { lab <- F}
  if ( (confint1>6) | (confint2>6))
  {lab <- F}
  if ( ((abs(line$HLA_type1copyNum_withBAFBin)<0.4) & (line$HLA_type1copyNum_withBAFBin_upper>0.8)) |   ((abs(line$HLA_type2copyNum_withBAFBin)<0.4) & (line$HLA_type2copyNum_withBAFBin_upper>0.8)))
  {lab <- F}
  if ( line$numMisMatchSitesCov < mincov )
  {lab <- F}
  return(lab)
}

analyseLohhla <- function(lohhla.master, clin.df){
  lohhla.signif <- lohhla.master[!is.na(lohhla.master$PVal_unique),]
  lohhla.signif <- lohhla.signif[lohhla.signif$PVal_unique<0.01,]
  
  lohhla.patients <- data.frame(row.names = unique(lohhla.master$region))
  lohhla.patients$ID <- sapply(row.names(lohhla.patients), function(x) substr(x,1,12))
  lohhla.patients$AI <- sapply(row.names(lohhla.patients), function(x) x %in% lohhla.signif$region)
  lohhla.patients$CN <- sapply(row.names(lohhla.patients),
                               function(x) lohhla.master[lohhla.master$region==x, 'Label'][1]!='NaN')
  lohhla.patients$LOH <- sapply(row.names(lohhla.patients),
                                function(x) ('LOH' %in% lohhla.master[lohhla.master$region==x, 'Label']))
  lohhla.patients$LOSS <- sapply(row.names(lohhla.patients),
                                 function(x) ('loss' %in% lohhla.master[lohhla.master$region==x, 'Label']))
  lohhla.patients$NORM <- sapply(row.names(lohhla.patients),
                                 function(x) sum(lohhla.master[lohhla.master$region==x, 'Label']!='none')==0)
  lohhla.patients$HIGH <- sapply(row.names(lohhla.patients),
                                 function(x) ('gain' %in% lohhla.master[lohhla.master$region==x, 'Label']))
  
  lohhla.patients$MSI <- clin.df[match(lohhla.patients$ID, clin.df$Patient), 'MSI']=='MSI-H'
  lohhla.patients$HYP <- clin.df[match(lohhla.patients$ID, clin.df$Patient), 'Hypermut']==1
  return(list(sig=lohhla.signif, pat=lohhla.patients))
}


dir <- '~/CRCdata/HLA_LOH/IBD/Rand_CN/'
fileList <- list.files(dir, pattern='*.10.DNA.HLAloss*')

lohhla.master <- data.frame(matrix(vector()))
for (sampleName in fileList){
  hlapred.failed <- NULL
  hlapred <- tryCatch(
    {read.table(paste0(dir,sampleName), sep='\t', header=T,
                stringsAsFactors = F)},
    error=function(e){
      message(paste0('File for ',substr(sampleName,1,12),
                     ' does not have content, sample is probably completely homozygous.'))
      return(NA)
    }
  )
  if (is.na(hlapred)){next}
  hlapred <- hlapred[!duplicated(hlapred[,c('region', 'HLA_A_type1')]),]
  lohhla.master <- rbind(lohhla.master, hlapred)
}


lohhla.master$QVal <- p.adjust(lohhla.master$PVal_unique, method='fdr')
lohhla.master$Label <- sapply(1:nrow(lohhla.master), function(i) labelHLAAI(lohhla.master[i,]))
lohhla.master$Rel <- sapply(1:nrow(lohhla.master), function(i) labelHLArel(lohhla.master[i,], 5))

lohhla.master <- read.table('~/Dropbox/Code/TCGA/CRC_LOHHLA_master_file.txt', stringsAsFactors = F, header=T)
# Filtered version: filter on the master table level

lohhla.master <- subset(lohhla.master, Rel==T)
anal <- analyseLohhla(lohhla.master, clin.df)
lohhla.patients <- anal$pat

lohhla.patients$B2M <- clin.df[match(lohhla.patients$ID, clin.df$Patient), 'B2M']>0
lohhla.patients$MUT <- sapply(lohhla.patients$ID, function(x) x %in% muthla.signif$individual)

#lohhla.patients <- subset(lohhla.patients, !(HYP & !MSI)) #Filter out hypermutated MSS cases (e.g. POLE)
lohhla.patients.sub <- na.omit(subset(lohhla.patients, CN==T))

t1 <- fisher.test(table(lohhla.patients$MSI, lohhla.patients$AI))
p1 <- ggplot(na.omit(lohhla.patients), aes(x=MSI, fill=AI)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t1$p.value,5)), x='MSI high', y='')
t2 <- fisher.test(table(lohhla.patients.sub$MSI, lohhla.patients.sub$LOH))
p2 <- ggplot(na.omit(lohhla.patients.sub), aes(x=MSI, fill=LOH)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t2$p.value,5)), x='MSI high', y='')
t3 <- fisher.test(table(lohhla.patients.sub$MSI, lohhla.patients.sub$LOSS))
p3 <- ggplot(na.omit(lohhla.patients.sub), aes(x=MSI, fill=LOSS)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t3$p.value,5)), x='MSI high', y='')


#pdf('~/Dropbox/Code/TCGA/Figures/MSI_comp_all_Filtered.pdf', width=12, height=5)
#grid.arrange(p1, p2, p3, nrow=1)
#dev.off()


t4 <- fisher.test(table(lohhla.patients$MSI, lohhla.patients$B2M))
p4 <- ggplot(na.omit(lohhla.patients), aes(x=MSI, fill=B2M)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t4$p.value,5)), x='MSI high', y='')

t5 <- fisher.test(table(lohhla.patients$MSI, lohhla.patients$MUT))
p5 <- ggplot(na.omit(lohhla.patients), aes(x=MSI, fill=MUT)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t4$p.value,5)), x='MSI high', y='')

pdf('~/Dropbox/Code/TCGA/Figures/MSI_comp_filtered_HLA.pdf', width=15, height=5)
grid.arrange(p1, p2, p4, p5, nrow=1)
dev.off()

t4 <- fisher.test(table(lohhla.patients.sub$B2M, lohhla.patients.sub$MUT))
p4 <- ggplot(na.omit(lohhla.patients), aes(x=B2M, fill=MUT)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t4$p.value,5)), x='B2M', y='')
t3 <- fisher.test(table(lohhla.patients.sub$LOH, lohhla.patients.sub$MUT))
p3 <- ggplot(na.omit(lohhla.patients), aes(x=LOH, fill=MUT)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t3$p.value,5)), x='LOH', y='')
t2 <- fisher.test(table(lohhla.patients.sub$LOH, lohhla.patients.sub$B2M))
p2 <- ggplot(na.omit(lohhla.patients), aes(x=LOH, fill=B2M)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t2$p.value,5)), x='LOH', y='')

pdf('~/Dropbox/Code/TCGA/Figures/MSI_comp_filtered_mutations.pdf', width=11, height=5)
grid.arrange(p2, p3, p4, nrow=1)
dev.off()


###########################################################################
# Alternative escapes -----------------------------------------------------

# HLA mutations -----------------------------------------------------------

#dirList <- paste0('~/CRCdata/HLA_MUT/CRCmseq/Set',1:10,'_unk')
dirList <- '~/CRCdata/HLA_MUT/TCGA'

muthla.master <- data.frame(matrix(vector()))

correctHeader <- c('individual', 'contig', 'position', 'context', 'ref_allele', 'alt_allele',
                   'tumor_name', 'normal_name', 'score', 'dbsnp_site', 'covered', 'power', 'tumor_power',
                   'normal_power', 'normal_power_nsp', 'normal_power_wsp', 'total_reads', 'map_Q0_reads',
                   'init_t_lod', 't_lod_fstar', 't_lod_fstar_forward', 't_lod_fstar_reverse', 'tumor_f',
                   'contaminant_fraction', 'contaminant_lod', 't_q20_count', 't_ref_count', 't_alt_count',
                   't_ref_sum', 't_alt_sum', 't_ref_max_mapq', 't_alt_max_mapq', 't_ins_count', 't_del_count',
                   'normal_best_gt', 'init_n_lod', 'normal_f', 'n_q20_count', 'n_ref_count', 'n_alt_count',
                   'n_ref_sum', 'n_alt_sum', 'power_to_detect_positive_strand_artifact',
                   'power_to_detect_negative_strand_artifact', 'strand_bias_counts', 'tumor_alt_fpir_median',
                   'tumor_alt_fpir_mad', 'tumor_alt_rpir_median', 'tumor_alt_rpir_mad', 'observed_in_normals_count',
                   'failure_reasons', 'judgement','exon_or_intron', 'exon_number', 'intron_number', 'protein_change')

muthla.master <- data.frame(matrix(vector()))

noMutPreds <- c()
for (dir in dirList){
  sampleList <- list.files(dir, pattern='*')
  for (sample in sampleList){
    hlaList <- list.files(paste0(dir, '/', sample), pattern='*.mutect.unfiltered.annotated')
    if (length(hlaList)==0){noMutPreds = c(noMutPreds, sample)}
    for (hla in hlaList){
      muthla <- tryCatch(
        {read.table(paste0(dir,'/', sample,'/', hla), stringsAsFactors = F, header=F, skip=1, sep='\t')},
        error=function(e){return(NA)}
      )
      if (!is.na(muthla)){
        names(muthla) <- correctHeader
        muthla.master <- rbind(muthla.master, muthla)
      }
    }
    hlaList2 <- list.files(paste0(dir, '/', sample), pattern='*.strelka_indels.unfiltered.annotated')
    for (hla in hlaList2){
    indels.tmp <- read.table(paste0(dir, '/', sample, '/', hla), stringsAsFactors = F, header=T, sep='\t')
    if(nrow(indels.tmp)>0){print('Woohoo')}
    }
  }
}

muthla.master <- read.table('~/Dropbox/Code/TCGA/CRC_MUTHLA_master_file.txt', stringsAsFactors = F)

muthla.signif <- subset(muthla.master, protein_change!=-1)
muthla.signif$nonsyn <- sapply(muthla.signif$protein_change, function(x) substr(x,3,3) != substr(x, nchar(x), nchar(x)))
muthla.signif <- subset(muthla.signif, nonsyn)


# Heterogeneity of loci ---------------------------------------------------

numAlleles <- table(lohhla.master$region)
pnAlleles <- ggplot(melt(table(numAlleles)),aes(x=numAlleles, y=value, fill=numAlleles)) + geom_bar(stat='identity') +
  labs(y='Number of patients', x='Number of heterozygous HLA loci')


divAlleles <- sapply(1:length(numAlleles),
                     function(x) (length(unique(lohhla.master[lohhla.master$region==(names(numAlleles)[x]),'Label']))-1)/numAlleles[x] )
divAlleles[divAlleles>1/3] <- 1; divAlleles[divAlleles==1/3]<-0.5
pdAlleles <- ggplot(melt(table(divAlleles)),aes(x=divAlleles, y=value, fill=divAlleles)) + geom_bar(stat='identity') +
  labs(y='Number of patients', x='Diversity of copy number of ABC loci')

divAlleles.sub <- divAlleles[row.names(lohhla.patients[lohhla.patients$LOH,])]
pdAllelesLOH <- ggplot(melt(table(divAlleles.sub)),aes(x=divAlleles.sub, y=value, fill=divAlleles.sub)) + geom_bar(stat='identity') +
  labs(y='Number of patients', x='Diversity of copy number of ABC loci in patients with LOH')

pdf('~/Dropbox/Code/TCGA/Figures/HLA_div_LOH.pdf', width=7, height=5)
print(pnAlleles)
print(pdAlleles)
print(pdAllelesLOH)
dev.off()

# Anywhere having AI on A & B but not C?
numAlleles.3 <- names(numAlleles[numAlleles==3])

alleleDF <- data.frame(row.names=numAlleles.3)

for (regionName in numAlleles.3){
  x <- lohhla.master[lohhla.master$region==regionName,c('KeptAllele','PVal_unique')]
  if (sum(is.na(x$KeptAllele))>0){next}
  alleleDF[regionName,'pA'] <- subset(x, startsWith(KeptAllele, 'hla_a'))$PVal_unique < 0.01
  alleleDF[regionName,'pB'] <- subset(x, startsWith(KeptAllele, 'hla_b'))$PVal_unique < 0.01
  alleleDF[regionName,'pC'] <- subset(x, startsWith(KeptAllele, 'hla_c'))$PVal_unique < 0.01
}
alleleDF <- na.omit(alleleDF)

dim(subset(alleleDF, pA & pB & pC))
dim(subset(alleleDF, !pA & !pB & pC))



# HLA ggplot heatmap ---------------------------------------------------------

# Get information for all alleles : generate melted format already
lohhla.df <- data.frame('Region'=character(), 'Allele' = character(), 'CopyNumber'=numeric(), 'p.value'=numeric(),stringsAsFactors = F)

for (i in 1:nrow(lohhla.master)){
  x <- lohhla.master[i,]
#r <- substr(x$region,1,nchar(x$region)-6)
  r <- x$region
allele <- paste0('HLA-',toupper(substr(x$HLA_A_type1,5,5)))
p <- x$PVal_unique
val <- ifelse(x$PVal_unique<0.01, min(x$HLA_type1copyNum_withBAFBin, x$HLA_type2copyNum_withBAFBin), NA)
lohhla.df[i,] <- c(r, allele, val-1.5, p)
}
lohhla.df$CopyNumber <- as.numeric(lohhla.df$CopyNumber)
lohhla.df$p.value <- as.numeric(lohhla.df$p.value)
lohhla.df <- lohhla.df[order(lohhla.df$Region),]

pLOH <- ggplot(lohhla.df, aes(y=Region, x=Allele, fill=CopyNumber)) + geom_tile() +
  scale_fill_gradientn(colours=c('red4', 'red4', 'red4', 'red4', 'brown3', 'mistyrose2', 'bisque', 'cornsilk3'), limits=c(-2.5, 2), na.value='grey70') +
  theme_bw()

# HLA modifications heatmap -----------------------------------------------

muthla.df <- data.frame('Region'=character(), 'Allele' = character(), 'HLA mutation'=numeric())

for (i in 1:nrow(head(muthla.signif))){
  x <- muthla.signif[i,]
  #r <- substr(x$region,1,nchar(x$region)-6)
  r <- x$individual
  allele <- paste0('HLA-',toupper(substr(x$contig,5,5)))
  val <- ifelse(x$nonsyn, 1, 0)
  muthla.df[i,] <- c(r, allele, val)
}



# Escape table ------------------------------------------------------------

lohhla.master$region <- gsub('.tumour', '', lohhla.master$region)
lohhla.ai <- subset(lohhla.master, PVal_unique<0.01)
lohhla.loh <- subset(lohhla.ai, Label=='LOH')

for (pat in unique(lohhla.ai$region)){ escape.df[escape.df$Patient==pat, 'AI'] <- paste0(lohhla.ai[lohhla.ai$region==pat,'LossAllele'], collapse=',') }

for (pat in unique(lohhla.loh$region)){ escape.df[escape.df$Patient==pat, 'HLA_LOH'] <- paste0(lohhla.loh[lohhla.loh$region==pat,'LossAllele'], collapse=',') }

for (pat in unique(muthla.signif$individual)){ escape.df[escape.df$Patient==pat, 'HLA_MUT'] <- paste0(muthla.signif[muthla.signif$individual==pat,'contig'], collapse=',') }

escape.df$B2M_MUT <- clin.df[match(escape.df$Patient, clin.df$Patient), 'B2M']
escape.df$PDL1 <- exprpdl[match(gsub('-', '.',escape.df$Patient), names(exprpdl))] > mean(exprpdl.norm)+2*sd(exprpdl.norm)
escape.df$CTLA4 <- exprctla[match(gsub('-', '.',escape.df$Patient), names(exprctla))] > mean(exprctla.norm)+2*sd(exprctla.norm)
escape.df$CYT <- cyts[match(gsub('-', '.',escape.df$Patient), names(cyts))]
escape.df$MSI <- clin.df[match(escape.df$Patient, clin.df$Patient), 'MSI']

escape.df <- read.table('~/Dropbox/Code/TCGA/CRC_escape_master_file.txt',stringsAsFactors = F, sep='\t', header=T)
escape.df <- escape.df[!is.na(escape.df$PDL1),]
escape.df <- subset(escape.df, FULL_INFO)

escape.df$Escape <- NA
escape.df[ (is.na(escape.df$HLA_LOH) & !is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT==0 & !escape.df$PDL1 & !escape.df$CTLA4), 'Escape' ] <- 'HLA_MUT'
escape.df[ (is.na(escape.df$HLA_LOH) & is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT==0 & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT'
escape.df[ (!is.na(escape.df$HLA_LOH) & is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT==0 & !(escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'HLA_LOH'
escape.df[ (is.na(escape.df$HLA_LOH) & is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT>0 & !escape.df$PDL1 & !escape.df$CTLA4), 'Escape' ] <- 'B2M_MUT'
escape.df[ ((is.na(escape.df$HLA_LOH)) & is.na(escape.df$HLA_MUT) & (escape.df$B2M_MUT>0) & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT_&_B2M'
escape.df[ (!(is.na(escape.df$HLA_LOH)) & is.na(escape.df$HLA_MUT) & (escape.df$B2M_MUT==0) & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT_&_LOH'
escape.df[ ((is.na(escape.df$HLA_LOH)) & !is.na(escape.df$HLA_MUT) & (escape.df$B2M_MUT==0) & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT_&_HLAMUT'
escape.df[ (!(is.na(escape.df$HLA_LOH)) & (!is.na(escape.df$HLA_MUT) | (escape.df$B2M_MUT>0)) & !(escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'LOH_&_MUT'
escape.df[ ((!is.na(escape.df$HLA_LOH)) + (!is.na(escape.df$HLA_MUT)) + (escape.df$B2M_MUT>0) + (escape.df$PDL1 | escape.df$CTLA4))>2, 'Escape' ] <- 'COMBINATION'
escape.df[is.na(escape.df$Escape) & !(is.na(escape.df$AI)),'Escape'] <- 'AI'
escape.df[is.na(escape.df$Escape),'Escape'] <- 'NONE'

escape.df.msi <- subset(escape.df, MSI=='MSI-H')
escape.df.mss <- subset(escape.df, MSI=='MSI-L')
escape.df.pole <- subset(escape.df, MSI=='POLE')

pesc1 <- ggplot(data.frame(value=as.numeric(table(escape.df.mss$Escape)),var=names(table(escape.df.mss$Escape))), aes(x='', y=value ,fill=var)) +
  geom_bar(stat='identity') +
  coord_polar('y', start=0) +
  scale_fill_manual(values=c('#e1b0b0',
                             '#7577c9',
                             '#e2b01d', '#7abf9f', '#85bb59','#e9813d',
                             '#a77955',
                             '#d53e4f', '#3288bd',
                             'grey80')) +
  labs(fill='Escape mechanism', x='', y='', title=paste0('MSS tumours (n=',nrow(escape.df.mss),')')) + 
  theme_minimal() + theme(axis.text.x=element_blank(), panel.grid=element_blank(), text=element_text(size=16))

pesc2 <- ggplot(data.frame(value=as.numeric(table(escape.df.msi$Escape)),var=names(table(escape.df.msi$Escape))), aes(x='', y=value ,fill=var)) +
  geom_bar(stat='identity') +
  coord_polar('y', start=0) +
  scale_fill_manual(values=c('#7577c9',
                             '#e2b01d','#7abf9f', '#85bb59','#e9813d',
                             '#a77955',
                             '#d53e4f', '#3288bd',# '#aa4c9a',
                             'grey80')) +
  labs(fill='Escape mechanism', x='', y='', title=paste0('MSI tumours (n=',nrow(escape.df.msi),')')) + 
  theme_minimal() + theme(axis.text.x=element_blank(), panel.grid=element_blank(), text=element_text(size=16))


pesc3 <- ggplot(data.frame(value=as.numeric(table(escape.df.pole$Escape)),var=names(table(escape.df.pole$Escape))), aes(x='', y=value ,fill=var)) +
  geom_bar(stat='identity') +
  coord_polar('y', start=0) +
  scale_fill_manual(values=c('#e2b01d','#7abf9f', '#85bb59',
                             '#a77955',
                             #'#d53e4f',
                             '#aa4c9a',
                             'grey80')) +
  labs(fill='Escape mechanism', x='', y='', title=paste0('POLE tumours (n=',nrow(escape.df.pole),')')) + 
  theme_minimal() + theme(axis.text.x=element_blank(), panel.grid=element_blank(), text=element_text(size=16))


######################################################################################
# Gene expression ---------------------------------------------------------

tcga.tpm <- read.table('~/Dropbox/Code/TCGA/CRC_RNA_expression.tpm', stringsAsFactors = F, header=T)
tcga.normal.tpm <- read.table('~/Dropbox/Code/TCGA/Normal_colon_samples.tpm', stringsAsFactors = F, header=T, row.names=1)

# Cytolytic activity ------------------------------------------------------

gzma <- 'ENSG00000145649'
prf <- 'ENSG00000180644'
pdl <- 'ENSG00000120217'
ctla <- 'ENSG00000163599'

tcag <- c('ENSG00000108691', 'ENSG00000277632', 'ENSG00000275302', 'ENSG00000138755', 'ENSG00000169245',
          'ENSG00000153563', 'ENSG00000241106', 'ENSG00000242574', 'ENSG00000204252', 'ENSG00000113088',
          'ENSG00000163600', 'ENSG00000125347')
exprcyt <- log10(tcga.tpm[c(gzma, prf),]+1)
cyts <- apply(exprcyt, 2, function(x) sqrt(x[1]*x[2]))
ggplot(data.frame(value=cyts), aes(x=value)) + geom_density()

exprpdl <- sapply(tcga.tpm[pdl,], function(x) log10(x+1))

exprpdl.norm <- sapply(tcga.normal.tpm[pdl,], function(x) log10(x+1))
exprpdl.norm <- exprpdl.norm[exprpdl.norm<0.9] #get rid of clear outlier
ggplot(data.frame(value=exprpdl), aes(x=value)) + geom_density(colour='firebrick3') +
  geom_density(data=data.frame(value=exprpdl.norm), aes(x=value)) +
  geom_vline(xintercept=(mean(exprpdl.norm)+2*sd(exprpdl.norm)))

pdl.over <- names(exprpdl)[exprpdl > (mean(exprpdl.norm)+2*sd(exprpdl.norm))]

exprctla <- sapply(tcga.tpm[ctla,], function(x) log10(x+1))
exprctla.norm <- sapply(tcga.normal.tpm[ctla,], function(x) log10(x+1))
ggplot(data.frame(value=exprctla), aes(x=value)) + geom_density(colour='firebrick3') +
  geom_density(data=data.frame(value=exprctla.norm), aes(x=value)) +
  geom_vline(xintercept=(mean(exprctla.norm)+2*sd(exprctla.norm)))

ctla.over <- names(exprctla)[exprctla > (mean(exprctla.norm)+2*sd(exprctla.norm))]

#CTLA vs PDL
ggplot(data.frame(c=exprctla,p=exprpdl), aes(x=p, y=c)) + geom_point() +
  geom_hline(yintercept=(mean(exprctla.norm)+2*sd(exprctla.norm))) +
  geom_vline(xintercept=(mean(exprpdl.norm)+2*sd(exprpdl.norm)))

#All and CYT
totexpr <- data.frame(p=exprpdl, c=cyts, c4=exprctla,
                      MSI=clin.df[match(gsub('\\.', '-', names(exprpdl)), clin.df$Patient),'MSI'])
ggplot(totexpr[totexpr$MSI %in% c('MSI-H', 'MSI-L'),], aes(x=c, y=p, colour=MSI)) +
  geom_point(size=2) + stat_cor(mapping = aes(x=c,y=p,colour='')) +
  scale_color_manual(values=c('black', 'darkorange3', 'darkseagreen4')) + theme_bw() +
  guides(colour=F) + labs(x='CYT score (logTPM)', y='PD-L1 expression (logTPM)') +
  theme(text=element_text(size=16))


#T-cell associated genes
exprtcag <- log10(tcga.tpm[tcag,]+1)
tcavg <- apply(exprtcag, 2, mean)

lohhla.patients$CYT <- cyts[match(gsub('-','.',lohhla.patients$ID), names(cyts))]
ggplot(na.omit(lohhla.patients), aes(x=MSI, y=CYT, fill=MSI)) + geom_violin() +
  scale_y_continuous(trans='log10') + stat_compare_means()

lohhla.patients$PDL <- exprpdl[match(gsub('-','.',lohhla.patients$ID), names(exprpdl))]
ggplot(na.omit(lohhla.patients), aes(x=MSI, y=PDL, fill=MSI)) + geom_violin() +
  scale_y_continuous(trans='log10') + stat_compare_means()
lohhla.patients$CTLA <- exprpdl[match(gsub('-','.',lohhla.patients$ID), names(exprpdl))]
ggplot(na.omit(lohhla.patients), aes(x=MSI, y=CTLA, fill=MSI)) + geom_violin() +
  scale_y_continuous(trans='log10') + stat_compare_means()
lohhla.patients.mss <- subset(lohhla.patients, HYP==F)
ggplot(na.omit(lohhla.patients.mss), aes(x=LOH, y=PDL, fill=LOH)) + geom_violin() +
  scale_y_continuous(trans='log10') + stat_compare_means()

lohhla.patients.mss <- subset(lohhla.patients, HYP==F)
lohhla.patients.mss.norm <- subset(lohhla.patients, (HYP==F) & (HIGH==T | LOH==T)   )
ggplot(na.omit(lohhla.patients.mss), aes(x=LOH, y=CYT, fill=LOH)) + geom_violin() +
  scale_y_continuous(trans='log10') + stat_compare_means()



# Purity in CRC from Pierre -----------------------------------------------

dir <- '~/Desktop/tmp/'
fileList <- list.files(dir, pattern='*.stats-cna-baf*')

mseqPur <- data.frame(matrix(vector(), ncol=4))
for(sample in fileList){
tp <- read.table(paste0(dir,sample), sep=',', stringsAsFactors = F)
mseqPur <- rbind(mseqPur, tp[,c(3,6,4,5)])
}
mseqPur$V4 <- sapply(mseqPur$V4, function(x) substr(strsplit(x, ',')[[1]][1], 2, nchar(strsplit(x, ',')[[1]][1])))
row.names(mseqPur) <- mseqPur$V3; mseqPur <- mseqPur[,c(2,3,4)]
names(mseqPur) <- c('Ploidy', 'tumorPurity', 'tumorPloidy')
write.table(mseqPur, file='~/Desktop/tmp/solutions2.txt', sep='\t', quote=F)

###########################################################################
# Individual cases --------------------------------------------------------

case <- 'TCGA-CK-4951'
lohhla.master[lohhla.master$region==paste0(case,'.tumour'),]
clin.df[clin.df$Patient==case,]
if(case %in% muthla.master$individual){print(muthla.master[muthla.master$individual==case,])}


###########################################################################
# Paper figures -----------------------------------------------------------

pcrcexpr <- ggplot(totexpr[totexpr$MSI %in% c('MSI-H', 'MSI-L'),], aes(x=c, y=p, colour=MSI)) +
  geom_point(size=2) + stat_cor(mapping = aes(x=c,y=p,colour='')) +
  scale_color_manual(values=c('black', 'darkorange3', 'darkseagreen4')) + theme_bw() +
  guides(colour=F) + labs(x='CYT score (logTPM)', y='PD-L1 expression (logTPM)') +
  theme(text=element_text(size=16))
