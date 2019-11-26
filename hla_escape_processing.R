library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
source('rfunctions_postprocessing.R')


#############################################################################
# HLA types ---------------------------------------------------------------
#############################################################################
# Gather HLA types from file list
# Convert to format recognized by netMHCpan and substitute with nearest neighbour to
# obtain the type netMHCpan would use in its run (e.g. HLA-02:04 == HLA-02:01)
# Subset particular HLA types or shuffle them

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
  origHLA <- substr(x, i[[1]][1], i[[1]][1]+10)
  origHLA <- gsub(' ', '', origHLA)
  nn <- substr(x, i[[1]][2], i[[1]][2]+9)
  j <- regexpr("[0-9]\\.[0-9][0-9][0-9]", x) #match the format of distance
  dist <- substr(x, j, j+4)
  return(list(HLA=origHLA, NN=nn, Distance=dist))
}

#Read in, filter out possibly incorrect predictions and convert to netMHCpan format
hlaFile <- '~/CRCdata/TCGA_UCEC/hlatypes.txt'

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


# Generate HLA file for a particular type ---------------------------------

allele<-'HLA-C12:03'
nonAllele <- hlas.mapped!=allele
hlasOut <- hlasCorrect; hlasOut[nonAllele] <- NA #NA all other entries in the hla table
hlasOut <- hlasOut[rowSums(is.na(hlasOut))<6,]
outDir <- '~/TCGA/'
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
write.table(hlasOut, file=paste0('~/TCGA/hlatypes_total_shuffled_',i,'.txt'), sep='\t', quote=F)
}


#############################################################################
# LOHHLA analysis ---------------------------------------------------------
#############################################################################

dir <- '~/CRCdata/HLA_LOH/TCGA_UCEC/'
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

#adjust p-values and evaluate results
lohhla.master$QVal <- p.adjust(lohhla.master$PVal_unique, method='fdr')
lohhla.master$Label <- sapply(1:nrow(lohhla.master), function(i) labelHLAAI(lohhla.master[i,]))
lohhla.master$Rel <- sapply(1:nrow(lohhla.master), function(i) labelHLArel(lohhla.master[i,], 5))
lohhla.master$Label[lohhla.master$Label=='LOH' & !lohhla.master$Rel] <- 'AI'

# Filter out unreliable ones
#lohhla.master <- subset(lohhla.master, Rel==T)

lohhla.signif <- lohhla.master[!is.na(lohhla.master$PVal_unique),]
lohhla.signif <- lohhla.signif[lohhla.signif$PVal_unique<0.01,]

lohhla.patients <- analyseLohhla(lohhla.master)

# If mutations were also detected, compare LOH and mutation status
lohhla.patients$MUT <- sapply(lohhla.patients$ID, function(x) x %in% muthla.signif$individual)

lohhla.patients.sub <- na.omit(subset(lohhla.patients, CN==T))

# Measure and plot the association between subtypes and LOH/mutation
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


t3 <- fisher.test(table(lohhla.patients.sub$LOH, lohhla.patients.sub$MUT))
p3 <- ggplot(na.omit(lohhla.patients), aes(x=LOH, fill=MUT)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t3$p.value,5)), x='LOH', y='')


t5 <- fisher.test(table(lohhla.patients$MSI, lohhla.patients$MUT))
p5 <- ggplot(na.omit(lohhla.patients), aes(x=MSI, fill=MUT)) +
  geom_bar(stat='count', position='fill') + scale_fill_manual(values=c('skyblue4', 'firebrick3')) +
  scale_y_continuous(labels=percent_format()) +
  labs(title=paste0('Fisher exact test p-value: ',round(t4$p.value,5)), x='MSI high', y='')

# Heterogeneity of loci ---------------------------------------------------
# Investigate how often does LOH occur in 1/2/3 of the possible alleles
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


############################################################################
# HLA mutations -----------------------------------------------------------
############################################################################

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

# Read in all polysolver HLA mutation predictions
noMutPreds <- c()
noMut <- c()
noIndel <- c()

  hlaList <- list.files(paste0('~/CRCdata/HLA_MUT/TCGA_', canc), pattern='*.mutect.unfiltered.annotated')
    if (length(hlaList)==0){noMutPreds = c(noMutPreds, sample)}
    for (hla in hlaList){
      muthla <- tryCatch(
        {read.table(paste0('~/CRCdata/HLA_MUT/TCGA_',canc,'/',hla), stringsAsFactors = F, header=F, skip=1, sep='\t')},
        error=function(e){return(NA)}
      )
      if (!is.na(muthla)){
        names(muthla) <- correctHeader
        muthla.master <- rbind(muthla.master, muthla)
      }
      else{noMut = c(noMut, substr(hla,1,12))}
    }
    hlaList2 <- list.files(paste0('~/CRCdata/HLA_MUT/TCGA_', canc), pattern='*.strelka_indels.unfiltered.annotated')
    for (hla in hlaList2){
    indels.tmp <- read.table(paste0('~/CRCdata/HLA_MUT/TCGA_',canc,'/',hla), stringsAsFactors = F, header=T, sep='\t')
    if(nrow(indels.tmp)>0){print('Woohoo')}
    else{noIndel = c(noIndel, substr(hla,1,12))}
    }

write.table(muthla.master, file=paste0('~/Dropbox/Code/TCGA/',canc,'/',canc,'_MUTHLA_master_file.txt'),
                   sep='\t',row.names=F,quote=F)
# Pull out nonsynonymous, protein changing mutations
muthla.signif <- subset(muthla.master, protein_change!=-1)
muthla.signif$nonsyn <- sapply(muthla.signif$protein_change, function(x) substr(x,3,3) != substr(x, nchar(x), nchar(x)))
muthla.signif <- subset(muthla.signif, nonsyn)


############################################################################
# Gene expression ---------------------------------------------------------
############################################################################

tcga.tpm <- read.table(paste0('TCGA/',canc,'/',canc,'_RNA_expression.tpm'), stringsAsFactors = F, header=T,row.names=1)
tcga.normal.tpm <- read.table(paste0('TCGA/',canc,'/',canc,'_RNA_normal.tpm'), stringsAsFactors = F, header=T, row.names=1)

# Define immune escape genes in table
pdl <- 'ENSG00000120217'
ctla <- 'ENSG00000163599'

# Get PD-L1 expression and define over-expressed
exprpdl <- sapply(tcga.tpm[pdl,], function(x) log10(x+1))
exprpdl.norm <- sapply(tcga.normal.tpm[pdl,], function(x) log10(x+1))
exprpdl.norm <- exprpdl.norm[exprpdl.norm<0.9] #get rid of clear outlier
#plot normal PDL expression
ggplot(data.frame(value=exprpdl), aes(x=value)) + geom_density(colour='firebrick3') +
  geom_density(data=data.frame(value=exprpdl.norm), aes(x=value)) +
  geom_vline(xintercept=(mean(exprpdl.norm)+2*sd(exprpdl.norm)))

pdl.over <- names(exprpdl)[exprpdl > (mean(exprpdl.norm)+2*sd(exprpdl.norm))]

# Define CTLA-4 expression and define over-expressed
exprctla <- sapply(tcga.tpm[ctla,], function(x) log10(x+1))
exprctla.norm <- sapply(tcga.normal.tpm[ctla,], function(x) log10(x+1))
exprctla.norm <- exprctla.norm[exprctla.norm<0.9] #get rid of clear outlier
ggplot(data.frame(value=exprctla), aes(x=value)) + geom_density(colour='firebrick3') +
  geom_density(data=data.frame(value=exprctla.norm), aes(x=value)) +
  geom_vline(xintercept=(mean(exprctla.norm)+2*sd(exprctla.norm)))

ctla.over <- names(exprctla)[exprctla > (mean(exprctla.norm)+2*sd(exprctla.norm))]

#CTLA vs PDL
ggplot(data.frame(c=exprctla,p=exprpdl), aes(x=p, y=c)) + geom_point() +
  geom_hline(yintercept=(mean(exprctla.norm)+2*sd(exprctla.norm))) +
  geom_vline(xintercept=(mean(exprpdl.norm)+2*sd(exprpdl.norm)))


#############################################################################
# Compile escape table and make summary plots -----------------------------
#############################################################################


lohhla.master$region <- gsub('.tumour', '', lohhla.master$region)
lohhla.ai <- subset(lohhla.master, PVal_unique<0.01)
lohhla.loh <- subset(lohhla.ai, Label=='LOH')
muthla.uniq <- muthla.signif[!duplicated(paste0(muthla.signif$individual, muthla.signif$contig)),]
#escape.df <- data.frame(Patient = unique(lohhla.master$region),
#                        HLA_LOH = NA, HLA_MUT=NA, B2M_MUT=0, AI=NA,
#                        PDL1 = NA, CTLA4=NA,MSI=NA)
escape.df <- read.table(paste0('TCGA/',canc,'/',canc,'_escape_master_file.txt'),
                        sep='\t', header=T, stringsAsFactors = F)
clin.df <- read.table(paste0('TCGA/',canc,'/',canc,'_clinical_master_file.txt'),
          sep='\t', header=T, stringsAsFactors = F)

for (pat in unique(lohhla.ai$region)){ escape.df[escape.df$Patient==pat, 'AI'] <- paste0(lohhla.ai[lohhla.ai$region==pat,'LossAllele'], collapse=',') }

for (pat in unique(lohhla.loh$region)){ escape.df[escape.df$Patient==pat, 'HLA_LOH'] <- paste0(lohhla.loh[lohhla.loh$region==pat,'LossAllele'], collapse=',') }

for (pat in unique(muthla.uniq$individual)){ escape.df[escape.df$Patient==pat, 'HLA_MUT'] <- paste0(muthla.uniq[muthla.uniq$individual==pat,'contig'], collapse=',') }

escape.df$B2M_MUT <- clin.df[match(escape.df$Patient, clin.df$Patient), 'B2M']
escape.df$WDFY4_MUT <- clin.df[match(escape.df$Patient, clin.df$Patient), 'WDFY4']
escape.df$PDL1 <- exprpdl[match(gsub('-', '.',escape.df$Patient), names(exprpdl))] > mean(exprpdl.norm)+2*sd(exprpdl.norm)
escape.df$CTLA4 <- exprctla[match(gsub('-', '.',escape.df$Patient), names(exprctla))] > mean(exprctla.norm)+2*sd(exprctla.norm)
escape.df$MSI <- clin.df[match(escape.df$Patient, clin.df$Patient), 'MSI']

write.table(escape.df, paste0('TCGA/',canc,'/',canc,'_escape_master_file.txt'), sep='\t', row.names=F,quote=F)

# HLA ggplot heatmap ---------------------------------------------------------

#Only parts that have an LOH:
lohhla.master.sub <- subset(lohhla.master, startsWith(region, 'Set9') | startsWith(region, 'Set6') | startsWith(region, 'Set8') | startsWith(region, 'Set10') | startsWith(region, 'Set5') | startsWith(region, 'Polyp.9'))
lohhla.master.sub <- subset(lohhla.master.sub, !(region %in% c('Set10_92.mkdub')))

lohhla.df <- data.frame('Region'=character(), 'Allele' = character(), 'CopyNumber'=numeric(), 'p.value'=numeric(),stringsAsFactors = F)

for (i in 1:nrow(lohhla.master.sub)){
  x <- lohhla.master.sub[i,]
  #r <- substr(x$region,1,nchar(x$region)-6)
  r <- x$region
  allele <- paste0('HLA-',toupper(substr(x$HLA_A_type1,5,5)))
  p <- x$PVal_unique
  val <- ifelse(x$PVal_unique<0.05, min(x$HLA_type1copyNum_withBAFBin, x$HLA_type2copyNum_withBAFBin), 1)
  lohhla.df[i,] <- c(r, allele, val, p)
}
lohhla.df$CopyNumber <- as.numeric(lohhla.df$CopyNumber)<0.5
lohhla.df$p.value <- as.numeric(lohhla.df$p.value)
lohhla.df <- lohhla.df[order(lohhla.df$Region),]

pLOH <- ggplot(lohhla.df, aes(y=Allele, x=Region, fill=CopyNumber)) + geom_tile() +
  theme_bw() + scale_fill_manual(values=c('grey80','#d53e4f'), labels=c('No', 'Yes', '')) +
  labs(fill='HLA_LOH') + theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank(), text=element_text(size=16))
scale_fill_gradientn(colours=c('red4', 'red4', 'red4', 'red4', 'brown3', 'mistyrose2', 'bisque', 'cornsilk3'), limits=c(-2.5, 2), na.value='grey70')

