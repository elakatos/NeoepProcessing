
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
source('rfunctions_postprocessing.R')

############################################################################
# Read in patient information ---------------------------------------------
############################################################################

# Process clinical data/ MSI status ---------------------------------------

clin.df <- read.table(paste0('TCGA/',canc,'/',canc,'_clinical_master_file.txt'), sep='\t', header=T, stringsAsFactors = F)

# CRC
# CRC has hypermutation from TRONCO analysis as well as Mantis analysis
clin.df$Hypermut <- clin.df$Hypermut_tronco
clin.df[is.na(clin.df$Hypermut_tronco),'Hypermut'] <- 1* clin.df[is.na(clin.df$Hypermut),'NumMut']>2000

#MSI is taken based on mantis wherever no data is available from tronco and from tronco for all other points
clin.df[is.na(clin.df$MSI_tronco), 'MSI'] <- ifelse((clin.df[is.na(clin.df$MSI_tronco), 'MSI_mantis'] > 0.5), 'MSI-H', 'MSI-L')
clin.df[is.na(clin.df$MSI), 'MSI'] <- clin.df[is.na(clin.df$MSI), 'MSI_tronco']
clin.df[(!is.na(clin.df$MSI)) & (clin.df$MSI=='MSS'),'MSI'] <- 'MSI-L'
#Ones with contradicting information from tronco&mantis are NAd out
clin.df[(clin.df$MSI_tronco %in% c('MSI-L', 'MSS')) & (clin.df$MSI_mantis > 0.5),'MSI'] <- NA

#Strong signature evidence is taken to identify POLE or definite MSI cases
clin.df[clin.df$Sig6>1000, 'MSI'] <- 'MSI-H' #only very strong Sig6 signal is taken into account
clin.df[clin.df$Sig10>1500, 'MSI'] <- 'POLE'

#STAD
clin.df[,'Hypermut'] <- 1* clin.df[,'NumMut']>1200
clin.df[, 'MSI'] <- ifelse((clin.df[, 'MSI_mantis'] > 0.4), 'MSI', 'MSS')
clin.df[clin.df$Sig10>1500, 'MSI'] <- 'POLE'

#UCEC
clin.df[,'Hypermut'] <- 1* clin.df[,'NumMut']>800
clin.df[, 'MSI'] <- ifelse((clin.df[, 'MSI_mantis'] > 0.45), 'MSI', 'MSS')
clin.df[clin.df$Sig10>800, 'MSI'] <- 'POLE'
# ones with contradicting information are NAd out
clin.df[clin.df$Patient=='TCGA-AX-A05Y','MSI'] <- NA

# Process sample quality and immune status --------------------------------

#ploidy.df <- read.delim('TCGA/CRC_ploidy_master.txt',stringsAsFactors = F)
goodSamples <- scan(paste0('TCGA/',canc,'/',canc,'_goodSamples.txt'), what='character()')

tcga.tpm <- read.table(paste0('TCGA/',canc,'/',canc,'_RNA_expression.tpm'), stringsAsFactors = F, header=T)

#T-cell associated genes
tcag <- c('ENSG00000108691', 'ENSG00000277632', 'ENSG00000275302', 'ENSG00000138755', 'ENSG00000169245',
          'ENSG00000153563', 'ENSG00000241106', 'ENSG00000242574', 'ENSG00000204252', 'ENSG00000113088',
          'ENSG00000163600', 'ENSG00000125347')
exprtcag <- log10(tcga.tpm[tcag,]+1)
tcavg <- apply(exprtcag, 2, mean)

clin.df$TCellScore <- tcavg[match(clin.df$Patient, gsub('\\.','-',names(tcavg)))]


############################################################################
# Total burden analysis ---------------------------------------------------
############################################################################

dir = paste0('~/CRCdata/TCGA_',canc)
prefix='Total'

# Read in epitope table
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample','LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epTable <- subset(epTable, Novelty==1)

recoTable <- read.table(paste0(dir,'/Neopred_results/PredictedRecognitionPotentials.txt'),
                        stringsAsFactors = F, header=T)

# Filter according to RecognitionPotential table and sample quality
epTable$AntigenID <- paste0(epTable$Sample,':',epTable$peptide)
recoTable$AntigenID <- paste0(recoTable$Sample,':',recoTable$MutantPeptide)
# plot recognition potential distribution
ggplot(recoTable, aes(x=NeoantigenRecognitionPotential)) + geom_density() +
  theme_mypub() + scale_x_log10(limits=c(1e-6,1e3)) 

recoTable.imm <- subset(recoTable, NeoantigenRecognitionPotential>1e-1)
epTable.imm <- subset(epTable, AntigenID %in% recoTable.imm$AntigenID)

write.table(epTable.imm, file=paste0(dir,'/Neopred_results/',prefix,'.neoantigens.filtered.txt'),
            sep='\t', row.names=F, quote=F)


############################################################################
# generate CCFs of mutations

epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.filtered.txt'), header=T,
                      sep = '\t',stringsAsFactors = F)
cna.df <- subset(meta.df, Patient %in% unique(epTable$Sample))

allTotVAFs <- data.frame(matrix(vector(), ncol=6))

for(sample in unique(epTable$Sample)){
#sampleFileEx <- paste0(dir, '/avannotated/',sample,'.avannotated.exonic_variant_function')
#exonic <- readExonicFile(sampleFileEx)
vcf <- read.table(paste0(dir, '/VCF/',sample,'.vcf'), sep='\t', stringsAsFactors = F)
names(vcf)[c(1,2,11)] <- c('Chr', 'Start', 'Region_0')

#tmpVAF <- data.frame(Sample=sample, Chr=exonic$Chrom, Start=exonic$Start,  LineID=exonic$LineID, VAF=computeVafAD(exonic, 'Region_0'))
tmpVAF <- data.frame(Sample=sample, Chr=vcf$Chr, Start=vcf$Start,  VAF=computeVafAD(vcf, 'Region_0'))


cnafile <- cna.df[cna.df$Patient==sample,'FileName']
if (length(cnafile)==0){next}
cna <- read.table(paste0('~/CRCdata/TCGA_',canc,'/',canc,'_CNA/',cnafile), sep='\t', header=T, stringsAsFactors = F)
cna$Chromosome <- paste0('chr', cna$Chromosome)
cna <- subset(cna,Num_Probes > 30)
cna$Segment_Mean <- (((2^(cna$Segment_Mean))*2 -2 )/pur.df[pur.df$Patient==sample,'Purity'] + 2)
#cna$Segment_Mean <- (2^(cna$Segment_Mean))

tmpVAF$CN <- sapply(1:nrow(tmpVAF), function(i) getCNAofMut(cna,tmpVAF[i,]))
tmpVAF$CCF <- (tmpVAF$VAF*tmpVAF$CN)*(1/pur.df[pur.df$Patient==sample,'Purity'])
allTotVAFs <- rbind(allTotVAFs, tmpVAF)
}

write.table(allTotVAFs, file=paste0('TCGA/',canc,'/',canc,'_allVAF_master_file.txt'),sep='\t',quote=F, row.names=F)

# Also generate an exonic file table
allMutVAFs <- data.frame(matrix(vector(),ncol=6))
for(sample in unique(epTable$Sample)){
sampleFileEx <- paste0(dir, '/Neopred_results/avannotated/',sample,'.avannotated.exonic_variant_function')
exonic <- readExonicFile(sampleFileEx)

tmpVAF <- data.frame(Sample=sample, Chr=exonic$Chrom, Start=exonic$Start,  LineID=exonic$LineID)
tmpVAF[,c('VAF','CN','CCF')] <- allTotVAFs[match(paste0(tmpVAF$Sample,tmpVAF$Chr,tmpVAF$Start),paste0(allTotVAFs$Sample,allTotVAFs$Chr,allTotVAFs$Start)),c('VAF','CN','CCF')]
allMutVAFs <- rbind(allMutVAFs, tmpVAF)
}

allMutVAFs$Gene <- NA
#annotate exonic table with gene
for (samp in unique(allMutVAFs$Sample)){
  exn <- read.table(paste0(dir,'/Neopred_results/avannotated/',samp,'.avannotated.exonic_variant_function'),
                    stringsAsFactors = F, sep='\t')
  allMutVAFs[allMutVAFs$Sample==samp,]$Gene <- sapply(exn$V3, function(x) unlist(strsplit(x, ':'))[1])
}

write.table(allMutVAFs, file=paste0('TCGA/',canc,'/',canc,'_exonicVAF_master_file.txt'),sep='\t',quote=F, row.names=F)


# Add ccf info to epitope table
epTable[,c('VAF','CN','CCF')] <- allMutVAFs[match(paste0(epTable$Sample, epTable$LineID),paste0(allMutVAFs$Sample, allMutVAFs$LineID)),c('VAF','CN','CCF')]

write.table(epTable, file=paste0(dir,'/Neopred_results/',prefix,'.neoantigens.filtered.txt'),sep='\t',quote=F, row.names=F)
