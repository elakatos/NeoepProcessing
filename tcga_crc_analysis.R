
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
source('rfunctions_postprocessing.R')

############################################################################
# Read in patient information ---------------------------------------------
############################################################################

# Process clinical data/ MSI status ---------------------------------------

clin.df <- read.table('~/TCGA/CRC_clinical_master_file.txt', sep='\t', header=T, stringsAsFactors = F)

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

# Plots for sanity check
ggplot(clin.df[!is.na(clin.df$Hypermut_tronco),], aes(x=NumMut, fill=paste(Hypermut_tronco))) + geom_density(alpha=0.3, adjust=0.25) +scale_x_continuous(trans='log10') +
  geom_vline(xintercept=2000)
ggplot(clin.df[!is.na(clin.df$MSI_tronco),], aes(x=MSI_mantis, fill=paste0(MSI_tronco))) + geom_density(alpha=0.5) +
  geom_vline(xintercept=0.5)
ggplot(clin.df[!is.na(clin.df$MSI_tronco),], aes(x=Sig6+1, fill=paste0(MSI_tronco))) + geom_density(alpha=0.5) + scale_x_continuous(trans='log10') +
  geom_vline(xintercept=150)
ggplot(clin.df[!is.na(clin.df$Hypermut_tronco),], aes(x=Sig10+1, fill=paste0(MSI_tronco,Hypermut_tronco))) + geom_density(alpha=0.5) + scale_x_continuous(trans='log10') +
  geom_vline(xintercept=1500)
ggplot(clin.df, aes(x=Sig6, y=MSI_mantis, colour=MSI_tronco)) + geom_point() +
  geom_hline(yintercept=0.4, colour='grey50') + geom_vline(xintercept=1000, colour='grey50') + geom_vline(xintercept=150, colour='grey70')



# Process sample quality and immune status --------------------------------

ploidy.df <- read.delim('~/TCGA/CRC_ploidy_master.txt',stringsAsFactors = F)
goodSample <- scan('~/TCGA/CRC_goodSamples.txt', what='character()')

tcga.tpm <- read.table('~/TCGA/CRC_RNA_expression.tpm', stringsAsFactors = F, header=T)

#T-cell associated genes
tcag <- c('ENSG00000108691', 'ENSG00000277632', 'ENSG00000275302', 'ENSG00000138755', 'ENSG00000169245',
          'ENSG00000153563', 'ENSG00000241106', 'ENSG00000242574', 'ENSG00000204252', 'ENSG00000113088',
          'ENSG00000163600', 'ENSG00000125347')
exprtcag <- log10(tcga.tpm[tcag,]+1)
tcavg <- apply(exprtcag, 2, mean)


############################################################################
# Total burden analysis ---------------------------------------------------
############################################################################

dir = '~/CRCdata/TCGA_CRC'
prefix='Total'

# Read in epitope table
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.txt'), header=F,
                      sep = '\t',stringsAsFactors = F, fill=T)
names(epTable) <- c('Sample','LineID', 'Chrom', 'Start',
                    'RefAll', 'AltAll', 'Gene', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                    'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score','Affinity', 'Rank', 'Cand', 'BindLevel', 'Novelty')
epTable <- subset(epTable, Novelty==1)

epTable$AntigenID <- apply(epTable, 1,function(x) paste0(x['Sample'], ':',x['Identity'], ':',x['peptide'] ) )

# Read in filter and according to RecognitionPotential table and sample quality
recoTable <- read.table(paste0(dir,'/Neopred_results/PredictedRecognitionPotentials.txt'),
                        stringsAsFactors = F, header=T)
recoTable$AntigenID <- apply(recoTable, 1,function(x) paste0(x['Sample'], ':',x['Mutation'], ':',x['MutantPeptide'] ) )
recoTable.imm <- subset(recoTable, NeoantigenRecognitionPotential>1e-1)
epTable.imm <- subset(epTable, AntigenID %in% recoTable.imm$AntigenID)

epTable <- subset(epTable.imm, Sample %in% goodSamples)

# Generate summary table with mutation numbers
summaryTable <- processSummaryOfSampleSet(dir, epTable, prefix)

summaryTable$MutRatio <- summaryTable$Total/summaryTable$Total_MUT
summaryTable$NeoepRatio <- summaryTable$Neoep/summaryTable$Total

summaryTable$MSI <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'MSI']
summaryTable$Hyp <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Hypermut']
summaryTable$Purity <- ploidy.df[match(row.names(summaryTable), row.names(ploidy.df)),'tumorPurity']
summaryTable$B2M <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'B2M']>0
summaryTable$VS <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'VS']
summaryTable$Age <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Age']/365
summaryTable$Race <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Race']
summaryTable$Cancer <- clin.df[match(row.names(summaryTable), clin.df$Patient), 'Cancer']
summaryTable$Stage <- gsub('[abc]', '', clin.df[match(row.names(summaryTable), clin.df$Patient), 'Stage'])
summaryTable[summaryTable=='not reported'] <- NA

immTable <- read.table('~/TCGA/data_CRC_fullimmune.txt',stringsAsFactors = F, header=T, sep='\t')
summaryTable$IP <- immTable[match(row.names(summaryTable), immTable$sampleID), 'Immune.phenotype']

# Compare neoep-mutation ratio in subgroups
summaryTable.sub <- subset(summaryTable, MSI != 'POLE')
summaryTable.sub <- subset(summaryTable.sub, Total_MUT > 30)
mlab <- 'All non-POLE tumours'
ylb <- 'RATIO of neoantigen associated mutations'

#further possible comparisons:
#summaryTable.sub$MutRatio <- summaryTable.sub$NeoepRatio
#ylb <- 'AVG NUMBER of neo-epitopes from one mutation'
#summaryTable.sub$MutRatio <- summaryTable.sub$Total
#ylb <- 'TOTAL neo-epitope associated mutations'
#summaryTable.sub$MutRatio <- summaryTable.sub$Neoep
#ylb <- 'TOTAL neo-epitopes'
#summaryTable.sub$MutRatio <- summaryTable.sub$Total_MUT
#ylb <- 'TOTAL somatic missense mutations'

p0 <- ggplot(summaryTable.sub, aes(x=Total_MUT, y=MutRatio, colour=MSI)) + geom_point(size=2) + scale_x_continuous(trans='log10') +
  scale_color_manual(values=c('black', 'darkorange3','darkseagreen4')) +
  labs(colour='Hypermutated', x='Total somatic missense mutations', y=ylb) + theme_bw() + guides(color=F)
p1 <- ggplot(summaryTable.sub, aes(x=MSI, y=MutRatio, fill=as.factor(MSI))) + geom_violin() +
  stat_compare_means(comparisons = list(c('MSI-H', 'MSI-L')), aes(label = paste0("p = ", ..p.format..))) +
  guides(fill=F) + labs(y=ylb, title=mlab) + theme_bw() + scale_fill_manual(values=c('darkorange3', 'darkseagreen4', 'purple'))
p2 <- ggplot(summaryTable.sub, aes(x=Hyp, y=MutRatio, fill=as.factor(Hyp))) + geom_violin() +
  stat_compare_means(comparisons = list(c(1,0))) +
  guides(fill=F) + labs(y=ylb, x='Hypermutated', title=mlab)
p3 <- ggplot(summaryTable.sub, aes(x=B2M, y=MutRatio, fill=B2M)) + geom_violin() +
  stat_compare_means(comparisons = list(c('TRUE','FALSE'))) +
  guides(fill=F) + labs(y=ylb, x='Mutated in B2M', title=mlab)
p4 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$VS),], aes(x=VS, y=MutRatio, fill=as.factor(VS))) + geom_violin() +
  stat_compare_means(comparisons = list(c('alive','dead'))) +
  guides(fill=F) + labs(y=ylb, x='Vital status', title=mlab)
p5 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Stage),], aes(x=Stage, y=MutRatio, fill=as.factor(Stage))) + geom_violin() +
  guides(fill=F) + labs(y=ylb, title=mlab)
p6 <- ggplot(summaryTable.sub, aes(x=Age, y=MutRatio, colour=MSI)) + geom_point(size=2) + stat_cor(mapping=aes(x=Age,y=MutRatio,colour=Age),size=5) +
  guides(colour=F) + labs(y=ylb,title=mlab) + theme_bw() + scale_color_manual(values=c('darkseagreen4', 'darkorange3'))
p7 <- ggplot(summaryTable.sub[!is.na(summaryTable.sub$Race),], aes(x=Race, y=MutRatio, fill=Race)) + geom_violin() +
  stat_compare_means(comparisons = list(c('black or african american','white'), c('black or african american','asian'), c('asian','white'))) +
  guides(fill=F) + labs(y=ylb,title=mlab)

# Figures for publication : 
pcrcrat <- p1 + geom_boxplot(width=0.05, fill='grey80') + labs(y='Relative neoantigen burden', title='',x='') +
  theme(text=element_text(size=16))

pcrcage1 <- ggplot(summaryTable.sub, aes(x=Age, y=MutRatio, colour=MSI)) + geom_point(size=2) + stat_cor(mapping=aes(x=Age,y=MutRatio,colour=Age),size=5) +
  guides(colour=F) + labs(y='Relative neoantigen burden') + theme_bw() + scale_color_manual(values=c('darkorange3','darkseagreen4')) +
  scale_y_continuous(limits=c(0.0,0.25)) + scale_x_continuous(breaks=c(30,60,90)) +
  theme(text=element_text(size=16))


############################################################################
# Immune escape analysis --------------------------------------------------
############################################################################

# Read in compiled escape table
escape.df <- read.table('TCGA/CRC_escape_master_file.txt',stringsAsFactors = F, sep='\t', header=T)
escape.df <- escape.df[!is.na(escape.df$PDL1),]
escape.df <- subset(escape.df, FULL_INFO)
escape.df <- subset(escape.df, Patient %in% goodSamples)

# Label different escape types based on information in compiled table
escape.df$Escape <- NA
escape.df[ (is.na(escape.df$HLA_LOH) & !is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT==0 & !escape.df$PDL1 & !escape.df$CTLA4), 'Escape' ] <- 'HLA_MUT'
escape.df[ (is.na(escape.df$HLA_LOH) & is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT==0 & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT'
escape.df[ (!is.na(escape.df$HLA_LOH) & is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT==0 & !(escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'HLA_LOH'
escape.df[ (is.na(escape.df$HLA_LOH) & is.na(escape.df$HLA_MUT) & escape.df$B2M_MUT>0 & !escape.df$PDL1 & !escape.df$CTLA4), 'Escape' ] <- 'B2M_MUT'
escape.df[ ((is.na(escape.df$HLA_LOH)) & is.na(escape.df$HLA_MUT) & (escape.df$B2M_MUT>0) & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT_&_B2M'
escape.df[ (!(is.na(escape.df$HLA_LOH)) & is.na(escape.df$HLA_MUT) & (escape.df$B2M_MUT==0) & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT_&_LOH'
escape.df[ ((is.na(escape.df$HLA_LOH)) & !is.na(escape.df$HLA_MUT) & (escape.df$B2M_MUT==0) & (escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'CHECKPOINT_&_HLA'
escape.df[ (!(is.na(escape.df$HLA_LOH)) & (!is.na(escape.df$HLA_MUT) | (escape.df$B2M_MUT>0)) & !(escape.df$PDL1 | escape.df$CTLA4)), 'Escape' ] <- 'LOH_&_MUT'
escape.df[ ((!is.na(escape.df$HLA_LOH)) + (!is.na(escape.df$HLA_MUT)) + (escape.df$B2M_MUT>0) + (escape.df$PDL1 | escape.df$CTLA4))>2, 'Escape' ] <- 'COMBINATION'
escape.df[is.na(escape.df$Escape) & !(is.na(escape.df$AI)),'Escape'] <- 'AI'
escape.df[is.na(escape.df$Escape),'Escape'] <- 'NONE'

# Generate barplot of immune escapes
esc.count <- as.data.frame(t(table(escape.df$Escape, escape.df$MSI)))
esc.count$Var1 <- factor(esc.count$Var1, levels=c('MSI-L', 'MSI-H', 'POLE'))
levels(esc.count$Var1) <- c('MSS (n=280)', 'MSI (n=49)', 'POLE (n=10)')
esc.count <- esc.count[order(esc.count$Var1),]
esc.count$Var2 <- factor(esc.count$Var2, levels=c('B2M_MUT', 'HLA_MUT', 'CHECKPOINT','CHECKPOINT_&_B2M',
                                                  'CHECKPOINT_&_HLA','CHECKPOINT_&_LOH','HLA_LOH','LOH_&_MUT','COMBINATION','AI','NONE'))
esc.count$Freq[esc.count$Var1=='MSS (n=280)'] <- esc.count$Freq[esc.count$Var1=='MSS (n=280)']/280
esc.count$Freq[esc.count$Var1=='MSI (n=49)'] <- esc.count$Freq[esc.count$Var1=='MSI (n=49)']/49
esc.count$Freq[esc.count$Var1=='POLE (n=10)'] <- esc.count$Freq[esc.count$Var1=='POLE (n=10)']/10
esc.count$Var3 <- !(esc.count$Var2 %in% c('NONE', 'AI'))

#barchart with escape yes/no for each subtype
cc <- chisq.test(escape.df$Escape %in% c('NONE', 'AI'), escape.df$MSI)
pescbyn <- ggplot(esc.count, aes(x=Var1, y=Freq ,fill=Var3)) +
  geom_bar(stat='identity', position='fill', width=0.75) +
  labs(fill='Immune\nescape', x='', y='Proportion of samples') + 
  scale_fill_manual(values=c('grey70', '#c35071')) +
  theme_bw() + theme(text=element_text(size=16), axis.title.x=element_blank()) + scale_y_continuous(labels=percent_format(), limits=c(0,1.08)) +
  scale_x_discrete(labels=c('MSS (n=280)', 'MSI (n=49)', 'POLE (n=10)')) +
  geom_signif(y_position=c(1.03), xmin=c(1), xmax=c(3),
              annotation=cc$p.value, tip_length=0, size=0.5)

#generate significance annotation for barplot
signif <- sapply(unique(escape.df$Escape), function(x) chisq.test(escape.df$MSI, escape.df$Escape==x)$p.value)

esc.annot <- esc.count; esc.annot$Label <- NA

esc.annot$Label[esc.annot$Var2 %in% names(signif)[signif<0.05]] <- '*'
esc.annot$Label[esc.annot$Var2 %in% names(signif)[signif<0.01]] <- '**'
esc.annot$Label[esc.annot$Var2 %in% names(signif)[signif<0.001]] <- '***'


pescb2 <- ggplot(esc.annot, aes(x='', y=Freq ,fill=Var2)) + facet_wrap(Var1~.,nrow=3) +
  geom_bar(stat='identity', position='dodge', width=1) +
  labs(fill='Escape mechanism', x='', y='Proportion of samples') + 
  scale_fill_manual(values=c('#7577c9','#3288bd',
                             '#e2b01d', '#7abf9f', '#85bb59','#e9813d',
                             '#d53e4f', '#aa4c9a', '#a77955',
                             '#e1b0b0','grey80'),
                    labels=c('B2M mut','HLA mut',
                             'Checkpoint', 'Checkpoint & B2M', 'Checkpoint & HLA','Checkpoint & LOH',
                             'LOH','LOH & HLA/B2M', 'Combination', 'Allelic imbalance', 'NONE')) +
  theme_bw() + theme(text=element_text(size=16), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(labels=percent_format())

# check value to annotate Chi-squared test for checkpoint vs not
print(chisq.test(escape.df$MSI, grepl('CHECKPOINT',escape.df$Escape))$p.value)
print(chisq.test(escape.df$MSI, escape.df$Escape))

pescb_signif <- pescb2 + facet_wrap(.~Var1,nrow=3) + geom_text(aes(x=0.5+(2*as.numeric(Var2)-1)/22,y=Freq+0.005,label=Label)) +
  geom_signif(y_position=c(0.34), xmin=c(0.5+2/11), xmax=c(0.5+6/11),
              annotation=c("***"), tip_length=0, size=0.5) + theme(text=element_text(size=14), legend.text=element_text(size=11)) +
  scale_y_continuous(breaks=c(0, 0.15, 0.3), labels=percent_format(), limits=c(0, 0.38))


############################################################################
# Burden and infiltrate analysis ------------------------------------------
############################################################################
# Read in VAF/clonality of neoantigens ------------------------------------


allTotVAFs <- read.table('~/TCGA/CRC_allVAF_master_file.txt',
                         sep='\t', stringsAsFactors = F, header=T)
allMutVAFs <- read.table('~/TCGA/CRC_exonicVAF_master_file.txt',
                         sep='\t', stringsAsFactors = F, header=T)
allMutVAFs$mutID <- apply(allMutVAFs, 1,function(x) paste0(x['Sample'], ':',x['LineID'] ) )

# Assign CCF to recotable mutations
recoTable.imm$mutID <- epTable.imm[match(recoTable.imm$AntigenID,epTable.imm$AntigenID),'mutID']
recoTable.imm$CCF <- allMutVAFs[match(recoTable.imm$mutID, allMutVAFs$mutID),'CCF']
recoTable.imm$CCF[recoTable.imm$CCF > 1] <- 1
recoTable.imm <- subset(recoTable.imm,Sample %in% goodSamples)


# Subclonal antigen burden vs T-cell score --------------------------------

recoTable.sclonal <- subset(recoTable.imm, CCF < 0.6 & CCF > 0.3)

numClonal <- data.frame(Sample = intersect(goodSamples,recoTable$Sample))
numClonal$Count <- sapply(numClonal$Sample, function(s) sum(recoTable.sclonal$Sample==s))
numClonal$tc <- tcavg[match(numClonal$Sample, gsub('\\.', '-',names(tcavg)))]

#take only samples for whom escape information is available
numClonal2 <- numClonal[numClonal$Sample %in% escape.df$Patient,]
numClonal2$MSI <- escape.df[match(numClonal2$Sample,escape.df$Patient), 'MSI']
numClonal2$escape <- escape.df[match(numClonal2$Sample,escape.df$Patient), 'Escape']
numClonal2$escape <- numClonal2$escape %in% c('AI', 'NONE')
numClonal2$tc3 <- cut(numClonal2$tc, seq(min(numClonal2$tc),max(numClonal2$tc), length.out=4))

#plot violin of subclonal burden
ggplot(na.omit(numClonal2), aes(x=tc3, y=Count, fill=tc3)) + geom_violin() +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  theme_mypub() +
  labs(x='T-cell score', y='# of subclonal antigenic mutations') +
  scale_x_discrete(labels=c('low', 'medium', 'high')) + guides(fill=F) +
  stat_compare_means(comparisons = list(c('(0.295,0.85]','(1.4,1.96]'), c('(0.85,1.4]', '(1.4,1.96]'))) +
  scale_y_continuous(trans='log1p', breaks=c(0, 20,500), limits=c(0, 1500))

#for immune escaped and not mutations
ggplot(na.omit(numClonal2), aes(x=!escape, y=Count, fill=tc3)) + geom_violin() +
  scale_fill_manual(values=c('#c35071', 'grey70')) +
  theme_mypub() +
  labs(x='Immune escape', y='# of subclonal antigenic mutations') +
  scale_x_discrete(labels=c('No','Yes')) + guides(fill=F) +
  stat_compare_means(comparisons = list(c(FALSE,TRUE))) +
  scale_y_continuous(trans='log1p', breaks=c(0, 20,500), limits=c(0, 1500))



# Is there a difference in large subclone mutational burden for MSI -------

escape.df$Total_neoep <- NA; escape.df$Subclonal_neoep <- NA

for (samp in goodSamples){
  vaftmp <- subset(allMutVAFs, Sample == samp)
  escape.df[escape.df$Patient==samp,'Total_neoep'] <- sum((vaftmp$mutID %in% epTable.imm$mutID), na.rm=T)
  escape.df[escape.df$Patient==samp,'Subclonal_neoep'] <- sum((vaftmp$mutID %in% epTable.imm$mutID) & (vaftmp$CCF < 0.6) & (vaftmp$CCF > 0.3), na.rm=T)
  escape.df[escape.df$Patient==samp,'Subclonal_all'] <- sum((vaftmp$CCF < 0.6) & (vaftmp$CCF > 0.3), na.rm=T)
}

# Focus on MSI samples
x <- subset(escape.df, MSI=='MSI-H')

#plot subclonal burden for MSI samples
psc <- ggplot(x, aes(x=(Escape %in% c('NONE','AI')), y=Subclonal_neoep, fill=(Escape %in% c('NONE','AI')))) +
  geom_boxplot() +
  stat_compare_means(comparisons=list(c('FALSE', 'TRUE')),label.x.npc = 'centre', method.args= list(alternative='greater')) + guides(fill=F) +
  scale_fill_manual(values=c('#c35071', 'grey70')) + scale_x_discrete(labels=c('Yes (n=42)', 'No (n=7)')) +
  labs(y='Neoantigens in large subclone', x='Immune escape') +
  theme_bw() + theme(text=element_text(size=12)) + scale_y_continuous(limits=c(0, 75))



# Distribution of antigenicity values -------------------------------------

recoTable$tc3 <- numClonal2[match(recoTable$Sample, numClonal2$Sample),'tc3']
recoTable.imm$tc3 <- numClonal2[match(recoTable.imm$Sample, numClonal2$Sample),'tc3']
epTable.imm$tc3 <- numClonal2[match(epTable.imm$Sample, numClonal2$Sample),'tc3']

#plot all recognitionPotential values with cut-off used in article
ggplot(recoTable, aes(x=NeoantigenRecognitionPotential)) + geom_density(fill='#9b8049', alpha=0.7) +
  scale_x_log10(limits=c(1e-6,1e3)) +
  theme_mypub() +
  labs(x='Recognition Potential', y='Density') + geom_vline(xintercept = 1e-1, colour='firebrick', size=1.2, linetype='dashed')


# Get distribution from individual density fits
temp <- data.frame(matrix(vector(),nrow=length(Amids)))
h_breaks <- seq(-0.31, 2.5, by=0.1)
for (samp in unique(epTable.imm$Sample)){
  if (nrow(epTable.imm[epTable.imm$Sample==samp,])>3){
    d <- density(epTable.imm[epTable.imm$Sample==samp,]$invRank, from=-0.305, to=2.45, n=301, adjust=0.75)
    temp <- cbind(temp, d$y); names(temp)[ncol(temp)] <- samp
  }}
temp$invRank <- Amids
allInvRank <- melt(temp,id='invRank')
allInvRank$tc3 <- numClonal2[match(allInvRank$variable, numClonal2$Sample),'tc3']

#aggregate all fits to obtain mean and SD values to plot with shaded regions
invRankStats <- aggregate(allInvRank[,c('value')], list(allInvRank$tc3, allInvRank$invRank), mean)
names(invRankStats) <- c('tc3', 'invRank', 'avg')
invRankStats$SD <- aggregate(allInvRank[,c('value')], list(allInvRank$tc3, allInvRank$invRank), sd)$x

#plot mean antigenicity curve and SD around it
ggplot(invRankStats, aes(x=invRank, y=avg, group=tc3, colour=tc3)) +
  theme_mypub() +
  geom_ribbon(aes(ymin=avg-SD, ymax=avg+SD, fill=tc3), alpha=0.15, colour=NA) + geom_line(size=1.5) +
  scale_colour_manual(values=c('#7b69af','#67a9cf','#027757')) +
  scale_fill_manual(values=c('#7b69af','#67a9cf','#027757')) +
  labs(x='Normalised binding strength', y='Density') + guides(colour=F, fill=F) +
  scale_x_continuous(breaks=c(0, 1, 2), expand=c(0.02,0.02))

#plot on total population's violins too
ggplot(subset(epTable.imm, !is.na(tc3)), aes(x=tc3, y=invRank, fill=tc3)) + geom_violin() +
  theme_mypub() +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  guides(fill=F) + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(),
                         text=element_text(size=10)) +
  scale_y_continuous(breaks=c(0,1,2), limits=c(-0.3, 3.3)) + labs(y='Norm binding strength') +
  stat_compare_means(comparisons=list(c('(0.295,0.85]','(1.4,1.96]'),c('(0.85,1.4]' ,'(1.4,1.96]'),
                                      c('(0.295,0.85]','(0.85,1.4]')),
                     method.args=list(alternative='less'), size=2.5, label='p.signif') #+
geom_boxplot(width=0.05, fill='grey80')



# Pull together samples into VAF curve ------------------------------------

allMutVAFs$Gene <- NA
#annotate mutation table with gene
for (samp in unique(allMutVAFs$Sample)){
  exn <- read.table(paste0('~/CRCdata/TCGA_CRC/avannotated/',samp,'.avannotated.exonic_variant_function'),
                    stringsAsFactors = F, sep='\t')
  allMutVAFs[allMutVAFs$Sample==samp,]$Gene <- sapply(exn$V3, function(x) unlist(strsplit(x, ':'))[1])
}

#read in cell essential genes
cellEss <- scan('~/Documents/cell_essential_genes.txt', what='character()')

fmax=0.65; fmin=0.25
steps <- seq(fmax,fmin,by=(-1e-2))

# Get high quality samples by subsetting the ones that have best purity and ploidy
veryGoodSamples <- row.names(subset(ploidy.df, tumorPurity > 0.75 & tumorPloidy < 3.5))

checkSamples <- subset(escape.df, Escape %in% c('NONE') & MSI=='MSI-L')$Patient
checkSamples <- checkSamples[checkSamples %in% veryGoodSamples]

#cell essential neoantigens, total (including non-exonic), exonic and essential gene
vafEp <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples) & (mutID %in% epTable$mutID) & (Gene %in% cellEss) )$CCF)
vafTot <- na.omit(subset(allTotVAFs, (Sample %in% checkSamples))$CCF)
vafMut <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples))$CCF)
vafEss <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples) & (Gene %in% cellEss))$CCF)

cumvaf <- data.frame(invf = (1/steps),
                     cumvafEp=sapply(steps, function(x) sum(vafEp>=x)),
                     cumvafTot=sapply(steps, function(x) sum(vafTot>=x)),
                     cumvafEx=sapply(steps, function(x) sum(vafMut>=x)),
                     cumvafEss=sapply(steps, function(x) sum(vafEss>=x)))
cumvaf$cumvafEp <- cumvaf$cumvafEp - min(cumvaf$cumvafEp);cumvaf$cumvafEp <- cumvaf$cumvafEp/max(cumvaf$cumvafEp)
cumvaf$cumvafTot <- cumvaf$cumvafTot - min(cumvaf$cumvafTot);cumvaf$cumvafTot <- cumvaf$cumvafTot/max(cumvaf$cumvafTot)
cumvaf$cumvafEx <- cumvaf$cumvafEx - min(cumvaf$cumvafEx);cumvaf$cumvafEx <- cumvaf$cumvafEx/max(cumvaf$cumvafEx)
cumvaf$cumvafEss <- cumvaf$cumvafEss - min(cumvaf$cumvafEss);cumvaf$cumvafEss <- cumvaf$cumvafEss/max(cumvaf$cumvafEss)
cumvaf <- cumvaf[,c('invf','cumvafTot', 'cumvafEx', 'cumvafEss', 'cumvafEp')]

# Plot 1/f VAF curve
ggplot(melt(cumvaf, id='invf'), aes(x=invf*2, y=value, color=variable, shape=variable, size=variable)) +
  geom_line() +
  scale_color_manual(values=c('grey45','#3d9cca',"#8150a3", "#d63027"), labels=c('All', 'Exonic', 'Essential\ngene', 'Neoantigen\nin essential\ngene')) +
  labs(color='Mutations', x='Inverse allelic frequency 1/f', y='Cumulative frequency distribution') +
  theme_bw() + theme(text=element_text(size=14)) + guides(shape=F, size=F) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  scale_x_continuous(breaks=c(4,5,6),labels=c('1/0.25','1/0.2', '1/0.1667')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text=element_text(size=11)) +
  scale_size_manual(values=c(2,1,1,2)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

#check statistical difference between distributions
ks.test(vafTot[vafTot>fmin & vafTot<fmax], vafMut[vafMut>fmin & vafMut<fmax])
ks.test(vafTot[vafTot>fmin & vafTot<fmax], vafEp[vafEp>fmin & vafEp<fmax])
ks.test(vafMut[vafMut>fmin & vafMut<fmax], vafEp[vafEp>fmin & vafEp<fmax])
ks.test(vafMut[vafMut>fmin & vafMut<fmax], vafEss[vafEss>fmin & vafEss<fmax])
ks.test(vafEp[vafEp>fmin & vafEp<fmax], vafEss[vafEss>fmin & vafEss<fmax])


