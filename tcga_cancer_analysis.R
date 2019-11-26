
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
source('rfunctions_postprocessing.R')

############################################################################
# Read in patient information ---------------------------------------------
############################################################################

clin.df <- read.table(paste0('TCGA/',canc,'/',canc,'_clinical_master_file.txt'), sep='\t', header=T, stringsAsFactors = F)
goodSamples <- scan(paste0('TCGA/',canc,'/',canc,'_goodSamples.txt'), what='character()')
pur.df <- read.table('TCGA/ascat_acf_ploidy.tsv', header=T, sep='\t')
pur.df <- subset(pur.df, Cancer_Type_Code==canc)
#pur.df <- subset(pur.df, Cancer_Type_Code %in% c('COAD','READ')) #CRC
pur.df$Patient <- sapply(pur.df$Sample, function(x) gsub('\\.','-',substr(x,1,12)))
names(pur.df)[3] <- 'Purity'

clin.df <- subset(clin.df, Patient %in% goodSamples)

############################################################################
# Immune escape analysis --------------------------------------------------
############################################################################

# Read in compiled escape table
escape.df <- read.table(paste0('TCGA/',canc,'/',canc,'_escape_master_file.txt'),stringsAsFactors = F, sep='\t', header=T)
escape.df <- escape.df[!is.na(escape.df$PDL1),]
escape.df <- subset(escape.df, FULL_INFO) #CRC
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

# STAD has only one POLE sample, should be excluded
#escape.df <- subset(escape.df, MSI!='POLE')

# Generate barplot of immune escapes
esc.count <- as.data.frame(t(table(escape.df$Escape, escape.df$MSI)))
type.tab <- table(escape.df$MSI)
esc.count$Var1 <- factor(esc.count$Var1, levels=c('MSS', 'MSI', 'POLE'))
levels(esc.count$Var1) <- c(paste0('MSS (n=',type.tab['MSS'],')'), paste0('MSI (n=',type.tab['MSI'],')'),paste0('POLE (n=',type.tab['POLE'],')'))
esc.count <- esc.count[order(esc.count$Var1),]
esc.count$Var2 <- factor(esc.count$Var2, levels=c('B2M_MUT', 'HLA_MUT', 'CHECKPOINT','CHECKPOINT_&_B2M',
                                                  'CHECKPOINT_&_HLA','CHECKPOINT_&_LOH','HLA_LOH','LOH_&_MUT','COMBINATION','AI','NONE'))
esc.count$Freq[esc.count$Var1==paste0('MSS (n=',type.tab['MSS'],')')] <- esc.count$Freq[esc.count$Var1==paste0('MSS (n=',type.tab['MSS'],')')]/type.tab['MSS']
esc.count$Freq[esc.count$Var1==paste0('MSI (n=',type.tab['MSI'],')')] <- esc.count$Freq[esc.count$Var1==paste0('MSI (n=',type.tab['MSI'],')')]/type.tab['MSI']
esc.count$Freq[esc.count$Var1==paste0('POLE (n=',type.tab['POLE'],')')] <- esc.count$Freq[esc.count$Var1==paste0('POLE (n=',type.tab['POLE'],')')]/type.tab['POLE']
esc.count$Var3 <- !(esc.count$Var2 %in% c('NONE', 'AI'))

#barchart with escape yes/no for each subtype
cc <- chisq.test(escape.df$Escape %in% c('NONE', 'AI'), escape.df$MSI)
pescbyn <- ggplot(esc.count, aes(x=Var1, y=Freq ,fill=Var3)) +
  geom_bar(stat='identity', position='fill', width=0.75) +
  labs(fill='Immune\nescape', x='', y='Proportion of samples') + 
  scale_fill_manual(values=c('grey70', '#c35071')) +
  theme_mypub()+ scale_y_continuous(labels=percent_format(), limits=c(0,1.08)) +
  scale_x_discrete(labels=c(paste0('MSS\n(n=',type.tab['MSS'],')'), paste0('MSI\n(n=',type.tab['MSI'],')'),paste0('POLE\n(n=',type.tab['POLE'],')'))) +
  geom_signif(y_position=c(1.03), xmin=c(1), xmax=c(3),
              annotation=scientific(cc$p.value,digits=3), tip_length=0, size=0.5)
  #geom_signif(y_position=c(1.03), xmin=c(1), xmax=c(2), #STAD
  #            annotation=scientific(cc$p.value,digits=3), tip_length=0, size=0.5) 

#generate significance annotation for barplot
signif <- sapply(unique(escape.df$Escape), function(x) chisq.test(escape.df$MSI, escape.df$Escape==x)$p.value)

esc.annot <- esc.count; esc.annot$Label <- NA

esc.annot$Label[esc.annot$Var2 %in% names(signif)[signif<0.05]] <- '*'
esc.annot$Label[esc.annot$Var2 %in% names(signif)[signif<0.01]] <- '**'
esc.annot$Label[esc.annot$Var2 %in% names(signif)[signif<0.001]] <- '***'


pescb2 <- ggplot(esc.annot, aes(x='', y=Freq ,fill=Var2)) + facet_wrap(Var1~.,nrow=3) +
  geom_bar(stat='identity', position='dodge', width=1) +
  labs(fill='Escape mechanism', x='', y='Proportion of samples') + 
  theme_bw() + theme(text=element_text(size=16), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(labels=percent_format()) +
  scale_fill_manual(values=c('#7577c9','#3288bd',#CRC
                             '#e2b01d', '#7abf9f', '#85bb59','#e9813d',
                             '#d53e4f', '#aa4c9a', '#a77955',
                             '#e1b0b0','grey80'),
                    labels=c('B2M mut','HLA mut',
                             'Checkpoint', 'Checkpoint & B2M', 'Checkpoint & HLA','Checkpoint & LOH',
                             'LOH','LOH & HLA/B2M', 'Combination', 'Allelic imbalance', 'NONE'))
  #scale_fill_manual(values=c('#7577c9','#3288bd', #STAD
  #                           '#e2b01d', '#7abf9f', '#85bb59','#e9813d',
  #                           '#d53e4f', '#aa4c9a',
  #                           '#e1b0b0','grey80'),
  #                  labels=c('B2M mut','HLA mut',
  #                           'Checkpoint', 'Checkpoint & B2M', 'Checkpoint & HLA','Checkpoint & LOH',
  #                           'LOH','LOH & HLA/B2M', 'Allelic imbalance', 'NONE'))
#scale_fill_manual(values=c('#7577c9','#3288bd', #UCEC
#                           '#e2b01d', '#7abf9f', '#85bb59','#e9813d',
#                           '#d53e4f', '#a77955',
#                           '#e1b0b0','grey80'),
#                  labels=c('B2M mut','HLA mut',
#                           'Checkpoint', 'Checkpoint & B2M', 'Checkpoint & HLA','Checkpoint & LOH',
#                           'LOH','Combination', 'Allelic imbalance', 'NONE')) 


# check value to annotate Chi-squared test for checkpoint vs not
print(chisq.test(escape.df$MSI, grepl('CHECKPOINT',escape.df$Escape))$p.value)
print(chisq.test(escape.df$MSI, escape.df$Escape))

# CRC
pescb_signif <- pescb2 + facet_wrap(.~Var1,nrow=3) +
  geom_text(aes(x=0.5+(2*as.numeric(Var2)-1)/22,y=Freq+0.005,label=Label)) +
  geom_signif(y_position=c(0.34), xmin=c(0.5+2/11), xmax=c(0.5+6/11),
              annotation=c("***"), tip_length=0, size=0.5) +
  theme_mypub() +
  scale_y_continuous(breaks=c(0, 0.15, 0.3), labels=percent_format(), limits=c(0, 0.38))

# UCEC
pescb_signif <- pescb2 + facet_wrap(.~Var1,nrow=3) +
  geom_text(aes(x=0.5+(2*as.numeric(Var2)-1)/20,y=Freq+0.005,label=Label)) +
  geom_signif(y_position=c(0.48), xmin=c(0.5+2/10), xmax=c(0.5+6/10),
              annotation=c("***"), tip_length=0, size=0.5) +
  theme_mypub()+
  scale_y_continuous(breaks=c(0, 0.20, 0.4), labels=percent_format(), limits=c(0, 0.55))

# STAD
pescb_signif <- pescb2 + facet_wrap(.~Var1,nrow=3) +
  geom_text(aes(x=0.5+(2*as.numeric(Var2)-1)/20,y=Freq+0.005,label=Label))


# Immune infiltration compared to escape
escape.df$TC <- clin.df[match(escape.df$Patient, clin.df$Patient),'TCellScore']
ggplot(escape.df, aes(x=!(Escape %in% c('AI','NONE')), fill=!(Escape %in% c('AI','NONE')), y=TC)) +
  geom_violin() + labs(y='T-cell score') +
  theme_mypub() + theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels=c('No escape', 'Escape')) +
  scale_fill_manual(values=c('grey70','#c35071')) + guides(fill=F) +
  stat_compare_means()

ggplot(escape.df[escape.df$MSI=='MSS',], aes(x=!(Escape %in% c('AI','NONE')), fill=!(Escape %in% c('AI','NONE')), y=TC)) +
  geom_violin() + labs(y='T-cell score') +
  theme_mypub() + theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels=c('No escape', 'Escape')) +
  scale_fill_manual(values=c('grey70','#c35071')) + guides(fill=F) +
  stat_compare_means()

ggplot(escape.df, aes(x=MSI, y=TC, fill=MSI)) +
  geom_violin() + theme_mypub() +
  stat_compare_means()

sumTable <- read.table(paste0('TCGA/',canc,'/',canc,'_summary_table.tsv'), sep='\t',header=T, stringsAsFactors = F)
escape.df$TCS3 <- sumTable[match(escape.df$Patient,sumTable$Patient),'TCS3_gen']
escape.df$TCS3 <- factor(escape.df$TCS3, levels=c('low','medium','high'))

cc <- chisq.test(escape.df$Escape %in% c('NONE', 'AI'), escape.df$TCS3)
ggplot(escape.df, aes(x = TCS3, fill=!(Escape %in% c('AI','NONE')))) +
  geom_bar(pos='fill', width=0.75) + labs(x='T-cell score', y='Proportion of samples') +
  theme_mypub() +
  scale_y_continuous(labels=percent_format(), limits=c(0, 1.06)) +
  scale_fill_manual(values=c('grey70','#c35071')) + guides(fill=F) +
  geom_signif(aes(y=CYT),y_position=c(1.03), xmin=c(1), xmax=c(3),
              annotation=scientific(cc$p.value,digits=3), tip_length=0, size=0.5)

############################################################################
# Total burden analysis ---------------------------------------------------
############################################################################

dir = paste0('~/CRCdata/TCGA_',canc)
prefix='Total'

# Read in epitope table
epTable <- read.table(paste0(dir, '/Neopred_results/',prefix,'.neoantigens.filtered.txt'), header=T,
                      sep = '\t',stringsAsFactors = F)
epTable <- subset(epTable, Sample %in% goodSamples)
epTable$mutID <- paste0(epTable$Sample, ':', epTable$LineID)
epTable.dedup <- epTable[!duplicated(paste0(epTable$Sample, epTable$LineID)),]

# plot CCF distribution
ggplot(epTable, aes(x=CCF)) + geom_histogram(bins=100) + theme_mypub() +
  geom_vline(xintercept = c(0.6,1)) + scale_x_continuous(limits=c(0,3))
epTable.clonal <- subset(epTable, CCF>0.65) #STAD
epTable.clonal <- subset(epTable, CCF>0.6) #CRC
epTable.clonal <- subset(epTable, CCF>0.55) #UCEC

epTable$Clonal <- ifelse(epTable$CCF>=0.65, 'clonal', 'subclonal')


# Generate summary table
sumTable <- clin.df[,c('Patient','MSI','Cancer','TCellScore')]
sumTable$Neoantigen_mutation <- sapply(sumTable$Patient, function(x) sum(epTable.dedup$Sample==x))
fastas <- gsub('\\.tmp.10.fasta','',list.files(paste0(dir,'/Neopred_results/fastaFiles/'), pattern='*.tmp.10.fasta'))
sumTable$Total_missense_mutation <- sapply(sumTable$Patient, function(x) ifelse(x %in% fastas,
                                                                                getTotalMutFromFasta(paste0(dir,'/Neopred_results/'),x), NA ))
sumTable$Clonal_neoantigen <- sapply(sumTable$Patient, function(x) sum(epTable.clonal$Sample==x))

sumTable$TCS3 <- cut(sumTable$TCellScore,breaks=3,labels=c('low', 'medium','high'))
sumTable$Escape <- escape.df[match(sumTable$Patient, escape.df$Patient),'Escape']

sumTable[,c('Ploidy','Purity')] <- pur.df[match(sumTable$Patient, pur.df$Patient),c('Ploidy','Purity')]

write.table(sumTable, file=paste0('TCGA/',canc,'/',canc,'_summary_table.tsv'),
            sep='\t', row.names=F, quote=F)

# Read in if already generated
sumTable <- read.table(paste0('TCGA/',canc,'/',canc,'_summary_table.tsv'), sep='\t',header=T, stringsAsFactors = F)
sumTable$TCS3 <- factor(sumTable$TCS3_gen, levels=c('low','medium','high'))
# Plot relative burden
sumTable.sub <- subset(sumTable,Total_missense_mutation>30 )
ggplot(sumTable.sub[!is.na(sumTable.sub$TCS3),], aes(x=TCS3,y=Neoantigen_mutation/Total_missense_mutation,fill=TCS3)) +
  geom_violin() +
  geom_boxplot(width=0.05,fill='grey80') +
  theme_mypub() +
  labs(y='Relative neoantigen burden', x='T-cell score') +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  stat_compare_means(comparisons=list(c('medium','high'),c('low','medium'))) +
  guides(fill=F) #+geom_jitter(width=0.1, aes(colour=Cancer), size=0.2)

ggplot(sumTable.sub[sumTable.sub$MSI!='POLE',], aes(x=MSI,y=Neoantigen_mutation/Total_missense_mutation,fill=MSI)) +
  geom_violin() +
  geom_boxplot(width=0.05,fill='grey80') +
  theme_mypub() +
  scale_fill_manual(values=c('darkorange3','darkseagreen4','grey50')) +
  labs(y='Relative neoantigen burden', x='') +
  stat_compare_means() +
  guides(fill=F)

ggplot(sumTable.sub[!is.na(sumTable.sub$Escape),], aes(x=!(Escape %in% c('AI','NONE')),y=Neoantigen_mutation/Total_missense_mutation,fill=!(Escape %in% c('AI','NONE')))) +
  geom_violin() +
  geom_boxplot(width=0.05,fill='grey80') +
  theme_mypub() +
  scale_x_discrete(labels=c('No escape','Escape')) +
  scale_fill_manual(values=c('grey70','#c35071')) +
  labs(y='Relative neoantigen burden', x='') +
  stat_compare_means() +
  guides(fill=F)

############################################################################
# Burden and infiltrate analysis ------------------------------------------
############################################################################

# Subclonal antigen burden vs T-cell score --------------------------------

epTable.sclonal <- subset(epTable, CCF < 0.6) #CRC
epTable.sclonal <- subset(epTable, CCF < 0.65) #STAD
epTable.sclonal <- subset(epTable, CCF < 0.55) #UCEC

numClonal <- data.frame(Sample = intersect(epTable$Sample, goodSamples))

numClonal$Count <- sapply(numClonal$Sample, function(s) sum(epTable[epTable$Clonal=='subclonal',]$Sample==s, na.rm=T))
numClonal$Count2 <- sapply(numClonal$Sample, function(s) sum(epTable[epTable$Clonal=='clonal',]$Sample==s, na.rm=T))
numClonal$TCS3 <- sumTable[match(numClonal$Sample,sumTable$Patient),'TCS3']
numClonal$MSI <- sumTable[match(numClonal$Sample,sumTable$Patient),'MSI']
numClonal$Escape <- sumTable[match(numClonal$Sample,sumTable$Patient),'Escape']
numClonal$EscapeYN <- ifelse(numClonal$Escape %in% c('AI','NONE'), 'No escape','Escape')
numClonal$TCS3 <- factor(numClonal$TCS3, levels=c('low','medium','high'))

#plot violin of subclonal burden
ggplot(na.omit(numClonal), aes(x=TCS3, y=Count, fill=TCS3)) + geom_violin() +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  theme_mypub() +
  labs(x='T-cell score', y='# of subclonal antigenic mutations') +
  guides(fill=F) +
  stat_compare_means(comparisons=list(c('medium','high'),c('low','high'))) +
  scale_y_continuous(trans='log1p', breaks=c(0, 10,100,1000), limits=c(0, 2000)) +
  geom_dotplot(binaxis='y',stackdir = 'center', dotsize=0.35, binwidth=0.25, fill='black')

ggplot(na.omit(numClonal), aes(x=TCS3, y=log10(Count/Count2+0.01), fill=TCS3)) + geom_violin() +
  scale_fill_manual(values=c('#cac3df','#67a9cf','#02818a')) +
  theme_mypub() +
  labs(x='T-cell score', y='# of subclonal antigenic mutations') +
  guides(fill=F) +
  stat_compare_means(comparisons=list(c('medium','high'),c('low','high'), c('low','medium')))

# for immune escaped and not mutations
tab <- table(na.omit(numClonal)$EscapeYN)
ggplot(na.omit(numClonal), aes(x=EscapeYN, y=Count, fill=EscapeYN)) + geom_violin() +
  scale_fill_manual(values=c('#c35071', 'grey70')) +
  theme_mypub() +
  labs(x='Immune escape', y='# of subclonal antigenic mutations') + guides(fill=F) +
  stat_compare_means(comparisons=list(c('No escape','Escape'))) +
  scale_y_continuous(trans='log1p', breaks=c(0, 10,100,1000), limits=c(0,1200)) +
  scale_x_discrete(labels=c(paste0('Yes (n=',tab[1],')'),
                            paste0('No (n=',tab[2],')'))) +
  geom_dotplot(binaxis='y',stackdir = 'center', dotsize=0.35, binwidth=0.25, fill='black')

# Is there a difference in large subclone mutational burden for MSI

epTable.lgclonal <- subset(epTable, CCF < 0.6 & CCF > 0.3) 
numClonal$LScCount <- sapply(numClonal$Sample, function(s) sum(epTable.lgclonal$Sample==s))
numClonal.sub <- subset(numClonal, MSI=='MSI' & !is.na(Escape))
tab <- table(numClonal.sub$EscapeYN)
ggplot(numClonal.sub, aes(x=EscapeYN, y=LScCount, fill=EscapeYN)) +
  geom_violin() + guides(fill=F) +
  stat_compare_means(comparisons = list(c('No escape', 'Escape')),
                     method.args = list(alternative='less')) +
  scale_fill_manual(values=c('#c35071', 'grey70')) +
  scale_x_discrete(labels=c(paste0('Escape\n(n=',tab[1],')'), paste0('No escape\n(n=',tab[2],')'))) +
  labs(y='Neoantigens in large subclone', x='Immune escape') +
  theme_mypub() + scale_y_continuous(limits=c(0, 180)) +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=0.25, binwidth=10, fill='black')


# Pull together samples into VAF curve ------------------------------------
# Read in mutation data tables and cell essential genes
allTotVAFs <- read.table(paste0('TCGA/',canc,'/',canc,'_allVAF_master_file.txt'),
                         sep='\t', stringsAsFactors = F, header=T)
allMutVAFs <- read.table(paste0('TCGA/',canc,'/',canc,'_exonicVAF_master_file.txt'),
                         sep='\t', stringsAsFactors = F, header=T)
allMutVAFs$mutID <- paste0(allMutVAFs$Sample, ':', allMutVAFs$LineID)
cellEss <- scan('~/Documents/cell_essential_genes.txt', what='character()')

fmax=0.55; fmin=0.3 #STAD & CRC
fmax=0.45; fmin=0.2 #UCEC
steps <- seq(fmax,fmin,by=(-1e-2))

# Get high quality samples by subsetting the ones that have best purity and ploidy
veryGoodSamples <- subset(pur.df, Purity > 0.7 & Ploidy < 3.5)$Patient

checkSamples <- subset(escape.df, (Escape %in% c('NONE')) & MSI=='MSS')$Patient
checkSamples <- checkSamples[checkSamples %in% veryGoodSamples]

#cell essential neoantigens, total (including non-exonic), exonic and essential gene
vafEp <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples) & (mutID %in% epTable$mutID) & (Gene %in% cellEss) )$CCF)
vafEpTot <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples) & (mutID %in% epTable$mutID) )$CCF)
vafTot <- na.omit(subset(allTotVAFs, (Sample %in% checkSamples))$CCF)
vafMut <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples))$CCF)
vafEss <- na.omit(subset(allMutVAFs, (Sample %in% checkSamples) & (Gene %in% cellEss))$CCF)

cumvaf <- data.frame(invf = (1/steps),
                     cumvafEp=sapply(steps, function(x) sum(vafEp>=x)),
                     cumvafEpTot=sapply(steps, function(x) sum(vafEpTot>=x)),
                     cumvafTot=sapply(steps, function(x) sum(vafTot>=x)),
                     cumvafEx=sapply(steps, function(x) sum(vafMut>=x)),
                     cumvafEss=sapply(steps, function(x) sum(vafEss>=x)))
cumvaf$cumvafEp <- cumvaf$cumvafEp - min(cumvaf$cumvafEp);cumvaf$cumvafEp <- cumvaf$cumvafEp/max(cumvaf$cumvafEp)
cumvaf$cumvafEpTot <- cumvaf$cumvafEpTot - min(cumvaf$cumvafEpTot);cumvaf$cumvafEpTot <- cumvaf$cumvafEpTot/max(cumvaf$cumvafEpTot)
cumvaf$cumvafTot <- cumvaf$cumvafTot - min(cumvaf$cumvafTot);cumvaf$cumvafTot <- cumvaf$cumvafTot/max(cumvaf$cumvafTot)
cumvaf$cumvafEx <- cumvaf$cumvafEx - min(cumvaf$cumvafEx);cumvaf$cumvafEx <- cumvaf$cumvafEx/max(cumvaf$cumvafEx)
cumvaf$cumvafEss <- cumvaf$cumvafEss - min(cumvaf$cumvafEss);cumvaf$cumvafEss <- cumvaf$cumvafEss/max(cumvaf$cumvafEss)
cumvaf <- cumvaf[,c('invf','cumvafTot', 'cumvafEx', 'cumvafEss','cumvafEpTot', 'cumvafEp')]

# Plot 1/f VAF curve
ggplot(melt(cumvaf, id='invf'), aes(x=invf*2, y=value, color=variable, shape=variable, size=variable)) +
  geom_line() +
  scale_color_manual(values=c('grey45','#3d9cca',"#8150a3","#b36e8c", "#d63027"), labels=c('All', 'Exonic', 'Essential\ngene', 'Neoantigen','Neoantigen\nin essential\ngene')) +
  labs(color='Mutations', x='Inverse allelic frequency 1/f', y='Cumulative frequency distribution') +
  theme_bw() + theme(text=element_text(size=14)) + guides(shape=F, size=F) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  #scale_x_continuous(breaks=c(4,5,6),labels=c('1/0.25','1/0.2', '1/0.1667')) +
  scale_x_continuous(breaks=c(5,8,10),labels=c('1/0.2', '1/0.125', '1/0.1')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.text=element_text(size=11)) +
  scale_size_manual(values=c(2,1,1,1,2)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

#check statistical difference between distributions
ks.test(vafTot[vafTot>fmin & vafTot<fmax], vafMut[vafMut>fmin & vafMut<fmax])
ks.test(vafTot[vafTot>fmin & vafTot<fmax], vafEp[vafEp>fmin & vafEp<fmax])
ks.test(vafMut[vafMut>fmin & vafMut<fmax], vafEp[vafEp>fmin & vafEp<fmax])
ks.test(vafMut[vafMut>fmin & vafMut<fmax], vafEss[vafEss>fmin & vafEss<fmax])
ks.test(vafEp[vafEp>fmin & vafEp<fmax], vafEss[vafEss>fmin & vafEss<fmax])

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


