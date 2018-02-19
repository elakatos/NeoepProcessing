
dir <- '~/CRCdata/CRCmseq_Polyp'
summaryTable <- read.table(paste0(dir,'/Neopred_results/CRCmseq.neoantigens.summarytable.txt'), header=T, row.names=1)

barcolors = c('firebrick', 'darkorange3', 'goldenrod1')
pdf(paste0(dir, '_summary.pdf'), height = 5, width=8)
barplot(t(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')]/summaryTable$Total)), col=barcolors, legend=c('Clonal', 'Shared','Subclonal'), las=2)
barplot(t(as.matrix(summaryTable[,c('Total_WB', 'Total_SB')]/summaryTable$Total)), col=barcolors, legend=c('Weak binders', 'Strong binders'), las=2)
dev.off()

#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Generate summary table for MUTATIONS ------------------------------------


#aggregate(epTable$V2, by = list(epTable$V1), FUN=sum)
