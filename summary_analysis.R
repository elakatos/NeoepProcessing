
setwd('~/RNAseq/Neoepitopes/CRCmseq/')
fileName <- 'CRCmseq'
summaryTable <- read.table(paste0(fileName, '.neoantigens.summarytable.txt'), header=T, row.names=1)


barplot(t(as.matrix(summaryTable[,c('Clonal', 'Shared','Subclonal')]/summaryTable$Total)))
barplot(t(as.matrix(summaryTable[,c('Total_WB', 'Total_SB')]/summaryTable$Total)))

#Is the number of measured regions correlated with the ratio of clonal/subclonal epitopes?
cor.test(summaryTable$Subclonal/summaryTable$Total, rowSums(summaryTable[,1:16]>0))

#Is the number of epitopes detected correlated with the ratio?
cor.test(summaryTable$Shared/summaryTable$Total, log(summaryTable$Total))


# Generate summary table for MUTATIONS ------------------------------------

epTable <- read.table(paste0(fileName, '.neoantigens.txt'), header=F)



aggregate(epTable$V2, by = list(epTable$V1), FUN=sum)
