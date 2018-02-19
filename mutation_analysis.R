readAvinput <- function(sampleFile){
  avinput <- read.table(sampleFile, header=F)
  numtumors <- ncol(avinput)-18
  regionNames <- c('Normal', paste0(rep('Region',numtumors), 1:numtumors))
  names(avinput)[c(1,2,3,4,5, 18:(18+numtumors))] <- c('Chrom', 'Start', 'End', 'NormAll', 'AltAll', regionNames)
  return(avinput)
}




dir <- '~/CRCdata/CRCmseq_Set'

epTable <- read.table(paste0(dir, '/Neopred_results/CRCmseq.neoantigens.txt'), header=F)

sample = 'Set.01.snv'
sampleFile <- paste0(dir, '/avready',sample,'.avinput')
avinput <- readAvinput(sampleFile)
