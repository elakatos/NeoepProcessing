setwd('~/CRCdata/Simulation_results/')

dirList = c('18_05_23_3', '18_05_23_4', '18_05_23_5', '18_05_23_6')

for (dir in dirList){
  
  for (i in 1:200){
    Npost <- read.table(paste0(dir,'/postIT_',i,'.txt'), header=T, sep=',')
    sparsedRows <- seq(1, nrow(Npost), round(nrow(Npost)/1000))
    Npost <- Npost[sparsedRows, ]
    write.table(Npost, file=paste0(dir,'/postIT_sparse_',i,'.txt'), sep=',', row.names=F)
  }
  
}