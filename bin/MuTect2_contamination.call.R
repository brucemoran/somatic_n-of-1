#! R
argsIn <- commandArgs(trailingOnly = TRUE)
tab1 <- read.table(argsIn[1], header=T)
if(tab1[2] >= 1){
  write.table(tab1,
              file=paste0(argsIn[2], ".contamination.issue.table"),
              quote=F, row=F, col=T, sep="\t")
  print(paste0(argsIn[2], ": possible contamination"))
}

if(tab1[2] < 1){
  write.table(tab1,
              file=paste0(argsIn[2], ".contamination.no-issue.table"),
              quote=F, row=F, col=T, sep="\t")
}
