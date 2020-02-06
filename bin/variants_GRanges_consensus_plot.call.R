##script template for vcf consensus variant calls under nextflow;
##this is the associated R script from sh template, named as per s/\.R/\.sh/
##NB this takes input of dir, variant type supplied as 'input'

#! /usr/bin/R

##inputs
##[1] <- functions file
##[2] <- GERMLINE = ${params.germline}
##[3] <- VEPVCFPATTERN = $vartype".pass.vep"
##[4] <- INCLUDEDORDER, if more than 1 sample define order for plots
##[5] <- NAMECALLERS, variant callers to use in superset, 2 and only 2 required!

options(stringAsFactors=FALSE)
args <- commandArgs(trailingOnly = TRUE)
##if defined, used R_LIB_PATHS
source(args[1])

##germline sample ID
GERMLINE <- args[2]

##VEP vcf pattern
##raw VCF must be *raw.vcf
VEPVCFPATTERN <- args[3]
RAWVCFPATTERN <- "raw.vcf"

##parse VCFs
##all should exist in current dir, all output to same dir
vcfVec <- dir(pattern=VEPVCFPATTERN)
vcfList <- vcfVec[grep(paste0(VEPVCFPATTERN,"$"),vcfVec)]
rawVec <- dir(pattern=RAWVCFPATTERN)
rawList <- rawVec[grep(paste0(RAWVCFPATTERN,"$"),rawVec)]
inputList <- list(vcfList, rawList)
vcfExt <- gsub("-","_",
               paste(strSplitVec(
                              gsub("\\.vcf","",
                              vcfList[[1]]),
                              "\\.")[-c(1,2),1],
                    collapse="."))

##create order from vcfList
if(length(args)==4){
  if(length(args[4][grep("\\,", args[4])])>0){
    INCLUDEDORDER <- strSplitVec(args[4],",")[,1] ##comma delim string
  }
  else{
    INCLUDEDORDER <- vcfList %>%
                     lapply(.,function(f){strsplit(f,"\\.")[[1]][1]}) %>%
                     unlist() %>% unique()
  }
}

##operate over varianttype, raw
for(x in 1:2){
  inList <- inputList[[x]]
  samples <- strSplitVec(inList, "\\.")[1,]
  callers <- strSplitVec(inList, "\\.")[2,]
  outExt <- gsub("-","_",
                 paste(strSplitVec(
                                gsub("\\.vcf","",
                                inList[[1]]),
                                "\\.")[-c(1,2),1],
                      collapse="."))
  print(paste0("Working on: ", outExt))

  dirtest <- dir(pattern=paste0(outExt,".RData"))
  if(length(dirtest)==0){

    ##run function to make list of GRanges per caller
    grList <- base::as.list(unique(callers))
    names(grList) <- unique(callers)

    for (callern in 1:length(unique(callers))){
      caller <- unique(callers)[callern]
      print(paste0("Caller: ",caller))
      ##parse VCFs, raw or VEP annotated
      if(outExt == "raw"){
        grList[[caller]] <- lapply(inList[grep(caller,callers)],function(f){
          suppressWarnings(vcfParseGR(f, GERMLINE))
        })
      }
      else{
        grList[[caller]] <- lapply(inList[grep(caller,callers)],function(f){
          suppressWarnings(vcfVepAnnParseGR(f, GERMLINE))
        })
      }
      names(grList[[caller]]) <- samples[grep(caller,callers)]
    }

    dirtest <- dir(pattern=paste0(outExt,".RData"))
    ##assign output, save to dir
    assignedName <- paste0(outExt)
    assign(assignedName, value=grList)
    saveFile <- paste0(outExt,".RData")
    save(list=assignedName, file=saveFile)
  }
  else{
    print(paste0("Found: ", outExt, ".RData, this is used"))
  }
}

##load set of RData GRanges, not raw
print("Loading GRanges")
dirIn <- dir(pattern=".RData")
#dirIn <- dirIn[grep("raw",dirIn,invert=T)]
loadedGR <- lapply(dirIn,function(x){
  load(x,envir=.GlobalEnv)
  })
names(loadedGR) <- unlist(loadedGR)

##set GRnages into lists
rawList <- get(loadedGR$raw)
varList <- get(loadedGR[[vcfExt]])

##callers, samples used
callers <- names(varList)
NAMECALLERS <- c(callers[1], callers[2])
samples <- names(varList[[1]])

##ensure some variants in each sample to check on
anyVars <- c()
for(x in 1:length(varList)){
  for(y in 1:length(samples))
  anyVars <- c(anyVars, length(varList[[x]][[y]]))
}

if(all(anyVars > 0)){
  ##get GRanges superset per mutype
  GRsuper <- GRsuperSet(varList, nameCallers=NAMECALLERS)

  ##get list to plot from
  plotList <- atLeastTwo(varList,
                         GRsuper,
                         taga="HM")

  ##plot
  if(!is.list(plotList)){
    plotConsensusSingle(plotList, rawList, taga="HM")
  }
  if(is.list(plotList)){
    if(length(plotList)>1){
      if(length(plotList[[1]])!=0 & length(plotList[[2]])!=0){
        plotConsensusList(plotList, rawList, taga="HM", includedOrder=INCLUDEDORDER)
      }
      if(length(plotList[[1]])==0 & length(plotList[[2]])==0){
        print("No variants returned at HM, support across callers lacking")
      }
    }
  }

  ##run to get all impacts, print but not plot
  ##get GRanges superset per mutype
  GRsuperALL <- GRsuperSet(varList, impacts=c("HIGH","MODERATE","MODIFIER","LOW"), nameCallers=NAMECALLERS)

  ##get list to plot from
  plotListAll <- atLeastTwo(varList, GRsuperALL, taga="ALL")

  ##plot
  if(is.list(plotListAll)){
    if(length(plotListAll[[1]])!=0 & length(plotListAll[[2]])!=0){
      plotConsensusList(plotListAll, rawList, taga="ALL", includedOrder=INCLUDEDORDER)
    }
    if(length(plotListAll[[1]])==0 & length(plotListAll[[2]])==0){
      print("No variants returned at ALL, support across callers lacking")
    }
  }
  if(!is.list(plotListAll)){
    plotConsensusSingle(plotListAll, rawList, taga="ALL")
  }
}
if(all(anyVars == 0)){
  print(paste0("No variants found for: ", args[3], ", exited..."))
}
