#! /usr/bin/R

##functions for SNV_consensus generation and plotting

##load libraries
libs <- c("customProDB", "ensemblVEP", "org.Hs.eg.db", "GenomicRanges", "tidyverse", "plyr", "pheatmap", "data.table")

libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

##parse VCF into GR, for use on un-annotated VCFs
vcfParseGR <- function(vcfIn, germline){

  vcf <- readVcf(vcfIn)
  gr <- suppressWarnings(InputVcf(vcfIn))

  ##parse info
  infor <- info(header(vcf))

  ##somatic
  if(!is.null(germline)){
    somName <- names(gr)[names(gr)!=germline]
  }
  if(is.null(germline)){
    somName <- names(gr)
  }
  print(paste0("Working on: ",somName))
  som <- gr[[somName]]
  ##ensure an AF is there, pisces has VF instead (thanks pisces dev=D)
  if(! "AF" %in% names(mcols(som))) {
    AD <- as.numeric(unlist(mcols(som)["AD"]))
    AD1 <- as.numeric(unlist(mcols(som)["AD.1"]))
    tot <- AD+AD1
    mcols(som)$AF <- AD1/tot
  }
  seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]
  return(som)
}

#parse VEP annotated VCF into GR, uses CANONICAL transcript
vcfVepAnnParseGR <- function(vcfIn, germline){

  vcf <- readVcf(vcfIn)
  if(dim(vcf)[1] != 0){
    gr <- suppressWarnings(InputVcf(vcfIn))

    ##parse info
    infor <- info(header(vcf))

    ##VEP annotation naming
    annNames <- unlist(strSplitFun(infor[rownames(infor)=="ANN",]$Description,"\\|"))

    ##somatic
    if(!is.null(germline)){
      somName <- names(gr)[names(gr)!=germline]
    }
    if(is.null(germline)){
      somName <- names(gr)
    }
    print(paste0("Working on: ",somName))
    som <- gr[[somName]]
    seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]

    ##ensure an AF is there, pisces has VF instead (thanks pisces dev=D)
    if(! "AF" %in% names(mcols(som))){
      AD <- as.numeric(unlist(mcols(som)["AD"]))
      AD1 <- as.numeric(unlist(mcols(som)["AD.1"]))
      tot <- AD+AD1
      mcols(som)$AF <- AD1/tot
    }

    ##annotation by CANONICAL, and add to mcols
    somAnnDf <- t(as.data.frame(lapply(strSplitFun(som$ANN,"\\|"),function(ff){
      if(ff[annNames=="CANONICAL"]=="YES"){
        if(is.null(ff)){ff<-rep("",length(annNames))}
        if(length(ff)!=length(annNames)){
          lengExtra <- length(annNames)-length(ff)
          ff<-c(ff,rep("",lengExtra))}
        return(ff)}
        else{
          return(rep("",length(annNames)))
        }
      })))
    colnames(somAnnDf) <- annNames

    if(sum(dim(somAnnDf)) != 0){
      values(som) <- cbind(as.data.frame(mcols(som)),somAnnDf)
      som$ANN <- NULL
    }
    som <-unique(som)

    return(som)
  }
  else{
    print("No variants found")
    return(GRanges())
  }
}

##create single-letter HGVS protein annotation (VEP outputs 3-letter)
##take vector, gsub out aa3 for aa1
subHGVSp <- function(inVec){

  aa1 <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X", "D", "R", "C", "C", "C", "C", "C", "C", "C", "C", "H", "G", "H", "H", "H", "H", "H", "H", "D", "K", "K", "M", "K", "M", "C", "F", "Y", "S", "T")
  ##amino acid 3 letter to gsub HGVSp
  aa3 <- c("Ala","Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Aba", "Ash", "Cir", "Cme", "Cmt", "Csd", "Cso", "Csw", "Csx", "Cym", "Cyx", "Dde", "Glh", "Hid", "Hie", "Hip", "Hsd", "Hse", "Hsp", "Ias", "Kcx", "Lyn", "Mho", "Mly", "Mse", "Ocs", "Pff", "Ptr", "Sep", "Tpo")

  ##include * for Ter
  aa1 <-c(aa1,"*")
  aa3 <- c(aa3, "Ter")

  unlist(lapply(inVec,function(f){
    #check matches (should be none or two)
    a3 <- aa3[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    a1 <- aa1[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    ##beauty:
    #https://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
    if(length(a3)>0){
      names(a1) <- a3
      str_replace_all(f,a1)
    }
    else{
      return("")
    }
  }))
}

##create GRsupersets from list of variants
#varList = nested list of [[caller]][[samples1..n]]
##NB hardcode via =NULL
GRsuperSet <- function(varList, impacts=NULL, mcolsWant=NULL, nameCallers=NULL){

  ##set NULL vars
  if(is.null(impacts)){impacts <- c("HIGH","MODERATE")}
  if(is.null(mcolsWant)){mcolsWant <- c("AD", "AD.1", "AF", "Consequence", "IMPACT", "SYMBOL", "HGVSc", "HGVSp", "CLIN_SIG")}

  if(length(nameCallers) != 2){
    print("Require only 2 callers, no more and no less!")
    exit;
  }

  ##set up output
  GRsuper <- as.list(names(varList[[1]]))
  callerNames <- names(varList)

  ##read first entry
  call1 <- varList[[nameCallers[1]]]
  call2 <- varList[[nameCallers[2]]]

  #exclude MT, GL
  seqwant <- c(seq(from=1,to=22,by=1), "X")

  ##iterate over samples in callerset
  if(length(call1) > 1){
    for (x in seq_along(call1)){
      print(paste0("Working on: ",names(call1)[x]))

      calls1 <- call1[[x]]
      calls2 <- call2[[x]]

      ##test all wanted mcols exist, rename if "VF" not "AF" (Pisces)
      for(y in 1:2){
        if(length(mcolsWant[mcolsWant %in% names(mcols(call1[[y]]))]) != length(mcolsWant)){
          gsub("VF","AF", names(mcols(call1[[y]])))
        }
        if(length(mcolsWant[mcolsWant %in% names(mcols(call2[[y]]))]) != length(mcolsWant)){
          gsub("VF","AF", names(mcols(call2[[y]])))
        }
      }
      calls1$HGVSp1 <- subHGVSp(calls1$HGVSp)
      calls2$HGVSp1 <- subHGVSp(calls2$HGVSp)

      ##sets of call1, 2 and the difference
      gr11 <- calls1[calls1$IMPACT %in% impacts, names(mcols(calls1)) %in% mcolsWant]
      gr22 <- calls2[calls2$IMPACT %in% impacts, names(mcols(calls2)) %in% mcolsWant]
      gr12 <- suppressWarnings(setdiff(gr22, gr11))

      if(length(gr11)!=0 & length(gr12)!=0){
        mcols(gr11) <- mcols(gr11)[,mcolsWant]
        mcols(gr12) <- mcols(gr12)[,mcolsWant]
      }

      ##add 1 and difference (the superset of variants)
      GRsuper[[x]] <- suppressWarnings(c(gr11,gr12))
    }
  }
  else{
    ##only one sample, return this
    calls1 <- call1[[1]]
    calls2 <- call2[[1]]

    ##test all wanted mcols exist, rename if "VF" not "AF" (Pisces)
    if(length(mcolsWant[mcolsWant %in% names(mcols(call1[[1]]))]) != length(mcolsWant)){
        gsub("VF","AF", names(mcols(call1[[1]])))
    }
    if(length(mcolsWant[mcolsWant %in% names(mcols(call2[[1]]))]) != length(mcolsWant)){
        gsub("VF","AF", names(mcols(call2[[1]])))
    }
    calls1$HGVSp1 <- subHGVSp(calls1$HGVSp)
    calls2$HGVSp1 <- subHGVSp(calls2$HGVSp)

    ##sets of call1, 2 and the difference
    gr11 <- calls1[calls1$IMPACT %in% impacts, names(mcols(calls1)) %in% mcolsWant]
    gr22 <- calls2[calls2$IMPACT %in% impacts, names(mcols(calls2)) %in% mcolsWant]
    gr12 <- suppressWarnings(setdiff(gr22, gr11))

    if(length(gr11)!=0 & length(gr12)!=0){
      mcols(gr11) <- mcols(gr11)[,mcolsWant]
      mcols(gr12) <- mcols(gr12)[,mcolsWant]
    }

    ##add 1 and difference (the superset of variants)
    GRsuper[[1]] <- suppressWarnings(c(gr11,gr12))
  }
  return(GRsuper)
}

##find consensus in at least two callers, therefore in GRsuper
##this produces a set of positions to plot
atLeastTwo <- function(varList, GRsuper, taga, tmb=NULL){

  ##set TMB
  if(is.null(tmb)){
    tmb <- "null"
  }

  ##set vars
  callers <- names(varList)
  samples <- names(varList[[1]])

  ##iterate over list of callers
  if(length(samples) > 1){
    GRplots <- lapply(seq_along(samples), function(x){

        sample <- samples[x]
        print(paste0("Working on: ",sample))

        ##all possible combinations of intersects of callers
        ##output to clean GR
        upl <- GRanges()
        upl1 <- apply(t(combn(length(callers), m=2)), 1, function(xx){
          print(paste(callers[xx[1]]," vs. ",callers[xx[2]]))
          gr1 <- varList[[names(varList)[xx[1]]]][[x]]
          gr2 <- varList[[names(varList)[xx[2]]]][[x]]
          gri <- suppressWarnings(GenomicRanges::intersect(gr1,gr2))
          upl <- suppressWarnings(c(upl, gri))
        })
        upl2 <- upl1[[1]]
        if(length(upl1) > 1){
          for(xx in 2:length(upl1)){
            upl2 <- suppressWarnings(c(upl2, upl1[[xx]]))
          }
          GRplot <- suppressWarnings(GRsuper[[x]][GRsuper[[x]] %in% unique(upl2)])
        }
        if(length(upl1) == 1){
          GRplot <- intersect(GRsuper[[x]], unique(upl2))
        }

        #TMB
        fileOut <- paste0(sample, ".", taga, ".consensus.tab")
        vcfOut <- paste0(sample, ".", taga, ".pcgr.all.tab.vcf")

        if(tmb == "snv"){
            tmbOut <- exomeTumourMutationBurden(GRplot)
            fileOut <- paste0(sample, ".", taga, ".TMB_", tmbOut, "_SNV-Mb.consensus.tab")
        }
        write.table(GRplot, file=fileOut, quote=F, row=F, col=T, sep="\t")

        ##VCF output
        ###CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
        ##1	33238563	.	G	A	.	PASS	TVAF=0.145833333333;TDP=49;	GT:DPC:DPT:ADC:ADT	0/1:.:DP:.,.:AD,AD.1
        vcfGRplot <- as_tibble(as.data.frame(GRplot), rownames="r") %>%
          tidyr::separate(r,
                          c("#CHROM", "POS", "REF", "ALT"),
                          sep = "[:_\\/]") %>%
          dplyr::rename(AD1='AD.1') %>%
          dplyr::mutate(POS = as.integer(POS),
                        ID = ".",
                        QUAL = ".",
                        INFO = ".",
                        FILTER = "PASS",
                        FORMAT = "GT:DPC:DPT:ADC:ADT",
                        ADSUM = as.numeric(AD) + as.numeric(AD1),
                        sampleID = paste0("0/1:.:",
                                                      ADSUM,
                                                      ":.,.:",
                                                      AD,",",AD1)) %>%
          dplyr::select('#CHROM', POS, ID, REF, ALT, QUAL, FILTER, INFO,  FORMAT, sampleID) %>%
          dplyr::rename(!!sample:="sampleID")

        write_tsv(path=vcfOut, vcfGRplot)
        return(GRplot)
      })
      names(GRplots) <- samples
    }
    else{
      ##only one sample
      sample <- samples[1]
      print(paste0("Only one sample: ",sample))

      ##all possible combinations of intersects of callers
      ##output to clean GR
      upl <- GRanges()
      upl1 <- apply(t(combn(length(callers), m=2)), 1, function(xx){
        print(paste(callers[xx[1]]," vs. ",callers[xx[2]]))
        gr1 <- varList[[names(varList)[xx[1]]]][[1]]
        gr2 <- varList[[names(varList)[xx[2]]]][[1]]
        gri <- suppressWarnings(GenomicRanges::intersect(gr1,gr2))
        upl <- suppressWarnings(c(upl, gri))
      })
      upl2 <- upl1[[1]]
      if(length(upl1) > 1){
        for(xx in 2:length(upl1)){
          upl2 <- suppressWarnings(c(upl2, upl1[[xx]]))
        }
        GRplots <- suppressWarnings(GRsuper[[1]][GRsuper[[1]] %in% unique(upl2)])
      }
      if(length(upl1) == 1){
        GRplots <- intersect(GRsuper[[1]], unique(upl2))
      }

      #TMB
      fileOut <- paste0(sample,".",taga,".consensus.tab")
      vcfOut <- paste0(sample,".",taga,".pcgr.all.tab.vcf")

      if(tmb == "snv"){
          tmbOut <- exomeTumourMutationBurden(GRplot)
          fileOut <- paste0(sample,".",taga,".TMB_",tmbOut,"_SNV-Mb.consensus.tab")
      }
      write.table(GRplots, file=fileOut, quote=F, row=F, col=T, sep="\t")

      ##VCF output
      ###CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
      ##1	33238563	.	G	A	.	PASS	TVAF=0.145833333333;TDP=49;	GT:DPC:DPT:ADC:ADT	0/1:.:DP:.,.:AD,AD.1
      vcfGRplots <- as_tibble(as.data.frame(GRplots), rownames="r") %>%
        tidyr::separate(r,
                        c("#CHROM", "POS", "REF", "ALT"),
                        sep = "[:_\\/]") %>%
        dplyr::rename(AD1='AD.1') %>%
        dplyr::mutate(POS = as.integer(POS),
                      ID = ".",
                      QUAL = ".",
                      INFO = ".",
                      FILTER = "PASS",
                      FORMAT = "GT:DPC:DPT:ADC:ADT",
                      ADSUM = as.numeric(AD) + as.numeric(AD1),
                      sampleID = paste0("0/1:.:",
                                                    ADSUM,
                                                    ":.,.:",
                                                    AD,",",AD1)) %>%
        dplyr::select('#CHROM', POS, ID, REF, ALT, QUAL, FILTER, INFO,  FORMAT, sampleID) %>%
        dplyr::rename(!!sample:="sampleID")

      write_tsv(path=vcfOut, vcfGRplots)
    }
  return(GRplots)
}

##function to return overlapping variants for all samples
##also write VCF for PCGR based on example BRCA VCF from Github
writeConsensusAll <- function(plotList, plotDF, taga, includedOrder=NULL, cons=NULL){

  if(is.null(cons)){
    print("Require 'cons=c(\"consensus\", \"not_cons\")' to determine output type...")
  }

  if(!is.null(cons)){
    ##input names
    samples <- gsub("-","_",names(plotList))

    if(is.null(includedOrder)){
      includedOrder <- samples
    }

    ##split on colon, paste into label for plotList
    labelsOut <- unlist(lapply(strSplitFun(rownames(plotDF), ":"), function(f){
      paste0(f[1],":",f[2])
    }))

    ##parse required from plotList
    plotOut <- lapply(plotList, function(f){
      f[names(f) %in% labelsOut]
    })

    ##test any output
    outGR <- suppressWarnings(unique(do.call("c", unname(plotOut))))
    seqlevels(outGR) <- sort(seqlevels(outGR))
    outGR <- sort(outGR)
    outDF <- as.data.frame(outGR)

    if(dim(outDF)[1] != 0){
      ##write output
      fileOut <- paste0(taga,".", cons, ".tab")
      write.table(outDF, file=fileOut, quote=F, row=F, col=T, sep="\t")
    }
  }
}

##create two plots: all consensus, those in 2+ samples
plotConsensusList <- function(plotList, rawList, taga, includedOrder=NULL){

  ##remove hyphens
  samples <- gsub("-","_",names(plotList))
  if(is.null(includedOrder)){
    includedOrder <- samples
  }
  if(!is.null(includedOrder)){
    includedOrder <- gsub("-","_",includedOrder)
  }

  ##combined set of all samples
  combGR <- suppressWarnings(unique(do.call("c", unname(plotList))))
  seqlevels(combGR) <- sort(seqlevels(combGR))
  combGR <- sort(combGR)
  combDF <- as.data.frame(combGR)

  ##labels for plot
  hgvsp <- unlist(lapply(combGR$HGVSp1,function(f){
    strsplit(f,"\\.")[[1]][3]
  }))
  uniqLabels <- paste(names(combGR),
                     combDF$SYMBOL,
                     combDF$Consequence,
                     hgvsp, sep=":")

  ##take those positions, then query raw calls
  ##allows 'rescue' of those falling out from arbitrary filters
  ##enough support previously to allow re-entry
  rawAFs <- function(rawList){
    as.data.frame(lapply(rawList,function(ff){
      afs <- rep(0,length(combGR))
      lapply(ff,function(fff){
        seqlevels(fff) <- sort(seqlevels(fff))
        ffs <- sort(fff)
        ffsi <- ffs[ffs %in% combGR]
        afs[combGR %in% ffsi] <- as.numeric(mcols(ffsi)$AF)
        return(afs)
      })
    }))
  }
  plotDFrawout <- rawAFs(rawList)
  if(length(warnings())>1){
    plotDFrawout <- rawAFs(rawList)
  }
  colnames(plotDFrawout) <- unlist(lapply(names(rawList),function(f){
    paste(f,samples,sep=".")
  }))
  plotDFrawout <- do.call(cbind,lapply(samples, function(ss){
    apply(plotDFrawout[,grep(ss, colnames(plotDFrawout))],1,max)
  }))

  plotDFrawout <- as.data.frame(plotDFrawout)
  colnames(plotDFrawout) <- samples
  rownames(plotDFrawout) <- uniqLabels

  ##which samples to include, and order
  ##remove those with all 0 frequency
  plotDForder <- plotDFrawout[rowSums(plotDFrawout)!=0, colnames(plotDFrawout) %in% includedOrder]

  ##reduce any frequency >50% to 50% (somatic should not be >50%)
  ##and plot is purely representative
  plotDForder[plotDForder > 0.5] <- 0.5

  ##find all 0, count to allow separation, order
  notoDF <- do.call(cbind,list(apply(plotDForder,2,function(f){
    f!=0
  })))
  notoVec <- apply(notoDF,1,function(f){table(f)[[1]]})
  notoVec <- notoVec[!is.na(notoVec)]
  notoVec1 <- notoVec[notoVec==1]
  notoVec2 <- notoVec[notoVec>1]

  ##those in all samples
  plotDF1 <- plotDForder[!rownames(plotDForder) %in% names(notoVec1), includedOrder]
  plotDF2 <- plotDForder[!rownames(plotDForder) %in% names(notoVec2), includedOrder]

  ##write out those in plotDF
  writeConsensusAll(plotList=plotList, plotDF=plotDF1, taga=taga, includedOrder=includedOrder, cons="consensus")
  writeConsensusAll(plotList=plotList, plotDF=plotDF2, taga=taga, includedOrder=includedOrder, cons="all")

  ##ordering
  plotDF1$labels <- rownames(plotDF1)
  orderDF1 <- plotDF1[, includedOrder]
  rownames(orderDF1) <- plotDF1$labels
  plotDF1ordered <- orderDF1
  orderDF2 <- do.call(order, c(data.frame(plotDF2[, 1:(ncol(plotDF2)-1)], plotDF2[, ncol(plotDF2)])))

  plotVec <- c()
  plotTag <- c()

  if(dim(orderDF1)[1] == 0){
    print("No shared variants found...")

    ##no shared variants
    if(length(orderDF2)[1] != 0){
      plotDF2ordered <- plotDF2[orderDF2,]
      plotDF2ordered <- rbind(plotDF2ordered)
      plotVec <- list(plotDF2ordered)
      plotTag <- list("only-private")
      print("Plotting private variants...")
    }

    if(length(orderDF2)[1] == 0){
      print("No variants to plot")
      quit("n")
    }
  }

  if(dim(orderDF1)[1] != 0){
      plotDF2ordered <- plotDF2[orderDF2,]
      plotDF2ordered <- rbind(plotDF2ordered, plotDF1ordered)
      plotVec <- list(plotDF1ordered, plotDF2ordered)
      plotTag <- list("only-shared", "private-shared")
      print("Plotting shared, private variants...")
  }
  for (pl in 1:length(plotVec)){
    if(dim(plotVec[[1]])[2]!=0){
      plotLabels <- rep("",dim(plotVec[[pl]])[1])
      row_fontsize <- 1
      colz <- colorRampPalette(c("lightgrey","dodgerblue","blue"))
      if(dim(plotVec[[pl]])[1] < 120){
        if(dim(plotVec[[pl]])[1]<20){row_fontsize=8}
        if(dim(plotVec[[pl]])[1]<50){row_fontsize=6}
        if(dim(plotVec[[pl]])[1]>50 & dim(plotVec[[pl]])[1]<100){row_fontsize=4}
        if(dim(plotVec[[pl]])[1]>100){row_fontsize=2}
        plotLabels <- rownames(plotVec[[pl]])
      }

      pdf(paste0(taga,".",plotTag[[pl]],".pdf"),onefile=F)
      pheatmap(plotVec[[pl]][,c(1:length(includedOrder))],
         breaks=seq(from=0,to=0.5,length.out=101),
         color=colz(100),
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         clustering_distance_rows=NA,
         cellwidth=12,
         legend=TRUE,
         fontsize_row=row_fontsize,
         labels_row=plotLabels,
         border_color="lightgrey",
         gaps_col=c(1:length(includedOrder))
      )
      dev.off()

    }
  }
}

plotConsensusSingle <- function(plotList, rawList, taga, includedOrder=NULL){

  ##remove hyphens
  sample <- names(plotList)[1]

  ##combined set of all samples
  combGR <- suppressWarnings(unique(plotList))
  seqlevels(combGR) <- sort(seqlevels(combGR))
  combGR <- sort(combGR)
  combDF <- as.data.frame(as.data.table(combGR))

  ##labels for plot
  hgvsp <- unlist(lapply(combGR$HGVSp1,function(f){
    strsplit(f,"\\.")[[1]][3]
  }))
  uniqLabels <- paste(names(combGR),
                     combDF$SYMBOL,
                     combDF$Consequence,
                     hgvsp, sep=":")

  ##take those positions, then query raw calls
  ##allows 'rescue' of those falling out from arbitrary filters
  ##enough support previously to allow re-entry
  rawAFs <- function(rawList){
    as.data.frame(lapply(rawList,function(ff){
      afs <- rep(0,length(combGR))
      lapply(ff,function(fff){
        seqlevels(fff) <- sort(seqlevels(fff))
        ffs <- sort(fff)
        ffsi <- ffs[ffs %in% combGR]
        afs[combGR %in% ffsi] <- as.numeric(mcols(ffsi)$AF)
        return(afs)
      })
    }))
  }
  plotDFrawout <- rawAFs(rawList)
  if(length(warnings())>=1){
    plotDFrawout <- rawAFs(rawList)
  }
  colnames(plotDFrawout) <- unlist(lapply(names(rawList),function(f){
    paste(f,sample,sep=".")
  }))
  plotDFrawout <- do.call(cbind,lapply(sample, function(ss){
    apply(plotDFrawout[,grep(ss, colnames(plotDFrawout))],1,max)
  }))

  plotDFrawout <- as.data.frame(plotDFrawout)
  colnames(plotDFrawout) <- sample
  rownames(plotDFrawout) <- uniqLabels

  ##ordering
  plotDFordered <- as_tibble(plotDFrawout, rownames="row") %>%
                   dplyr::arrange(.[[2]]) %>%
                   base::as.data.frame()

  if(dim(plotDFordered)[1] != 0){
    plotVec <- data.frame(row.names=plotDFordered[,1],
                          tag=plotDFordered[,2])
    plotTag <- "variants"

    plotLabels <- rep("",times=dim(plotVec)[1])
    row_fontsize <- 1
    colz <- colorRampPalette(c("lightgrey","dodgerblue","blue"))
    if(dim(plotVec)[1] < 120){
      if(dim(plotVec)[1]<20){row_fontsize=8}
      if(dim(plotVec)[1]<50){row_fontsize=6}
      if(dim(plotVec)[1]>50 & dim(plotVec)[1]<100){row_fontsize=4}
      if(dim(plotVec)[1]>100){row_fontsize=2}
      plotLabels <- rownames(plotVec)
    }

    pdf(paste0(taga,".",plotTag,".consensus.onlyOverlap.pdf"),onefile=F)
    pheatmap(plotVec,
       breaks=seq(from=0,to=0.5,length.out=101),
       color=colz(100),
       cluster_rows=FALSE,
       cluster_cols=FALSE,
       clustering_distance_rows=NA,
       cellwidth=12,
       legend=TRUE,
       fontsize_row=row_fontsize,
       labels_row=plotLabels,
       labels_col=taga,
       border_color="lightgrey"
    )
    dev.off()
  }
}

exomeTumourMutationBurden <- function(GRplot){

  ##get exome for Illumina Nextera Rapid
  exomeBed <- fread("https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed", showProgress=FALSE, data.table=FALSE)

  ##triage chr tag, add header
  exomeBed[,1] <- gsub("chr","",exomeBed[,1])
  colnames(exomeBed) <- c("seqname", "start", "end")
  exomeGR <- makeGRangesFromDataFrame(exomeBed, ignore.strand=TRUE)
  exomeSize <- sum(width(exomeGR))/1000000

  ##overlap with input
  hits <- as.data.frame(findOverlaps(GRplot, exomeGR))
  ##output
  GRplotExome <- GRplot[hits$queryHits]
  TMB <- round(length(GRplotExome)/exomeSize,digits=1)
  return(TMB)
}
