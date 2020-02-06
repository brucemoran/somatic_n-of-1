#!/usr/bin/env R

##adapted from https://gitlab.com/pesk/glioma.nano-seq/blob/master/scripts/createCNprofile.R
libs <- c("QDNAseq", "QDNAseq.hg19", "PSCBS", "matrixStats", "ggplot2", "caTools", "dplyr")
libsLoaded <- lapply(libs,function(lib){suppressMessages(library(lib, character.only = TRUE))})
options(scipen=999)
argsIn <- commandArgs(trailingOnly = TRUE)
tumourbam <- argsIn[1]
bin.size <- argsIn[2]
bed <- gsub("bam", paste0("cna.", bin.size,"kb.bed"), tumourbam)
pdf.file <- gsub("bam", paste0("cna.", bin.size,"kb.pdf"), tumourbam)
rplot.file <- gsub("bam", paste0("cna.", bin.size,"kb.rplot"), tumourbam)
alpha <- 0.001

##plot function
CNplot <- function(r) {
  cn <- data.frame(chr = as.character(r@featureData@data$chromosome), start = r@featureData@data$start, end = r@featureData@data$end, counts = as.vector(r@assayData$copynumber), segmented = as.vector(r@assayData$segmented), use = r@featureData@data$use)
  cn$chr <- factor(cn$chr, levels = c(1:22,"X","Y"))
  cn <- subset(cn, use)
  transf <- function(x) log2(x + 1) - 1
  cn$counts <- transf(cn$counts)
  cn$segmented <- transf(cn$segmented)
  cn$movavg <- runmean(cn$counts, k=50)
  cn <- cn %>%
    group_by(chr) %>%
    arrange(start) %>%
    mutate(movavg = runmean(counts, k=50))
  ggplot(cn, aes(x = start, y = counts)) +
    geom_point(colour = "grey", size = 1) +
    facet_grid(~ chr, scales = "free_x", space = "free_x") +
    geom_line(aes(y = movavg), colour="#e41a1c") +
    geom_line(aes(y = segmented), color="#377eb8", size=1) +
    theme_classic() +
    theme(axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.minor.y=element_line(colour="black",size=0.5),
          panel.spacing.x=unit(0, "lines"),
          strip.text.x = element_text(size = 7),
          strip.background = element_blank(),
          panel.background = element_rect(colour = "grey", size=0.4)) +
    scale_y_continuous(minor_breaks = 0, limits = c(-1,1))
}

##read inputs
bins <- getBinAnnotations(binSize = as.numeric(bin.size))
r <- binReadCounts(bins, bamfiles = tumourbam, minMapq=20, isSecondaryAlignment=FALSE)
r <- applyFilters(r, residual = F, blacklist = F, mappability = F, bases = 100, chromosomes = "Y")
r <- estimateCorrection(r)
r <- correctBins(r, method = "none")
r <- normalizeBins(r, method = "mean")
r <- compareToReference(r, c(2, FALSE), force=T)
r <- normalizeBins(r, method = "mean", force=T)

x <- data.frame(r@featureData@data$chromosome, r@featureData@data$start, r@assayData$copynumber, r@featureData@data$use)
xx <- subset(x, r.featureData.data.use==T)
xx$r.featureData.data.use <- NULL
colnames(xx) <- c("chromosome","x","y")
xx$chromosome <- unlist(lapply(rownames(xx),function(f){strsplit(f,":")[[1]][1]}))
xx$chromosome <- gsub("X",23,xx$chromosome)
xx$chromosome <- gsub("Y",24,xx$chromosome)

gaps <- findLargeGaps(xx, minLength = 5e+06)
knownSegments <- gapsToSegments(gaps)
fit <- segmentByCBS(xx, knownSegments = knownSegments, alpha = as.numeric(alpha))

r <- segmentBins(r, transformFun="sqrt", alpha = as.numeric(alpha), min.width=2, undo.splits="none")
r <- callBins(r, method="cutoff")

save(r, file = rplot.file)
CNplot(r)
ggsave(pdf.file, width=225, height=45, units = "mm")

### simplify segments
segments=data.frame(chr = as.character(r@featureData@data$chromosome), start = r@featureData@data$start, end = r@featureData@data$end, seg.mean=r@assayData$segmented[,1], call=r@assayData$calls[,1])

shortSeg <- segments %>% arrange(chr,start) %>%
    mutate(identical=as.integer(seg.mean != lag(seg.mean))) %>%
    mutate(identical=ifelse(is.na(identical),1,identical)) %>%
    mutate(segmentNo = cumsum(identical)) %>%
    group_by(segmentNo) %>%
    summarise(chr=first(chr), start=min(start), end=max(end), seg.mean=log2(max(seg.mean)),call=first(call)) %>%
    select(-segmentNo)

# write .seg file for IGV
write.table(data.frame(sample=sub(".*(\\d{4}T).*", "\\1", tumourbam), shortSeg), file=bed, quote=F, row.names = F, sep="\t")
