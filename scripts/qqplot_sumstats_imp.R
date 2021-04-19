#!/usr/bin/env Rscript
manhattan_choose_col <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
                                                                                           "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
                                  genomewideline = -log10(5e-08), otherlines=FALSE,highlight = NULL, logp = TRUE, 
                                  annotatePval = NULL, low_alpha_snps=NULL, annotateTop = NULL, hl_col = "green3", hl_alpha = 0.5, annotate = NULL, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d[!d$SNP  %in% low_alpha_snps,], points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i] & !d$SNP  %in% low_alpha_snps, ], points(pos, 
                                                                                    logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (otherlines)
    trasnpos <- function(x) {loglog_p * log10(x) / log10(loglog_p)}
  abline(h = trasnpos(c(-log10(5e-10),-log10(5e-20),-log10(5e-30),-log10(5e-40),-log10(5e-50),-log10(5e-60),-log10(5e-70),-log10(5e-80))), col = alpha("grey68", 0.5))
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = alpha(hl_col, hl_alpha), 
                             pch = 20,  ...))
  }
  if (!is.null(low_alpha_snps)) {
    d.low_alpha_snps = d[which(d$SNP %in% low_alpha_snps), ]
    with(d.low_alpha_snps, points(pos, logp, col = "red",#alpha("grey68", 0.1), 
                                  pch = 20, ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625,
                                                labs = topHits$SNP,
                                                cex = 0.5), ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  
  if (!is.null(annotate)) {
    d.annotate = d[which(d$SNP %in% annotate), ]
    par(xpd = TRUE)
    with(d.annotate, textxy(pos, -log10(P),
                            offset = 0.625,
                            labs = d.annotate$rsid,
                            cex = 1), ...)
  }
  
  par(xpd = FALSE)
}

# # # # # 
library(qqman)
library(optparse)
library(data.table)
library(R.utils)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c","--chrcol"), type="character", default="CHR",
              help="chromosome column [default= %default]", metavar="character"),
  make_option(c("-p","--pval_col"), type="character", default="P",
              help="pvalue column [default= %default]. This can be a comma separated list and plots will be generated for each of these", metavar="character"),
  make_option(c("-b","--bp_col"), type="character", default="BP",
              help="bp column [default= %default]", metavar="character"),
  make_option(c("-l","--loglog_pval"), type="integer", default=10,
              help="-log10 p-val threshold for using log-log scale in manhattan plot [default= %default]", metavar="integer"),
  make_option(c("-y","--loglog_ylim"), type="integer", default=324,
              help="-log10 p-val limit for y-axis of log-log manhattan [default= %default]", metavar="integer"),
  make_option(c("--mlog10p"), type="logical", default=FALSE,
              help="whether the p-values are -log10 or not [default= %default]", metavar="logical"),
  make_option(c("-m","--minrep_col"), type="character",
              help="if given then chr:bp:ref:alt identifier assumed and chr and bp are read from there [default= %default]", metavar="character"),
  make_option(c("-s","--snp_col"), type="character",
              help="if given then snp identifiers are read from there [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

file <- opt$options$file
print(paste("reading file:", file))

data <- fread(file, header=T)

options(bitmapType='cairo')

print(str(opt))
bp_col <- opt$options$bp_col
chr_col <- opt$options$chrcol

print(summary(data))
print( summary( data[[chr_col]] ) )
#colnames(data) <- toupper( colnames(data) )

pcols <- unlist(strsplit(opt$options$pval_col,","))

output_prefix=file
title_text <- sub(".txt.*", "", file)

if( !is.null(opt$options$out)) {
  output_prefix=opt$options$out
}


if(! is.null(opt$options$minrep_col ) ) {
  print("getting BP and CHR from minrepid")
  split <- strsplit(as.character(data[[opt$options$minrep_col]]), ":")
  data[[bp_col]] <- unlist( lapply( split, function(x) as.numeric(x[2]) ))
  data[[chr_col]] <- unlist( lapply( split, function(x) x[1] ))
}

print(append(pcols,c(bp_col,chr_col)))

if( any( ! append(pcols,c(bp_col,chr_col)) %in% colnames(data)   )) {
  stop( paste0("All required columns do not exist in the data: ", paste(pcols,sep=",", collapse=""),",", bp_col, ",",chr_col,  collapse="" ))
}


print(summary(as.factor(data[[chr_col]])))

data[[chr_col]] <- gsub("chr","",data[[chr_col]])
data[[chr_col]] <- gsub("X|chrX","23",data[[chr_col]])
data[[chr_col]] <- gsub("Y|chrY","24",data[[chr_col]])
data[[chr_col]] <- gsub("MT|chrMT|M|chrM","25",data[[chr_col]])

data[[chr_col]] <- as.numeric(data[[chr_col]])
data <- data[ !is.na(data[[chr_col]]) ]

quants <- c(0.7,0.5,0.1,0.01, 0.001)


for( pcol in pcols) {
  subdata <- data[ !is.na(data[[pcol]]) & is.numeric( data[[pcol]]  ) ]
  if (opt$options$mlog10p) {
    data[[pcol]] <- ifelse(10^-data[[pcol]] < 5e-324, 5e-324, 10^-data[[pcol]])
  }
  lambda  <- round(  quantile(  (qchisq(1-subdata[[pcol]], 1) ), probs=quants ) / qchisq(quants,1), 3)
  png( paste(output_prefix,"_", pcol ,"_qqplot.png", sep="" ))
  qq(subdata[[pcol]])
  title(c(title_text, paste("\n", "\nlambda ", quants, ": ", lambda, sep="" )), line=-0.2)
  dev.off()
  
  sink( paste(output_prefix,"_",  pcol ,"_qquantiles.txt", sep="" ) )
  cat( paste( quants, ":", lambda, sep=""))
  sink()
  
  print("subsetting p-vals < 0.01 for manhattan...")
  
  subdata <- subdata[ subdata[[pcol]]<0.01 & subdata[[pcol]]>0 ]
  print( paste0("Plotting manhattan with ", nrow(subdata), " variants") )
  print( summary(subdata[[pcol]] ))
  png( paste(output_prefix,"_",pcol,"_manhattan.png", sep=""), width=1000, height=400)
  logs <- -log10(subdata[[pcol]])
  sd <- subdata[subdata[["raiss.imputed"]]==1]
  snps = sd[[opt$options$snp_col]]
  n_snps = nrow(data)
  n_imp = nrow(data[data[["raiss.imputed"]]==1])
  tt = paste0(title_text, " - variants: ", n_snps, ", imputed: ", n_imp)
  manhattan( data.table(subdata[,c(bp_col,pcol,chr_col,opt$options$snp_col), with=F]) , chr=chr_col, bp=bp_col, p=pcol, ylim=c( 2,max(logs)+1), main=tt, suggestiveline=F, snp=opt$options$snp_col,
             highlight=snps, hl_col = "cyan4", hl_alpha = 0.4 )
  dev.off()
  
  
  print("plotting log-log manhattan")
  loglog_p <- opt$options$loglog_pval
  logs <- ifelse(logs < loglog_p, logs, loglog_p * log10(logs) / log10(loglog_p))
  loglog_ylim <- opt$options$loglog_ylim
  if (loglog_ylim==0) {
    max_loglog <- max(logs)
  } else {
    max_loglog <- loglog_p * log10(loglog_ylim) / log10(loglog_p)
  }
  subdata[["p_scaled"]] <- 10^(-logs)
  tick_pos <- round(seq(1, max_loglog, length.out=max_loglog))
  tick_lab <- sapply(tick_pos, function(pos) { round(ifelse(pos < loglog_p, pos, loglog_p^(pos/loglog_p))) })
  png( paste(output_prefix,"_",pcol,"_manhattan_loglog.png", sep=""), width=1000, height=400)
  manhattan( data.table(subdata[,c(bp_col,"p_scaled",chr_col,opt$options$snp_col), with=F]) , chr=chr_col, bp=bp_col, p="p_scaled", ylim=c( 2,max_loglog), yaxt="n", main=title_text, suggestiveline=FALSE, snp=opt$options$snp_col,
             highlight=snps, hl_col = "cyan4", hl_alpha = 0.4 )
  axis(2, at = tick_pos, labels=tick_lab, las=2)
  dev.off()
}
