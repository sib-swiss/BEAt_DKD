#' Load file of experimental design
#'
#' @param filename Tab-delimited file formatted as
#' E.g.: file          cond      sample    ...
#'       mutant1.tsv   mutant    mutant1
#'       mutant2.tsv   mutant    mutant2
#'       wildtype1.tsv wildtype  wildtype1
#'       wildtype2.tsv wildtype  wildtype2
#' @param separator A chacracter to form contrast pair. By default, "-".
loadExpDesign <- function(filename) {
    exp.design <- read.table(filename, header=TRUE, sep="\t", colClasses="character", stringsAsFactors=F)
    rownames(exp.design) <- exp.design$sample
    #invisible(assert_that(length(grep("^[0-9]", unique(exp.design$cond))) == 0,
    #                      msg="Conditions starting by number are not allowed!"))
    return (exp.design)
}

#' Load file of contrasts to be studied
#'
#' @param filename Tab-delimited file formatted as
#' E.g.: cond    condref
#'       mutant1 wildtype
#'       mutant2 wildtype
#' @param separator A chacracter to form contrast pair. By default, "-".
loadContrasts <- function(filename, separator="-") {
    conts <- read.table(filename, header=TRUE, sep="\t", stringsAsFactors=F)
    invisible(assert_that(length(grep("^[0-9]", unique(c(t(conts))))) == 0,
                          msg="Conditions starting by number are not allowed!"))
    return (paste0(conts$cond, separator, conts$condref))
}

#' Load read counts
#' 
#' @param exp.design Data frame of experimental design.
#' @param file File name of count matrix. If file is provided, files reported in exp.design will be ignored. Default: NA.
#' @param header Logical value indicating whether there are headers in read count files from exp.design. Default: FALSE.
#' @return Matrix of read counts (gene x sample)
loadCountData <- function(exp.design, file = NA, header = !is.na(file), header.pattern = NA, header.replacement= NA) {
    if (is.na(file)) {
        rawCounts <- mclapply(exp.design$file, mc.cores=min(detectCores(), length(exp.design$file)), function(filename) {
            tab <- read.table(filename, header=header, sep='\t', row.names=1)
            tab <- tab[! rownames(tab) %in% c( "__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique"), , drop=F]
            return (tab)
        })
        assert_that(sum(unlist(mclapply(rawCounts, mc.cores=min(detectCores(), length(rawCounts)), function(x) {
            sum(rownames(x)!=rownames(rawCounts[[1]]))
        })))==0, msg="Genes in count files do not have the same order!")
        rawCounts <- do.call(cbind, rawCounts)
        colnames(rawCounts) <- exp.design$sample
    } else {
        rawCounts <- read.table(file, header=header, sep='\t', stringsAsFactors=F, row.names=1)
        colnames(rawCounts) <- gsub(colnames(rawCounts), pattern=header.pattern, replacement=header.replacement)
        ## reduce to conditions in exp.design
        rawCounts <- rawCounts[, grep(paste(paste0("_", exp.design$cond), collapse="|"), colnames(rawCounts), value=T)]
        batches <- sub(colnames(rawCounts), pattern=paste(paste0("_", exp.design$cond, "_.*$"), collapse="|"), replacement="")
        nobatches <- sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement="")
        date <- sub(colnames(rawCounts), pattern=".*_", replacement="")
        cond <- sub(sub(colnames(rawCounts), 
        				pattern=paste(paste0(batches, "_"), collapse="|"), replacement=""),
            		pattern=paste(paste0("_", date, "$"), collapse="|"), replacement="")
        pathology <- sub(cond, pattern=".*_", replacement="")
        cell <- sub(cond, pattern="_.*$", replacement="")
        ir <- sub(sub(cond, 
        			  pattern=paste(paste0("_", unique(pathology)), collapse="|"), replacement=""),
            	  pattern=paste(paste0(unique(cell), "_"), collapse="|"), replacement="")
        batch <- sapply(strsplit(batches, split="_"), function(x) x[1])
        batch.other <- sapply(strsplit(batches, split="_"), function(x) x[2])
        individual <- date
        colnames(rawCounts) <- nobatches
        exp.design <- data.frame(file=colnames(rawCounts),
                                 date=date,
                                 cond=cond,
                                 pathology=pathology,
                                 cell=cell,
                                 ir=ir,
                                 sample=colnames(rawCounts),
                                 run=colnames(rawCounts),
                                 batch=batch,
                                 batches=batches,
                                 batch.other=batch.other,
                                 individual=individual,
                                 row.names=colnames(rawCounts),
                                 stringsAsFactors=F)
        exp.design <- exp.design[with(exp.design, order(cell, cond)),]
 		rawCounts <- rawCounts[, exp.design$sample]
        # if (!is.null(exp.design$sample) && all(exp.design$sample %in% colnames(rawCounts))) {
            # rawCounts <- rawCounts[, exp.design$sample]
        # } else if (!is.null(exp.design$run) && all(exp.design$run %in% colnames(rawCounts))) {
            # rawCounts <- rawCounts[, exp.design$run]
        # } else if (sum(sapply(exp.design$sample, function(x) length(grep(x, colnames(rawCounts)))==1)) >= min(ncol(rawCounts), nrow(exp.design))) {
            # batches <- sub(colnames(rawCounts), pattern=paste(paste0("_", exp.design$sample), collapse="|"), replacement="")
            # nobatches <- sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement="")
            # exp.design <- exp.design[exp.design$sample %in% nobatches, ]
            # order.idx <- match(exp.design$sample, nobatches, nomatch=NA)
            # assert_that(all(exp.design$sample == nobatches[order.idx])) #TODO: do not accept NA now: some exp.design$sample not found in rawCounts
            # rawCounts <- rawCounts[, order.idx]
            # batches <- batches[order.idx]
            # nobatches <- nobatches[order.idx]
            # batch <- sapply(strsplit(batches, split="_"), function(x) x[1])
            # batch.other <- sapply(strsplit(batches, split="_"), function(x) x[2])
            # names(exp.design)[names(exp.design)=="batch"] <- "batch.old"
            # exp.design$batch <- batch
            # exp.design$batches <- batches
            # exp.design$batch.other <- batch.other
            # exp.design$sample <- exp.design$run <- rownames(exp.design) <- colnames(rawCounts)
        # } else if (sum(sapply(exp.design$run, function(x) length(grep(x, colnames(rawCounts)))==1)) >= min(ncol(rawCounts), nrow(exp.design))) {
            # batches <- sub(colnames(rawCounts), pattern=paste(paste0("_", exp.design$run), collapse="|"), replacement="")
            # nobatches <- sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement="")
            # exp.design <- exp.design[exp.design$run %in% nobatches, ]
            # order.idx <- match(exp.design$run, nobatches, nomatch=NA)
            # assert_that(all(exp.design$run == nobatches[order.idx])) #TODO: do not accept NA now: some exp.design$sample not found in rawCounts
            # rawCounts <- rawCounts[, order.idx]
            # batches <- batches[order.idx]
            # nobatches <- nobatches[order.idx]
            # batch <- sapply(strsplit(batches, split="_"), function(x) x[1])
            # batch.other <- sapply(strsplit(batches, split="_"), function(x) x[2])
            # names(exp.design)[names(exp.design)=="batch"] <- "batch.old"
            # exp.design$batch <- batch
            # exp.design$batches <- batches
            # exp.design$batch.other <- batch.other
            # exp.design$sample <- exp.design$run <- rownames(exp.design) <- colnames(rawCounts)
        # } else if (all(sapply(exp.design$sample, function(x) length(grep(x, colnames(rawCounts)))==1))) {
            # batches <- sub(colnames(rawCounts), pattern=paste(paste0("_", exp.design$sample), collapse="|"), replacement="")
            # nobatches <- sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement="")
            # order.idx <- sapply(exp.design$sample, function(x) which(nobatches==x))
            # assert_that(all(exp.design$sample == nobatches[order.idx]))
            # rawCounts <- rawCounts[, order.idx]
            # batches <- batches[order.idx]
            # nobatches <- nobatches[order.idx]
            # batch <- sapply(strsplit(batches, split="_"), function(x) x[1])
            # batch.other <- sapply(strsplit(batches, split="_"), function(x) x[2])
            # names(exp.design)[names(exp.design)=="batch"] <- "batch.old"
            # exp.design$batch <- batch
            # exp.design$batches <- batches
            # exp.design$batch.other <- batch.other
            # exp.design$sample <- exp.design$run <- rownames(exp.design) <- colnames(rawCounts)
        # } else if (all(sapply(exp.design$run, function(x) length(grep(x, colnames(rawCounts)))==1))) {
            # batches <- sub(colnames(rawCounts), pattern=paste(paste0("_", exp.design$run), collapse="|"), replacement="")
            # nobatches <- sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement="")
            # order.idx <- sapply(exp.design$run, function(x) which(nobatches==x))
            # assert_that(all(exp.design$run == nobatches[order.idx]))
            # rawCounts <- rawCounts[, order.idx]
            # batches <- batches[order.idx]
            # nobatches <- nobatches[order.idx]
            # batch <- sapply(strsplit(batches, split="_"), function(x) x[1])
            # batch.other <- sapply(strsplit(batches, split="_"), function(x) x[2])
            # names(exp.design)[names(exp.design)=="batch"] <- "batch.old"
            # exp.design$batch <- batch
            # exp.design$batches <- batches
            # exp.design$batch.other <- batch.other
            # exp.design$sample <- exp.design$run <- rownames(exp.design) <- colnames(rawCounts)
        # } else {
            # warning(paste0("Some sample and run in exp.design are not found in header of ", file, ". exp.design has been regenerated!"))
            # stop()
            # batches <- sub(colnames(rawCounts), pattern=paste(paste0("_", exp.design$cond, "_.*$"), collapse="|"), replacement="")
            # nobatches <- sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement="")
            # date <- sub(colnames(rawCounts), pattern=".*_", replacement="")
            # cond <- sub(sub(colnames(rawCounts), pattern=paste(paste0(batches, "_"), collapse="|"), replacement=""),
            				# pattern=paste(paste0("_", date, "$"), collapse="|"), replacement="")
            	# pathology <- sub(cond, pattern=".*_", replacement="")
            	# cell <- sub(cond, pattern="_.*$", replacement="")
            	# ir <- sub(sub(cond, pattern=paste(paste0("_", unique(pathology)), collapse="|"), replacement=""),
            			  # pattern=paste(paste0(unique(cell), "_"), collapse="|"), replacement="")
            	# batch <- sapply(strsplit(batches, split="_"), function(x) x[1])
            # batch.other <- sapply(strsplit(batches, split="_"), function(x) x[2])
            # exp.design <- data.frame(file=colnames(rawCounts),
            							 # date=date,
            							 # cond=cond,
            							 # pathology=pathology,
            							 # cell=cell,
            							 # ir=ir,
            						     # sample=colnames(rawCounts),           					
            							 # run=colnames(rawCounts),
            							 # batch=batch,
            							 # batches=batches,
            							 # batch.other=batch.other,           							 
            							 # row.names=colnames(rawCounts),
            							 # stringsAsFactors=F)
            	# exp.design <- exp.design[with(exp.design, order(cell, cond)),]
        # }
    }
    row2remove <- c("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual")
    rawCounts <- rawCounts[!str_detect(rownames(rawCounts), row2remove),]
    #rawCounts[is.na(rawCounts)] <- 0 # set NA to 0 ?
	rawCounts <- rawCounts[which(!is.na(rowSums(rawCounts))),]
    return (list(counts=as.matrix(rawCounts),
                 exp.design=exp.design))
}

#' Barplot of total count per sample
#' 
#' @param counts Matrix of read counts
#' @param exp.design Data frame of experimental design
#' @param filename Figure file name
#' @param out Logical variable for mode of plotting
barplotTC <- function(counts, exp.design, row.names.col, ylab="Total Read Count (in million)", filename="barplotTC.pdf", out=F) {
    df <- data.frame(sample=factor(colnames(counts), levels=colnames(counts)),
                     TC=colSums(counts)/10^6,
                     group=factor(exp.design[colnames(counts), 'group'],
                                  levels=unique(exp.design[colnames(counts), 'group'])),
                     color=factor(exp.design[colnames(counts), 'color'],
                                  levels=unique(exp.design[colnames(counts), 'color'])),
                     type=exp.design[colnames(counts), 'cell']
    )
    texts <- lapply(unique(df$type), function(x) {
        textGrob(x, gp=gpar(fontsize=14, fontface="bold", col=row.names.col[x]))
    })
    
    p <- ggplot(df,
                aes(x=sample, y=TC, fill=group)) +
        geom_bar(stat="identity") +
        xlab("Sample") +
        ylab(ylab) +
        scale_fill_manual(name="Condition", values=levels(df$color)) +
        theme_light() +
        theme(title           = element_text(size=18,face='bold'),
              legend.text     = element_text(size=12),
              legend.position = "top",
              axis.text       = element_text(size=10),
              axis.text.x     = element_text(angle=90, vjust=0.5,
                                             colour=row.names.col[df$type])
        )
    for (i in 1:length(texts)) {
        p <- p +
            annotation_custom(texts[[i]],
                              xmin=ifelse(i>1, sum(table(df$type)[unique(df$type)[1:(i-1)]]), 0)+table(df$type)[unique(df$type)[i]]/2+1/2,
                              xmax=ifelse(i>1, sum(table(df$type)[unique(df$type)[1:(i-1)]]), 0)+table(df$type)[unique(df$type)[i]]/2+1/2,
                              #xmin=(2*i-1)*nrow(df)/(2*length(texts)),
                              #xmax=(2*i-1)*nrow(df)/(2*length(texts)),
                              ymin=-max(df$TC)/30,
                              ymax=-max(df$TC)/30)
    }
    p <- p + coord_cartesian(ylim=c(0,max(df$TC)), clip="off")
    print(p)
    if (out) ggsave(p, filename=filename)
}

#' Keep genes with cpm >= cpm.threhold in at least num.threshold samples
#'
#' @param counts Matrix of read counts
#' @param exp.design Data frame of experimental design
#' @param cpm.threshold Lower threshold of cpm
#' @param num.threshold Lower threshold of number of samples. Default: NA, median count in at least one condition higher than cpm.threshold
#' @return Matrix of kept read counts
removeLowCpm <- function(counts, exp.design, cpm.threshold=1, num.threshold=NA) {
    if (!is.na(num.threshold)) {
        count.threshold <- rowSums(cpm(counts)>=cpm.threshold) >= num.threshold
    } else {
        count.threshold <- apply(cpm(counts), 1, function(x) max(tapply(x, factor(exp.design[colnames(counts), "group"]), median)) >= cpm.threshold) 
    }
    return (counts[count.threshold, , drop=F])
}

#' Boxplot of count distribution per sample
#'
#' @param counts Matrix of read counts
#' @param exp.design Data frame of experimental design
#' @param filename Figure file name
#' @param out Logical variable for mode of plotting
boxplotCounts <- function(counts, exp.design, row.names.col, ylab="Read count distribution", filename='boxplotCounts.pdf', out=F) {
    df <- melt(counts)
    colnames(df) <- c("gene", "sample", "count")
    df <- data.frame(df,
                     group=factor(exp.design[df$sample, 'group'],
                                  levels=unique(exp.design[df$sample, 'group'])),
                     color=factor(exp.design[df$sample, 'color'],
                                  levels=unique(exp.design[df$sample, 'color'])),
                     type=exp.design[df$sample, 'cell']
    )
    texts <- lapply(unique(exp.design$cell), function(x) {
        textGrob(x, gp=gpar(fontsize=14, fontface="bold", col=row.names.col[x]))
    })
    p <- ggplot(df,
                aes(x=sample, y=count, fill=group)) +
        geom_boxplot() +
        xlab("Sample") +
        ylab(ylab) +
        scale_fill_manual(name="Condition", values=levels(df$color)) +
        theme_light() +
        theme(title           = element_text(size=18,face='bold'),
              legend.text     = element_text(size=12),
              legend.position = "top",
              axis.text       = element_text(size=10),
              axis.text.x     = element_text(angle=90, vjust=0.5,
                                             colour=row.names.col[exp.design[levels(df$sample), "cell"]])
        )
    for (i in 1:length(texts)) {
        p <- p +
            annotation_custom(texts[[i]],
                              #xmin=ifelse(i>1, sum(table(df$type)[unique(df$type)[1:(i-1)]]), 0)+table(df$type)[unique(df$type)[i]]/2+1/2,
                              #xmax=ifelse(i>1, sum(table(df$type)[unique(df$type)[1:(i-1)]]), 0)+table(df$type)[unique(df$type)[i]]/2+1/2,
                              xmin=(2*i-1)*ncol(counts)/(2*length(texts)),
                              xmax=(2*i-1)*ncol(counts)/(2*length(texts)),
                              ymin=min(counts)-(max(counts)-min(counts))/30,
                              ymax=min(counts)-(max(counts)-min(counts))/30)
    }
    p <- p + coord_cartesian(ylim=c(min(counts), max(counts)), clip="off")
    print(p)
    if (out) ggsave(p, filename=filename)
}

#' Density plot of count distribution per sample
#'
#' @param counts Matrix of read counts
#' @param exp.design Data frame of experimental design
#' @param filename Figure file name
#' @param out Logical variable for mode of plotting
densityPlot <- function(counts, exp.design, count.type="Count", filename='densityCounts.pdf', legend.position="none", out=F) {
    df <- melt(counts)
    colnames(df) <- c("gene", "Sample", "count")
    p <- ggplot(df, aes(x=count, color=Sample)) +
        geom_density() +
        xlab(count.type) +
        ylab("Density") +
        theme(title           = element_text(size=18,face='bold'),
              legend.text     = element_text(size=10),
              legend.position = legend.position,
              axis.text       = element_text(size=12)
              )
    print(p)
    if (out) ggsave(p, filename=filename)
}

#' Density plot of count distribution per sample
#'
#' @param counts Matrix of read counts
#' @param exp.design Data frame of experimental design
#' @param filename Figure file name
#' @param out Logical variable for mode of plotting
PCAPlot <- function(counts, exp.design, title="All samples", filename='PCAPlot.pdf', out=F, colour.var='group', shape.var='batch', width=7, height=7) {
    constant.idx <- which(apply(counts, 1, function(x) length(unique(x))==1))
    if (length(constant.idx) > 0) counts[constant.idx,] <- counts[constant.idx,]*(1+(1e-3)*rnorm(ncol(counts), mean=0, sd=1))
    pca <- prcomp(t(counts), center=T, scale.=T)
    pca.dims <- list(c(1,2), c(3,2), c(1,3))
    pca.df <- data.frame(pca$x[, sort(unique(unlist(pca.dims)))],
                         exp.design[rownames(pca$x),]
                         )
    pca.plot <- ggplot(data=pca.df) +
        theme_grey() +
            theme(
                title            = element_text(size=18, face="bold"),
                axis.text        = element_text(size=12, face="plain"),
                axis.title       = element_text(size=14, face="plain"),
                legend.title     = element_text(size=13, face="bold"),
                legend.text      = element_text(size=12, face="plain"),
                legend.direction = "vertical",
                legend.position  = "none"
                )
    pca.plots <- lapply(1L:length(pca.dims), function(i) {
        p <- pca.plot +
            #ggtitle(ifelse(i==1, title, "")) +
            geom_point(aes_string(x=paste0('PC', pca.dims[[i]][1]),
                                  y=paste0('PC', pca.dims[[i]][2]),
                                  colour=colour.var, shape=shape.var), size=3) + #if (length(unique(pca.df$cell)) > 1) 'cell' else 'batch'
            geom_text(aes_string(x=paste0('PC', pca.dims[[i]][1]),
                                 y=paste0('PC', pca.dims[[i]][2]),
                                 label="individual"), size=2)+#, check_overlap=F) +
            xlab(paste0("PC",
                        pca.dims[[i]][1],
                        " (",
                        round((summary(pca)$importance)[2, pca.dims[[i]][1]], 2), ")")) +
            ylab(paste0("PC",
                        pca.dims[[i]][2],
                        " (",
                        round((summary(pca)$importance)[2, pca.dims[[i]][2]], 2), ")")) +
            scale_colour_manual(name="Condition", values=as.vector(unique(pca.df$color))) +
            scale_shape_manual(name=str_to_title(shape.var), #if (length(unique(pca.df$cell)) > 1) "cell" else "batch", 
                               values=if (length(unique(pca.df[,shape.var])) == 4) c(15, 16, 17, 6) else 
                                   if (length(unique(pca.df[,shape.var])) == 3) c(15, 16, 17) else 
                                       if (length(unique(pca.df[,shape.var])) == 2) c(15, 16) else
                                           if (length(unique(pca.df[,shape.var])) == 1) c(16) else 
                                               c(15:18, 0:14)[1:length(unique(pca.df[,shape.var]))],
            guide=guide_legend(title.position='top')) +
                theme(legend.position = "none")
        return (p)
    })
    legend <- gtable_filter(ggplot_gtable(ggplot_build(pca.plots[[1]] +
                                                       theme(legend.position="bottom"))),
                            "guide-box")
    ag <- ggarrange(plotlist=c(pca.plots, list(legend)), ncol=2, nrow=2)
    print(annotate_figure(ag, top = text_grob(title, color = "black", face = "bold", size = 16)))
    #if (out) ggsave(ag, filename=filename, width=width, height=height)
    if (out) ggsave(annotate_figure(ag, top = text_grob(title, color = "black", face = "bold", size = 16)), 
                    filename=filename, width=width, height=height)
}

#' Distance between points
#'
#' @param x A matrix with points in columns
#' @param method Method of distance: "pearson", "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Default: "pearson".
myDist <- function(x, method="pearson") {
    if (method %in% c("pearson")) {
        as.dist(1 - cor(x, use='p', method=method))
    } else {
        dist(t(x), method=method)
    }
}

#####################BACKUP



## scatter plot for two samples
## input: log normalized counts
## output: scatter plot
makeScatterPlot2 <- function(normdd, fc=1, RNAseq=F) {
    suppressPackageStartupMessages(require(MASS))
    m <- normdd[,1] - normdd[,2]
    xlim <- c(2,14)
    ylim <- c(2,14)
    main <- "Normalized expression levels lib1 vs. lib2"
    if (RNAseq) {
        xlim <- c(min(normdd), max(normdd))
        ylim <- c(min(normdd), max(normdd))
        main <- "RNA-seq Counts (transformed)"
    }
    ty.top1 <- 14/15*ylim[2]
    ty.top2 <- 13/15*ylim[2]
    ty.bottom <- ylim[1]+(ylim[2]-ylim[1])/15
    tx.left1 <- xlim[1]+(xlim[2]-xlim[1])/6
    tx.left2 <- xlim[1]+2*(xlim[2]-xlim[1])/3
    #print (colnames(normdd))
    #print (cor(normdd[,1], normdd[,2]))
    eqscplot(normdd[,1], normdd[,2], pch=16, cex=0.1, xlab=colnames(normdd)[1], ylab=colnames(normdd)[2], cex.axis=1.5, cex.lab=1.5, main=main, xlim=xlim, ylim=ylim)
    #qqplot(normdd[,1], normdd[,2], pch=16, cex=1, xlab=colnames(normdd)[1], ylab=colnames(normdd)[2], cex.axis=1.5, cex.lab=1.5, main=main, xlim=xlim, ylim=ylim)
    text(tx.left1,ty.top1, paste("Cor. coef: ", format(cor(normdd[,1], normdd[,2]), digits=4), sep=""), cex=1.5)
    text(tx.left1,ty.top2, paste("SD of diff: ", format(sd(m), digits=4), sep=""), cex=1.5)
    #text(tx.left2,ty.bottom, paste("Probe sets with ", format(2^fc, digits=2), "-fold change", sep=""), col=2, cex=1.5)
    #abline(fc,1, lty=2, col=2, lwd=2)
    #abline(-fc,1, lty=2, col=2, lwd=2)
    abline(0,1, col=4, lwd=2)
                                        #identify those with more than 2 fold change
                                        #fd=abs(m) > fc
                                        #points(normdd[fd,1], normdd[fd,2], pch=16, cex=0.1)
}

## scatter plot for all samples
## input: log normalized counts
## output: scatter plot
makeQCScatterPlots <- function(counts, OUT="", RNAseq=F, out) {
    file.Nb <- 1
    numrow <- 1
    numcol <- 2
    pdfunit <- 500
    if (out) {
        png(file=paste(OUT, "_", "pairs_scatter_plots_", file.Nb, ".png", sep=""), width=pdfunit*numcol, height=pdfunit*numrow)
    }
    count <- 0

    par(mfrow=c(numrow,numcol))
    for (i in 1:(ncol(counts)-1)) {
        for (j in (i+1):ncol(counts)){
            #if ((i==24 & j==26) || (i==25 & j==27)) {
                makeScatterPlot2(counts[,c(i,j)], RNAseq=RNAseq)
                count <- count+1
                if (count%%(numrow*numcol) == 0) {
                    file.Nb <- file.Nb + 1
                    if (out) {
                        invisible(dev.off())
                        png(file=paste(OUT, "_", "pairs_scatter_plots_", file.Nb, ".png", sep=""), width=pdfunit*numcol, height=pdfunit*numrow)
                    }
                    par(mfrow=c(numrow,numcol))
                }
            #}
        }
    }
    if (out) {
        invisible(dev.off())
    }
}

## sample clustering and PCA
## input: normalized counts
## output: clustering and PCA plot
makeQCClusteringAndPCA <- function(counts, conds, OUT, pcacolors, genes_nb = NULL, method = "single", is.trx = FALSE, legend.pos = "topright", out) {
    suppressPackageStartupMessages(require(maptools))
    suppressPackageStartupMessages(require(dendextend))
    
    numcon <- length(unique(colnames(counts)))
    numrep <- rle(colnames(counts))$lengths
    
    feat <- "genes"
    if (is.null(is.trx) | !is.trx) {
        feat <- "genes"
    } else {
        feat <- "transcripts"
    }

    ymfrow <- ifelse(is.null(genes_nb), 1, 2 - is.null(genes_nb[[2]]))
    if (ncol(counts) <= 27) {
        mfrow <- c(ymfrow, 2)
    } else {
        mfrow <- c(ymfrow, 1)
    }
    
    if (out) {
        pdf(file=paste(OUT,"_SampleCluster_PCA.pdf", sep=""),
            width=ifelse(ymfrow==2, 13, 13),
            height=ifelse(ymfrow==2, 15, 7.5), title="cluster")
    }

    par(mfrow=mfrow, pty="s", mar=c(7,5.5,5,2)+0.1)

    ## standard deviation for each gene
    sds <- apply(counts, 1, sd)
    
    if (length(genes_nb) == 0) {
        ## order of decreasing
        o <- order(sds, decreasing=T)
        #genes_nb <- c(nrow(counts), round(nrow(counts)/2), round(nrow(counts)/4), round(nrow(counts)/10))
        genes_nb <- c(nrow(counts), 1100, 1000, 100)#4000
        #genes_nb <- c(1000, 100)
        titles <- NULL
        for (i in 1:length(genes_nb)) {
            titles <- c(titles, sprintf(paste("%i", feat), as.integer(genes_nb[i])))
        }
        names(genes_nb) <- titles
        top.Genes <- lapply(1:length(genes_nb), function(k) {counts[o[1:genes_nb[k]], ]})
        names(top.Genes) <- c("Whole", rep("DE", length(top.Genes)-1))
    } else {
        top.Genes <- sapply(genes_nb, function(x) {
            if (length(x) > 0) {
                return (counts[intersect(x, rownames(counts)), , drop=F])
            } else {
                return (NULL)
            }
        }, USE.NAMES=T)
    }
    #print(genes_nb)
    #print(length(top.Genes))
    #print(dim(top.Genes[[1]]))
    ## #### A. BENHAMIAE PROJECT
    ## offset <- list(cbind(c(3,0,0,  0,0,0,  0,0,0,  0,0,0,  0,0,0),
    ##                      c(0,0,-3,  0,0,0,  0,0,0,  0,0,0,  0,0,0)),
    ##                cbind(c(0,0,-3,  0,2,0,  0,0,0,  0,0,0,  0,0,0),
    ##                      c(0.5,-2,0,  0,0,0,  0,0,0,  0,2,0,  0,0,0))
    ##                )
    ## groupoffset <- list(cbind(c(0,0,0,0,0), c(-14,12,0,14,-12)),
    ##                     cbind(c(0,0,0,0,0), c(-3,3,0,-3,-3)))
    ## textpos <- list(c(3,4,1, 2,2,3, 2,4,3, 4,3,2, 4,3,1),
    ##                 c(4,4,1, 2,3,2, 2,4,2, 1,3,2, 4,3,1)
    ##                 )
    ## ## textpos <- list(c(4,2,4,4,4,4), NULL
    ## ##                 )                   # in vitro
    ## ## textpos <- list(c(4,2,2,2), NULL
    ## ##                 )                   # galleria
    ## ## textpos <- list(c(4,4,2,2,4,4), NULL
    ## ##                 )                   # mouse
    ## ## textpos <- list(c(4,2,4,4,2,2,2,2), NULL
    ## ##                 )                   # galleria 12
    ## ## textpos <- list(rep(3,16), NULL
    ## ##                 )                   # all new
    labels <- list(c("A", "B"),
                   c("C", "D"),
                   c("E", "F"),
                   c("G", "H")
                   )
    pca.dims <- list(c(1, 2), c(1,3), c(2,3))
    cex.lab  <- ifelse (ncol(counts) <= 27, 1.5, 1)
    cex.axis <- ifelse (ncol(counts) <= 27, 1.5, 1)
    for (k in 1:length(top.Genes)) {
        if (! is.null(top.Genes[[k]])) {
            d <- as.dist(1 - cor(top.Genes[[k]]))
            #d <- dist(t(top.Genes[[k]]))
            hc <- hclust(d, method=method)
            dg <- as.dendrogram(hc)
            labels_colors(dg) <- pcacolors[hc$order]
            #desc.title <- sprintf(paste("Hierachical clustering with minimized variance, %i", feat), as.integer(nrow(top.Genes[[k]])))
            desc.title <- sprintf(paste("%i", feat), as.integer(nrow(top.Genes[[k]])))
            plot(dg, cex=1.6, cex.main=2, sub="", main="",
                 xlab="", ylab="Height", cex.lab=cex.lab, cex.axis=cex.axis)
            
            pca <- prcomp(t(top.Genes[[k]]), center=T, scale.=T)
            importance <- summary(pca)$importance
            print(importance)
            #desc.title <- sprintf(paste("Principal Components (1 & 2), %i", feat), as.integer(nrow(top.Genes[[k]])))
            desc.title <- sprintf(paste0(names(top.Genes)[k], ": %i ", feat), as.integer(nrow(top.Genes[[k]])))
            mtext(desc.title, side=3, line=3, outer=F, at=16.5, cex=2, font=2)
            
            lapply(pca.dims, function(pca.dim) {
                axlimits <- c(min(c(pca$x[, pca.dim[1]], pca$x[, pca.dim[2]])) * 1.2,
                              max(c(pca$x[, pca.dim[1]], pca$x[, pca.dim[2]])) * 1.2)

                propofvar1 <- importance[2, pca.dim[1]]
                propofvar2 <- importance[2, pca.dim[2]]
                xlab <- paste("PC ", pca.dim[1], " (fraction of total variance: ", signif(propofvar1, digits=2), ")", sep='')
                ylab <- paste("PC ", pca.dim[2], " (fraction of total variance: ", signif(propofvar2, digits=2), ")", sep='')

                                        #mtext(names(genes_nb)[k], side=3, line=3, outer=F, at=16, cex=cex.lab, font=2)
                                        #mtext(labels[[k]][1], side=3, line=0, outer=F, at=-1.5, cex=cex.lab, font=2)
                                        #mtext(labels[[k]][2], side=3, line=0, outer=F, at=16, cex=cex.lab, font=2)
                                        #print(rep(c(1:numcon), times=numrep))
                plot(pca$x[, pca.dim],
                     col=pcacolors,
                                        #col=rep(1:numcon, each=numrep),
                                        #col=rep(c("red", "blue", "darkgreen", "darkmagenta", "darkorange"), each=numrep),
                                        #col=rep(c("black", "red", "darkred", "blue", "darkblue", "green", "darkgreen", "magenta", "darkmagenta"), each=numrep),
                                        #col=rep(c("darkred", "darkblue", "darkgreen", "darkmagenta"), each=numrep),
                                        #col=rep(c(1:numcon), times=numrep),
                     pch=19, cex=1, cex.main=2, xlim=axlimits, ylim=axlimits,
                     xlab=xlab, ylab=ylab, cex.lab=cex.lab, cex.axis=cex.axis)
                point.labels <- setNames(rep("", ncol(counts)), colnames(counts))
                ## sel.labels <- c("1SC_C1.3", "1SC_H1.2", "1V_C1.3", "1V_H1.3",
                ##                 "20SC_C1.2", "20SC_H1.3", "20V_C1.2", "20V_H1.2",
                ##                 "8SC_C1.1", "8SC_H1.2", "8V_C1.1", "8V_H1.1")
                ##point.labels[sel.labels] <- gsub(sel.labels, pattern='[0-9]\\.[0-9]', replacement='')
                point.labels <- colnames(counts)
                pointLabel(pca$x[, pca.dim], labels=point.labels, allowSmallOverlap=F, cex=0.95, offset=0, col=pcacolors)
                                        #legend("bottomright", c("naive", "SC5314.8h", "101.8h", "SC5314.1d", "101.1d", "SC5314.3d", "101.3d", "SC5314.7d", "101.7d"),
                                        #       pch=19, col=c("black", "red", "darkred", "blue", "darkblue", "green", "darkgreen", "magenta", "darkmagenta"))
                ## legend("bottomright", c("glu.ev", "glu.wt", "xyl.ev", "xyl.wt"),
                ##        pch=19, col=c("darkred", "darkblue", "darkgreen", "darkmagenta"))
                #legend(legend.pos, conds, pch=19, col=unique(pcacolors))
            })
            
            
            ind <- min(k, 2)
            if (k > 2) {
                warning("offset or textpos would be adjusted after the second plot.")
            }
            ## text(pca$x[, pca.dim]+ offset[[ind]], labels=colnames(counts),
            ##      pos=textpos[[ind]],
            ##      adj=0, cex=2, offset=0.8)

            ## grouppos <- t(sapply(1:numcon, function(i) {
            ##     apply(pca$x[(numrep*(i-1)+1):(numrep*i), pca.dim], 2, mean)
            ## }))
            ## print(grouppos)
            
            ## #### A. BENHAMIAE PROJECT
            ## text(grouppos+groupoffset[[ind]], labels=c("Gp8", "Gp14", "K", "S", "Sa"),
            ##      adj=0, cex=2.5, offset=0.8)
            
        }
    }
    if (out) {
        invisible(dev.off())
    }
    pdf(file=paste(OUT,"_PCA_screeplots.pdf",sep=""), width=14, height=10, title="cluster")
    par(mfrow=c(2,2))
    for (k in length(top.Genes)) {
        if (! is.null(top.Genes[[k]])) {
            pca <- prcomp(t(top.Genes[[k]]), center=T, scale.=T)
            desc.title <- sprintf(paste("PCA screeplot, %i", feat), as.integer(nrow(top.Genes[[k]])))
            plot(pca, type="lines", main=desc.title)
        }
    }
    if (out) {
        invisible(dev.off())
    }
}

## histogram of p-values
## input: res, OUT_histoPName
## output: histogram
histoP <- function(res, OUT_histoPName, cond, out) {
    if (out) {
        pdf(file=paste(OUT_histoPName, ".pdf", sep=''))
    }
    
    hist(res, nclass=50, xlab="p-value", main=cond, col="skyblue")
    
    if (out) {
        invisible(dev.off())
    }
}

## MAplot of DE genes
## input: res, alpha, OUT_MAplotDEName 
## output: MAplot
MAplotDE <- function(A, M, p, main, fdr=0.05, OUT_MAplotDEName, out) {
    if (out) {
        pdf(file=paste(OUT_MAplotDEName, ".pdf", sep=''))
    }

    plot(A, M, main=main, xlab="Mean expression", ylab="log2FC", pch=16, col=ifelse(p<=fdr, "red", "black"), cex=ifelse(p<=fdr, 0.5, 0.5))
    abline(h=0, col="red")
    
    if (out) {
        invisible(dev.off())
    }
}

## Venn diagram plot
## input: results from decideTests, include=c("both", "down", "up"), output file name
## output: Venn diagram
vennplot <- function(res, include="both", colors, OUT_VennDiagramName, out) {
    if (out) {
        pdf(file=paste(OUT_VennDiagramName, ".pdf", sep=''))
    }
    
    vennDiagram(res, include, counts.col=colors, circle.col=c("red","blue","green","magenta"), cex=c(0.9, 0.9, 0.8))
    if (length(include)==1) {
       text(-2,-2, include, col=colors)
    }
    
    if (out) {
        invisible(dev.off())
    }
}

## heatmap plot
## input: counts (RPKM, normalized,...), output file name
## output: heatmap
heatplot <- function(counts, OUT_heatmapName, out) {
    suppressPackageStartupMessages(require(gplots))
    if (out) {
        pdf(file=paste(OUT_heatmapName, ".pdf", sep=''), width=8, height=6)
    }

    heatmap.2(counts, labRow="", labCol=NULL, scale="row", density.info="density", trace="none")

    if (out) {
        invisible(dev.off())
    }
}


## get arguments
## input: nomVariable, valeur
set <- function(nomVariable, valeur) {
    if (is.character(valeur)) {
        commande <- paste(nomVariable, "<<-\"", valeur, "\"", sep="")
    }
    else {
        commande <- paste(nomVariable, "<<-", valeur, sep="")
    }
    expr <- try(parse(text=commande), TRUE)
    eval(expr)
}

## initialize arguments
initArgs <- function() {
    nombreArg   <<- 0
    nomArgs     <<- list()
    valDefaut   <<- list()
    attributs   <<- list()
    description <<- list()
}

## set argument values
addArgs <- function(nom, def, args, desc) {
    nombreArg                <<- nombreArg + 1
    nomArgs[[nombreArg]]     <<- nom
    valDefaut[[nombreArg]]   <<- def
    attributs[[nombreArg]]   <<- args
    description[[nombreArg]] <<- desc
}

## help
printHelp <- function() {
    cat("RNAseq data analysis with limma
Command: limma2.R [options] designFile
Options:
")
    for(i in 1:length(nomArgs)) {
        cat(paste("\t",format(paste(attributs[[i]],':',sep=''), width=30),description[[i]]," (Default: ",sep=""))
        if(!is.na(valDefaut[[i]])) {
            cat(paste(valDefaut[[i]],")\n",sep=''))
        }
        else {
            cat(paste("<OUT>_",outputFilesSuffix[i],")\n",sep=''))
        }
    }
}
