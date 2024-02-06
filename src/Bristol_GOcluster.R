##########################################################################
## parse gmt RUN ONCE to create mapping GMT ID - GO ID ----
FIRSTRUN <- FALSE
if (FIRSTRUN) {
    options(stringsAsFactors = F)
    rm(list=ls())
    DIR <- '/scratch/local/permanent/ttran/BEAt_DKD/'
    setwd(paste0(DIR, '/data/'))
    
    ## gmt downloaded from MSigDB
    library(xml2)
    library(rvest) 
    library(stringr)
    #library(ramiGO)
    
    # invisible(sapply(list.files(path="download", pattern="c5.", full.names=T), function(x) {
    #     tab <- read.delim2(x, header = F, sep='\r', stringsAsFactors = F)[,1]
    #     tab.list <- strsplit(tab, split='\t')
    #     gmt.map <- do.call(rbind, lapply(tab.list, function(y) { # mclapply does not work!
    #         webpage <- read_html(y[2])
    #         tds <- html_text(html_nodes(webpage, "td"))
    #         goid <- grep("^GO:", tds, value=T)
    #         return (c(y[1], goid))
    #     }))
    #     write.table(gmt.map, file=paste0(gsub(x, pattern='download/', replacement='clean/gmt/'), "_map.tsv"), 
    #                 quote=F, sep='\t', col.names=F, row.names=F)
    # }))
    
    ## miRNA
    library(GO.db)
    goterms <- do.call(rbind, mclapply(as.list(GOTERM), mc.cores=detectCores(), function(x) {
        c(Term=Term(x), GOID=GOID(x))
    }))
    rownames(goterms) <- goterms[, "Term"]
    dir.create("clean/gmt/miRNA", showWarnings=T)
    invisible(sapply(list.files(path="download/miRNA", pattern="GO_[_A-Z]+validated_miRTarBase_all", full.names=T), function(x) {
        tab <- read.delim2(x, header = F, sep='\r', stringsAsFactors = F)[,1]
        tab.list <- strsplit(tab, split='\t')
        tab.term <- sapply(tab.list, function(y) y[1])
        tab.x <- goterms[intersect(tab.term, rownames(goterms)), ] ## due to different annotation versions
        tab.x[,"Term"] <- gsub(tab.x[,"Term"], pattern="\\+", replacement="PLUS")
        tab.x[,"Term"] <- gsub(tab.x[,"Term"], pattern="'",   replacement="")
        tab.x[,"Term"] <- gsub(tab.x[,"Term"], pattern="-| |\\.|\\(|\\)|/|,|:|>|=|#|\\[|\\]", replacement="_")
        tab.x[,"Term"] <- gsub(tab.x[,"Term"], pattern="__", replacement="_")
        tab.x[,"Term"]  <- paste0("GO_", toupper(tab.x[,"Term"]))
        write.table(tab.x[, c("Term", "GOID")], file=paste0(gsub(x, pattern='download/', replacement='clean/gmt/'), "_map.tsv"),
                    quote=F, sep='\t', col.names=F, row.names=F)
    }))
}
##-----

LIBRARIES <- c("GOSemSim", "GO.db",  "org.Hs.eg.db",
               "WGCNA", "ggplot2", "circlize", "dendextend", "gplots", "ComplexHeatmap",
               "assertthat", "parallel")
invisible(sapply(LIBRARIES, function(lib) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
}))

options(stringsAsFactors = F)
rm(list=ls())
DIR <- '/scratch/local/permanent/ttran/BEAt_DKD/'
#DIR <- '~/Work/BEAt-DKD/'
PROJECT <- "Endothelial"
font.size <- 12
GOSIZELIM <- 500
PLIM <- 0.05
QLIM <- 0.1
WEBSAFECOLORS <- labels2colors(1:340)
WEBSAFECOLORS <- WEBSAFECOLORS[!grepl(WEBSAFECOLORS, pattern='[0-9]|white|yellow$|^pale')]
WEBSAFECOLORS <- setdiff(WEBSAFECOLORS, c("lightcyan", "ivory", "slateblue", "mistyrose", "honeydew",
                                          "cornsilk", "aliceblue", "lavenderblush", "snow", "oldlace", 
                                          "linen", "blanchedalmond", "lavender", "moccasin", "mintcream"))
COLORGROUPS <- list()
COLORGROUPS[["BEAtDKDcell"]] <- setNames(
    c("coral", "gray40", "turquoise3", "gold2", "gold2", "turquoise3", "brown"),
    c("ABN",   "GEnC",   "K29",        "PT34",  "HMEC",  "ECFC",       "RAEC"))
noclusterCell <- F
setwd(paste0(DIR, '/data/clean/gmt'))
gmtfiles <- grep("all", list.files(pattern="^c5"), invert=T, value=T)
gmt.id <- setNames(mclapply(gmtfiles, mc.cores=3, function(x) {
    read.table(x, row.names=1)
}), toupper(gsub(gmtfiles, pattern="c5.|.v7.1.symbols.gmt_map.tsv", replacement = "")))

godatas <- sapply(c('BP', 'MF', 'CC'), function(x) {
    godata(OrgDb="org.Hs.eg.db", ont=x, computeIC = F)
})
titles <- c(BP="GO Biological Process", CC="GO Cellular Component", MF="GO Molecular Function")
if (PROJECT %in% c("Bristol", "Helsinki", "Helsinki_Basal")) {
    cells <- c("ABN", "GEnC", "K29", "PT34")
} else if (PROJECT %in% "Lund") {
    cells <- c("ECFC")
} else if (PROJECT %in% "Endothelial") {
    cells <- c("GEnC", "ECFC", "HMEC", "RAEC")
} 

if (PROJECT %in% c("Bristol", "Helsinki")) {
    ############################################
    ## GO clustering
    setwd(paste0(DIR, '/limma/output/', PROJECT, '/pathway_combined'))
    #setwd('~/Work/BEAt-DKD/Pathway_enrichment/')
    enrfiles <- grep("^GO_.+withOPLS_stat.tsv", list.files(), value=T)
    
    #' param filtering List of index vectors of the same length
    #' param filtering.lim Vector of limits, with the same length as filtering vectors
    filterPathway <- function(tab, filtering = NULL, filtering.lim = NULL, groups=NULL) {
        if (is.null(filtering)) return (tab)
        colnames.1 <- gsub(colnames(tab)[filtering[[1]]], pattern="_pvalue$|_qvalues$|_NES$", replacement="")
        ## check if filtering is correctly formatted
        if (length(filtering) > 1) {
            correct.order <- all(sapply(2:length(filtering), function(i) {
                colnames.i <- gsub(colnames(tab)[filtering[[i]]], pattern="_pvalue$|_qvalues$|_NES$", replacement="")
                return (all(colnames.i==colnames.1))
            }))
            stopifnot(correct.order)
        }
        ## filtering on each row of tab
        ind.f <- do.call(cbind, mclapply(1:length(filtering[[1]]), mc.cores=length(filtering[[1]]), function(j) {
            apply(do.call(cbind, mclapply(1:length(filtering), mc.cores=length(filtering), function(i) {
                ifelse(is.na(tab[, filtering[[i]][j]]), NA, as.numeric(tab[, filtering[[i]][j]]) <= filtering.lim[i])
            })), 1, all)
        }))
        colnames(ind.f) <- colnames.1
        if (!is.null(groups)) {
            ind.group.f <- apply(do.call(cbind, mclapply(groups, mc.cores=length(groups), function(gr) {
                #apply(ind.f[, grep(gr, colnames(ind.f))], 1, all)
                apply(ind.f[, grep(gr, colnames(ind.f))], 1, function(x) {
                    #length(which(x==T)) > 1
                    ## more than 1 filtering satisfaction from either multiblock or individual
                    sum(na.omit(x[!grepl("loadings", grep(gr, colnames(ind.f), value=T))])) > 1 || sum(na.omit(x[grepl("loadings", grep(gr, colnames(ind.f), value=T))])) > 1 
                })
            })), 1, any) # in at least 1 cell type
        } else {
            ind.group.f <- apply(ind.f, 1, any)
        }
        ind.group.f[is.na(ind.group.f)] <- FALSE
        return (tab[ind.group.f, ])
    }
    
    res <- mclapply(enrfiles, mc.cores=length(enrfiles), function(x) {
        ontology <- strsplit(x, split='_')[[1]][2]
        tab.init <- read.table(x, header=T, sep='\t', quote='')
        if (PROJECT=="Bristol") {
            tab.init <- tab.init[, c(1:3,
                                     grep("RNAseq", colnames(tab.init)),
                                     grep("Prot", colnames(tab.init)),
                                     grep("loadings_genes", colnames(tab.init)),
                                     grep("loadings_proteins", colnames(tab.init)))]
        }
        if (PROJECT=="Helsinki") {
            tab.init <- tab.init[, c(1:3, 
                                     grep("_RNAseq", colnames(tab.init)),
                                     grep("_miRNAseq", colnames(tab.init)), 
                                     grep("_Prot", colnames(tab.init)))]
        }
        tab.init <- tab.init[, c(1:3, 
                                 grep("ABN", colnames(tab.init)), 
                                 grep("GEnC", colnames(tab.init)), 
                                 grep("K29", colnames(tab.init)), 
                                 grep("PT34", colnames(tab.init)))]
        tab.init <- tab.init[, c(1:3, 
                                 grep("NES$", colnames(tab.init)), 
                                 grep("pvalue", colnames(tab.init)), 
                                 grep("qvalue", colnames(tab.init)), 
                                 grep("core_enrichment", colnames(tab.init)))]
        
        tab <- filterPathway(tab.init, filtering = list(p=grep("pvalue$", colnames(tab.init)),
                                                        q=grep("qvalues$", colnames(tab.init))),
                             filtering.lim = c(p=PLIM, q=QLIM),
                             groups = cells)
        rownames(tab) <- gsub(tab$ID, pattern="^GO", replacement=paste0("GO_", ontology))
        tab$ID <- gmt.id[[ontology]][tab$ID,1]
        
        gos_ss <- sapply(tab$ID, function(y) {
            sapply(tab$ID, function(z) {
                if (exists(".GOSemSimEnv")) rm(".GOSemSimEnv")
                goSim(y, z, semData=godatas[[ontology]], measure='Wang')
            })
        })
        colnames(gos_ss) <- rownames(tab)
        rownames(gos_ss) <- rownames(tab)
        gos_dist <- 1 - gos_ss
        gos_hc <- hclust(as.dist(gos_dist), method="ward.D")
        hcut <- 1.4
        clus <- cutree(gos_hc, h=hcut)
        clunum <- length(unique(clus))
        colors <- WEBSAFECOLORS[1:clunum]
        dend <- as.dendrogram(gos_hc) %>% set("branches_k_color", value=colors, k=clunum) #
        
        ##- MDS ----
        mds <- data.frame(cmdscale(gos_dist))
        colnames(mds) <- c("x", "y")
        mds$go <- rownames(mds)
        mds$col <- factor(get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))])
        p <- ggplot(mds, aes(x=x, y=y)) + 
            geom_point(aes(col=col), show.legend = F) + 
            geom_text(aes(label=go, colour=col), size=1.5, check_overlap = T, show.legend = F,
                      vjust=-1) +
            scale_colour_manual(values=levels(mds$col)) +
            xlim(-0.7, 0.7) +
            xlab("Coordinate 1") + ylab("Coordinate 2") +
            theme_bw()
        ggsave(p, filename=gsub(x, pattern=".tsv", replacement = '_mds.pdf'), width=8, height=8)
        ##-----
        
        ##- clustering ----
        pdf(gsub(x, pattern=".tsv", replacement = ifelse(noclusterCell, '.pdf', '_clusteredCell.pdf')), 
            width=22, height=switch(ontology,
                                    BP=nrow(tab)/4,
                                    CC=nrow(tab)/3,#4
                                    MF=nrow(tab)/2.8))
        # ha <- rowAnnotation(foo=anno_barplot(as.matrix(tab[, c("Cluster.frequency", "Nocluster.frequency")]),
        #                                      axis_param=list(gp=gpar(fontsize=font.size-4), labels_rot=0),
        #                                      gp=gpar(fill=c(ifelse(UPREG, "red", "green"), "grey50"))),
        #                     annotation_label="Number of genes (log2)", 
        #                     annotation_legend_param=list(
        #                         legend_position="left",
        #                         legend_direction="horizontal"
        #                     ),
        #                     annotation_name_gp=gpar(fontsize=font.size-2),
        #                     show_legend=T,
        #                     #gp=gpar(fontsize=font.size),
        #                     #annotation_name_side="bottom",
        #                     width=unit(6, "cm"))
        #hb <- HeatmapAnnotation(foo=anno_block(gp = gpar(fill = 1:4)))
        hm.tab <- as.matrix(tab[, grep("_NES$", colnames(tab)), drop=F])
        hm.tab[is.na(hm.tab)] <- 0
        colnames(hm.tab) <- sub(sub(sub(colnames(hm.tab), pattern="RNAseq|loadings_genes", replacement = "transcriptomics"),
                                    pattern="Prot|loadings_proteins", replacement = "proteomics"),
                                pattern="_NES", replacement = "")
        col.split <- sapply(strsplit(colnames(hm.tab), split="_"), function(hmcol) hmcol[1])
        
        colnames(hm.tab) <- sub(sub(sub(sub(colnames(hm.tab),
                                            pattern="ABN", replacement="Pod"),
                                        pattern="GEnC", replacement="GEC"),
                                    pattern="K29", replacement="MC"),
                                pattern="PT34", replacement="PTC")
        hm <- Heatmap(hm.tab,
                      column_title=paste0("Selected ", titles[ontology], ": enriched (p < ", PLIM, 
                                          ", q < ", QLIM, ") in at least one cell type\n * p < 0.05, ** p < 0.01"),
                      #height = unit(nrow(tab)*0.03*ifelse(grepl("mbf1", x), 2, 1), "npc"),
                      column_title_gp = gpar(fontsize=15, fontface='bold'),
                      column_split = if (noclusterCell) NULL else col.split,#rep(cells, each=4),
                      heatmap_legend_param=list(title="Normalized Enrichement Score",
                                                title_position="leftcenter-rot",
                                                legend_direction="vertical",
                                                legend_height=unit(2,'inches'),
                                                title_gp=gpar(fontsize=font.size),
                                                labels_gp=gpar(fontsize=font.size)),
                      cell_fun=function(j, i, x, y, width, height, col) {
                          #grid.text(sprintf("%.2f", as.matrix(tab[, grep("_NES$", colnames(tab)), drop=F])[i, j]), x, y, gp=gpar(fontsize=10))
                          grid.text(if (is.na(as.matrix(tab[, grep("_pvalue$", colnames(tab)), drop=F])[i, j])) ""
                                    else if (as.matrix(tab[, grep("_pvalue$", colnames(tab)), drop=F])[i, j] < 1e-3) "***"
                                    else if (as.matrix(tab[, grep("_pvalue$", colnames(tab)), drop=F])[i, j] < 1e-2) "**"
                                    else if (as.matrix(tab[, grep("_pvalue$", colnames(tab)), drop=F])[i, j] < 0.05) "*"
                                    else "", 
                                    x, y, gp=gpar(fontsize=font.size))
                      },
                      col=blueWhiteRed(50),
                      # col=colorRamp2(
                      #     c(0, 5, max(tab[, "Corrected.P.value", drop=F])),
                      #     c("white", "blue", "darkblue"), space="RGB"),
                      cluster_rows = dend,
                      cluster_columns = noclusterCell,
                      row_dend_width = unit(4, "cm"),
                      width = unit(20, "cm"),
                      show_column_names = T,
                      row_names_gp = gpar(col=get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))],
                                          fontsize=font.size),
                      row_names_max_width = max_text_width(tab$ID, #paste(rep(" ", max.length), collapse=''),
                                                           gp=gpar(fontsize=font.size)),
                      column_names_gp = gpar(col=COLORGROUPS[["BEAtDKDcell"]][col.split]) #rep(cells, each=4)
                      #right_annotation = ha
        )
        draw(hm, 
             heatmap_legend_side="left",
             #annotation_legend_list=list(Legend(labels="A", labels_gp=gpar(fontsize=2*font.size, fontface="bold"))),
             #annotation_legend_side="left",
             #labels=c(paste0(ifelse(UPREG, "Up", "Down"), "-regulated genes"), "Total number of genes"),
             #labels=c(paste0("Fold-change ratio", ifelse(UPREG, " < 0.5", " > 2")), "Total number of genes"),
             #legend_gp=gpar(fill=c(ifelse(UPREG, "red", "green"), "grey50")))), 
             #legend_title_gp=gpar(fontsize=font.size+10),
             #labels_gp=gpar(fontsize=font.size-2),
             ## margins
             padding = unit(c(2, -8, 1, 0), "inches"),
             newpage=F)
        dev.off()
        ##-----
        
        ##- genes for each GO cluster ----
        genes.clus <- setNames(lapply(1:clunum, function(clu) {
            clu.goid <- names(clus[clus==clu])
            if (any(grepl("IMMUNE", clu.goid))) { #plum
                clu.goid <- c(clu.goid, "GO_BP_LEUKOCYTE_MIGRATION", names(clus[clus %in% unique(clus[grep("VIRUS", names(clus))])]))
                clu.goid <- setdiff(clu.goid, "GO_BP_VIRAL_GENE_EXPRESSION")
            } else if (any(grepl("MITOCHONDRIAL_TRANSLATION", clu.goid))) { #red
                clu.goid <- c(clu.goid, names(clus[clus %in% unique(clus[grep("ELECTRON_TRANSPORT|RESPIRATORY_CHAIN|MITOCHONDRION_ORGANIZATION", names(clus))])]))
            } else if (any(grepl("RNA_PROCESSING", clu.goid))) { #turquoise
                clu.goid <- c(clu.goid, names(clus[clus %in% unique(clus[grep("RIBOSOME", names(clus))])]))
                clu.goid <- setdiff(clu.goid, "GO_BP_MEMBRANE_BIOGENESIS")
            } else if (any(grepl("HEART_DEVELOPMENT", clu.goid))) { #lightgreen
                clu.goid <- c(clu.goid, names(clus[clus %in% unique(clus[grep("MORPHOGENESIS", names(clus))])]))
            }
            
            tab.clu <- tab[clu.goid, grep("core_enrichment", colnames(tab), value=T)]
            Reduce(union, apply(tab.clu, 1, function(gclu) {
                Reduce(union, strsplit(na.omit(gclu), split=','))
            }))
        }), colors)
        ##-----
        
        return (list(genes.clus=genes.clus, mds=mds))
    })
    
    names(res) <- gsub(enrfiles, pattern='_HelsinkiIR_withOPLS_stat.tsv', replacement = '')
    save(res, file='goclustering.RData')
} else if (PROJECT %in% c("Helsinki_Basal")) {
    ############################################
    ## EV cell-type enriched
    setwd(paste0(DIR, '/limma/output/', PROJECT, '/pathway_combined'))
    enrfiles <- grep("^GO_.+_stat.tsv", list.files(), value=T)
    #' param filtering List of index vectors of the same length
    #' param filtering.lim Vector of limits, with the same length as filtering vectors
    filterPathway <- function(tab, filtering = NULL, filtering.lim = NULL, groups=NULL) {
        if (is.null(filtering)) return (tab)
        colnames.1 <- gsub(colnames(tab)[filtering[[1]]], pattern=".pvalue$|.qvalue$|.GeneRatio$", replacement="")
        ## check if filtering is correctly formatted
        if (length(filtering) > 1) {
            correct.order <- all(sapply(2:length(filtering), function(i) {
                colnames.i <- gsub(colnames(tab)[filtering[[i]]], pattern=".pvalue$|.qvalue$|.GeneRatio$", replacement="")
                return (all(colnames.i==colnames.1))
            }))
            stopifnot(correct.order)
        }
        ## filtering on each row of tab
        ind.f <- do.call(cbind, mclapply(1:length(filtering[[1]]), mc.cores=length(filtering[[1]]), function(j) {
            apply(do.call(cbind, mclapply(1:length(filtering), mc.cores=length(filtering), function(i) {
                ifelse(is.na(tab[, filtering[[i]][j]]), NA, as.numeric(tab[, filtering[[i]][j]]) <= filtering.lim[i])
            })), 1, all)
        }))
        colnames(ind.f) <- colnames.1
        if (!is.null(groups)) {
            ind.group.f <- apply(do.call(cbind, mclapply(groups, mc.cores=length(groups), function(gr) {
                apply(ind.f[, grep(gr, colnames(ind.f))], 1, function(x) {
                    ## more than 1 filtering satisfaction from either multiblock or individual
                    sum(na.omit(x[!grepl("loadings", grep(gr, colnames(ind.f), value=T))])) > 1 || sum(na.omit(x[grepl("loadings", grep(gr, colnames(ind.f), value=T))])) > 1 
                })
            })), 1, any) # in at least 1 cell type
        } else {
            ind.group.f <- apply(ind.f, 1, any)
        }
        ind.group.f[is.na(ind.group.f)] <- FALSE
        return (tab[ind.group.f, ])
    }
    
    res <- mclapply(enrfiles, mc.cores=length(enrfiles), function(x) {
        ontology <- strsplit(x, split='_')[[1]][2]
        tab.init <- read.table(x, header=T, sep='\t', quote='')

        tab.init <- tab.init[, c(1:3, 
                                 grep("ABN", colnames(tab.init)), 
                                 grep("GEnC", colnames(tab.init)), 
                                 grep("K29", colnames(tab.init)), 
                                 grep("PT34", colnames(tab.init)))]
        tab.init <- tab.init[, c(1:3, 
                                 grep("_RNAseq", colnames(tab.init)),
                                 grep("_miRNAseq", colnames(tab.init)), 
                                 grep("_Prot", colnames(tab.init)))]
        tab.init <- tab.init[, c(1:3, 
                                 grep("GeneRatio", colnames(tab.init)),
                                 grep("BgRatio", colnames(tab.init)),
                                 grep("pvalue", colnames(tab.init)),
                                 grep("qvalue", colnames(tab.init)),
                                 grep("geneID", colnames(tab.init)))]

        tab <- filterPathway(tab.init, filtering = list(p=grep("pvalue$", colnames(tab.init)),
                                                        q=grep("qvalue$", colnames(tab.init))),
                             filtering.lim = c(p=PLIM, q=QLIM),
                             groups = cells)
        rownames(tab) <- gsub(tab$ID, pattern="^GO", replacement=paste0("GO_", ontology))
        tab$ID <- gmt.id[[ontology]][tab$ID,1]
        if (nrow(tab) > 1) {
            gos_ss <- sapply(tab$ID, function(y) {
                sapply(tab$ID, function(z) {
                    if (exists(".GOSemSimEnv")) rm(".GOSemSimEnv")
                    goSim(y, z, semData=godatas[[ontology]], measure='Wang')
                })
            })
            colnames(gos_ss) <- rownames(tab)
            rownames(gos_ss) <- rownames(tab)
            gos_dist <- 1 - gos_ss
            gos_hc <- hclust(as.dist(gos_dist), method="ward.D")
            hcut <- 1.4
            clus <- cutree(gos_hc, h=hcut)
            clunum <- length(unique(clus))
            colors <- WEBSAFECOLORS[1:clunum]
            dend <- as.dendrogram(gos_hc) %>% set("branches_k_color", value=colors, k=clunum) #
            
            ##- MDS ----
            mds <- data.frame(cmdscale(gos_dist))
            colnames(mds) <- c("x", "y")
            mds$go <- rownames(mds)
            mds$col <- factor(get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))])
            p <- ggplot(mds, aes(x=x, y=y)) + 
                geom_point(aes(col=col), show.legend = F) + 
                geom_text(aes(label=go, colour=col), size=1.5, check_overlap = T, show.legend = F,
                          vjust=-1) +
                scale_colour_manual(values=levels(mds$col)) +
                xlim(-0.7, 0.7) +
                xlab("Coordinate 1") + ylab("Coordinate 2") +
                theme_bw()
            ggsave(p, filename=gsub(x, pattern=".tsv", replacement = '_mds.pdf'), width=8, height=8)
            ##-----
            
            ##- clustering ----
            pdf(gsub(x, pattern=".tsv", replacement = ifelse(noclusterCell, '.pdf', '_clusteredCell.pdf')), 
                width=23, height=switch(ontology,
                                        BP=nrow(tab)/4,
                                        CC=nrow(tab)/3,
                                        MF=nrow(tab)/0.8))

            hm.tab <- as.matrix(tab[, grep(".GeneRatio$", colnames(tab)), drop=F])
            hm.tab <- apply(hm.tab, c(1,2), function(hmt) eval(parse(text=hmt)))
            hm.tab[is.na(hm.tab)] <- 0
            col.omics <- sapply(strsplit(colnames(hm.tab), split="_"), function(hmcol) hmcol[2])
            col.cells <- sapply(strsplit(colnames(hm.tab), split="_"), function(hmcol) hmcol[1])
            hm <- Heatmap(hm.tab,
                          column_title=paste0("Selected ", titles[ontology], "\n enriched (p < ", PLIM, 
                                              ", q < ", QLIM, ") in at least two omics on one cell type\n * p < 0.05, ** p < 0.01, *** p < 0.001"),
                          #height = unit(nrow(tab)*0.03*ifelse(grepl("mbf1", x), 2, 1), "npc"),
                          column_title_gp = gpar(fontsize=15, fontface='bold'),
                          column_split = if (noclusterCell) NULL else col.cells,#col.omics,#
                          heatmap_legend_param=list(title="Gene Ratio",
                                                    title_position="leftcenter-rot",
                                                    legend_direction="vertical",
                                                    legend_height=unit(2,'inches'),
                                                    title_gp=gpar(fontsize=font.size),
                                                    labels_gp=gpar(fontsize=font.size)),
                          cell_fun=function(j, i, x, y, width, height, col) {
                              #grid.text(sprintf("%.2f", as.matrix(tab[, grep("_NES$", colnames(tab)), drop=F])[i, j]), x, y, gp=gpar(fontsize=10))
                              grid.text(if (is.na(as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j])) ""
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 1e-3) "***"
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 1e-2) "**"
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 0.05) "*"
                                        else "", 
                                        x, y, gp=gpar(fontsize=font.size))
                          },
                          #col=blueWhiteRed(50),
                          col = colorpanel(50, low="white", high="red"),
                          cluster_rows = dend,
                          cluster_columns = noclusterCell,
                          row_dend_width = unit(4, "cm"),
                          width = unit(20, "cm"),
                          show_column_names = T,
                          row_names_gp = gpar(col=get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))],
                                              fontsize=font.size),
                          row_names_max_width = max_text_width(tab$ID, #paste(rep(" ", max.length), collapse=''),
                                                               gp=gpar(fontsize=font.size)),
                          column_names_gp = gpar(col=COLORGROUPS[["BEAtDKDcell"]][col.cells])
                          #right_annotation = ha
            )
            draw(hm, 
                 heatmap_legend_side="left",
                 #annotation_legend_list=list(Legend(labels="A", labels_gp=gpar(fontsize=2*font.size, fontface="bold"))),
                 #annotation_legend_side="left",
                 #labels=c(paste0(ifelse(UPREG, "Up", "Down"), "-regulated genes"), "Total number of genes"),
                 #labels=c(paste0("Fold-change ratio", ifelse(UPREG, " < 0.5", " > 2")), "Total number of genes"),
                 #legend_gp=gpar(fill=c(ifelse(UPREG, "red", "green"), "grey50")))), 
                 #legend_title_gp=gpar(fontsize=font.size+10),
                 #labels_gp=gpar(fontsize=font.size-2),
                 ## margins
                 padding = unit(c(2, -8, 1, 0), "inches"),
                 newpage=F)
            dev.off()
            ##-----
        }
        return (NULL)
        # ##- genes for each GO cluster ----
        # genes.clus <- setNames(lapply(1:clunum, function(clu) {
        #     clu.goid <- names(clus[clus==clu])
        #     if (any(grepl("IMMUNE", clu.goid))) { #plum
        #         clu.goid <- c(clu.goid, "GO_BP_LEUKOCYTE_MIGRATION", names(clus[clus %in% unique(clus[grep("VIRUS", names(clus))])]))
        #         clu.goid <- setdiff(clu.goid, "GO_BP_VIRAL_GENE_EXPRESSION")
        #     } else if (any(grepl("MITOCHONDRIAL_TRANSLATION", clu.goid))) { #red
        #         clu.goid <- c(clu.goid, names(clus[clus %in% unique(clus[grep("ELECTRON_TRANSPORT|RESPIRATORY_CHAIN|MITOCHONDRION_ORGANIZATION", names(clus))])]))
        #     } else if (any(grepl("RNA_PROCESSING", clu.goid))) { #turquoise
        #         clu.goid <- c(clu.goid, names(clus[clus %in% unique(clus[grep("RIBOSOME", names(clus))])]))
        #         clu.goid <- setdiff(clu.goid, "GO_BP_MEMBRANE_BIOGENESIS")
        #     } else if (any(grepl("HEART_DEVELOPMENT", clu.goid))) { #lightgreen
        #         clu.goid <- c(clu.goid, names(clus[clus %in% unique(clus[grep("MORPHOGENESIS", names(clus))])]))
        #     }
        #     
        #     tab.clu <- tab[clu.goid, grep("core_enrichment", colnames(tab), value=T)]
        #     Reduce(union, apply(tab.clu, 1, function(gclu) {
        #         Reduce(union, strsplit(na.omit(gclu), split=','))
        #     }))
        # }), colors)
        # ##-----
        # 
        # return (list(genes.clus=genes.clus, mds=mds))
    })
    
} else if (PROJECT %in% c("Lund")) {
    ############################################
    ## GO clustering
    setwd(paste0(DIR, '/limma/output/', PROJECT, '/pathway_combined'))
    #setwd('~/Work/BEAt-DKD/Pathway_enrichment/')
    enrfiles <- grep("^GO_.+withOPLS_stat.tsv", list.files(), value=T)
    
    #' param filtering List of index vectors of the same length
    #' param filtering.lim Vector of limits, with the same length as filtering vectors
    filterPathway <- function(tab, filtering = NULL, filtering.lim = NULL, groups = NULL) {
        if (is.null(filtering)) return (tab)
        colnames.1 <- gsub(colnames(tab)[filtering[[1]]], pattern=".pvalue$|.qvalues$|.NES$", replacement="")
        ## check if filtering is correctly formatted
        if (length(filtering) > 1) {
            correct.order <- all(sapply(2:length(filtering), function(i) {
                colnames.i <- gsub(colnames(tab)[filtering[[i]]], pattern=".pvalue$|.qvalues$|.NES$", replacement="")
                return (all(colnames.i==colnames.1))
            }))
            stopifnot(correct.order)
        }
        ## filtering on each row of tab
        ind.f <- do.call(cbind, mclapply(1:length(filtering[[1]]), mc.cores=length(filtering[[1]]), function(j) {
            apply(do.call(cbind, mclapply(1:length(filtering), mc.cores=length(filtering), function(i) {
                ifelse(is.na(tab[, filtering[[i]][j]]), NA, as.numeric(tab[, filtering[[i]][j]]) <= filtering.lim[i])
            })), 1, all)
        }))
        colnames(ind.f) <- colnames.1
        if (!is.null(groups)) {
            ind.group.f <- apply(do.call(cbind, mclapply(groups, mc.cores=length(groups), function(gr) {
                apply(ind.f[, grep(gr, colnames(ind.f))], 1, function(x) {
                    ## more than 1 filtering satisfaction from either multiblock or individual
                    #sum(na.omit(x[!grepl("loadings", grep(gr, colnames(ind.f), value=T))])) > 1 || sum(na.omit(x[grepl("loadings", grep(gr, colnames(ind.f), value=T))])) > 1 
                    ## more than 3 filtering satisfactions
                    sum(na.omit(x[grep(gr, colnames(ind.f), value=T)])) > 3
                })
            })), 1, any) # in at least 1 cell type
        } else {
            ind.group.f <- apply(ind.f, 1, any)
        }
        ind.group.f[is.na(ind.group.f)] <- FALSE
        return (tab[ind.group.f, ])
    }
    contrs <- c("InsulinResistant_versus_Basal", "MannitolCtrl_versus_Basal", "MannitolCtrl_versus_InsulinResistant")
    mclapply(contrs, mc.cores=length(contrs), function(cont) {
        res <- mclapply(enrfiles, mc.cores=length(enrfiles), function(x) {
            ontology <- strsplit(x, split='_')[[1]][2]
            tab.init <- read.table(x, header=T, sep='\t', quote='')
            
            tab.init <- tab.init[, c(1:3, 
                                     grep(cont, colnames(tab.init)))]
            tab.init <- tab.init[, c(1:3,
                                     grep("ECFC_RNAseq", colnames(tab.init)),
                                     grep("RNA_loadings", colnames(tab.init)),
                                     grep("ECFC_Prot", colnames(tab.init)),
                                     grep("Prot_loadings", colnames(tab.init)))]
            tab.init <- tab.init[, c(1:3, 
                                     grep("NES$", colnames(tab.init)), 
                                     grep("pvalue", colnames(tab.init)), 
                                     grep("qvalue", colnames(tab.init)), 
                                     grep("core_enrichment", colnames(tab.init)))]
            
            tab <- filterPathway(tab.init, filtering = list(p=grep("pvalue$", colnames(tab.init)),
                                                            q=grep("qvalues$", colnames(tab.init))),
                                 filtering.lim = c(p=PLIM, q=QLIM),
                                 groups = cells)
            rownames(tab) <- gsub(tab$ID, pattern="^GO", replacement=paste0("GO_", ontology))
            tab$ID <- gmt.id[[ontology]][tab$ID,1]
            
            gos_ss <- sapply(tab$ID, function(y) {
                sapply(tab$ID, function(z) {
                    if (exists(".GOSemSimEnv")) rm(".GOSemSimEnv")
                    goSim(y, z, semData=godatas[[ontology]], measure='Wang')
                })
            })
            colnames(gos_ss) <- rownames(tab)
            rownames(gos_ss) <- rownames(tab)
            gos_dist <- 1 - gos_ss
            gos_hc <- hclust(as.dist(gos_dist), method="ward.D")
            hcut <- 1.4
            clus <- cutree(gos_hc, h=hcut)
            clunum <- length(unique(clus))
            colors <- WEBSAFECOLORS[1:clunum]
            dend <- as.dendrogram(gos_hc) %>% set("branches_k_color", value=colors, k=clunum) #
            
            ##- MDS ----
            mds <- data.frame(cmdscale(gos_dist))
            colnames(mds) <- c("x", "y")
            mds$go <- rownames(mds)
            mds$col <- factor(get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))])
            p <- ggplot(mds, aes(x=x, y=y)) + 
                geom_point(aes(col=col), show.legend = F) + 
                geom_text(aes(label=go, colour=col), size=1.5, check_overlap = T, show.legend = F,
                          vjust=-1) +
                scale_colour_manual(values=levels(mds$col)) +
                xlim(-0.7, 0.7) +
                xlab("Coordinate 1") + ylab("Coordinate 2") +
                theme_bw()
            ggsave(p, filename=gsub(x, pattern=".tsv", replacement=paste0('_', cont, '_mds.pdf')), width=8, height=8)
            ##-----
            
            ##- clustering ----
            pdf(gsub(x, pattern=".tsv", replacement=ifelse(noclusterCell, '.pdf', paste0('_', cont, '_clusteredCell.pdf'))), 
                width=25, height=4+ifelse(nrow(tab)>200, nrow(tab)/5, ifelse(nrow(tab)>100, nrow(tab)/4, ifelse(nrow(tab)>25, nrow(tab)/3.5, 7))))
                    # switch(ontology,
                    #                     BP=nrow(tab)/3.6,
                    #                     CC=nrow(tab)/3,#4
                    #                     MF=nrow(tab)/2.1))
            hm.tab <- as.matrix(tab[, grep(".NES$", colnames(tab)), drop=F])
            hm.tab[is.na(hm.tab)] <- 0
            col.split <- sapply(strsplit(colnames(hm.tab), split="_"), function(hmcol) hmcol[1])
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern=paste0(paste(cells, collapse="|"), "_"), replacement = "")
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern=".NES", replacement = "")
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern=cont, replacement = "")
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern="^_|_$", replacement = "")
            hm <- Heatmap(hm.tab,
                          column_title=paste0("Selected ", titles[ontology], "\n enriched (p < ", PLIM, 
                                              ", q < ", QLIM, ") in all omics from both individual and integration analyses\n * p < 0.05, ** p < 0.01, *** p < 0.001"),
                          column_title_gp = gpar(fontsize=15, fontface='bold'),
                          column_split = if (noclusterCell) NULL else col.split,
                          heatmap_legend_param=list(title="Normalized Enrichement Score",
                                                    title_position="leftcenter-rot",
                                                    legend_direction="vertical",
                                                    legend_height=unit(2,'inches'),
                                                    title_gp=gpar(fontsize=font.size),
                                                    labels_gp=gpar(fontsize=font.size)),
                          cell_fun=function(j, i, x, y, width, height, col) {
                              grid.text(if (is.na(as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j])) ""
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 1e-3) "***"
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 1e-2) "**"
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 0.05) "*"
                                        else "", 
                                        x, y, gp=gpar(fontsize=font.size))
                          },
                          col=blueWhiteRed(50),
                          cluster_rows = dend,
                          cluster_columns = noclusterCell,
                          row_dend_width = unit(4, "cm"),
                          width = unit(20, "cm"),
                          show_column_names = T,
                          row_names_gp = gpar(col=get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))],
                                              fontsize=font.size),
                          row_names_max_width = max_text_width(tab$ID, #paste(rep(" ", max.length), collapse=''),
                                                               gp=gpar(fontsize=font.size)),
                          column_names_gp = gpar(col=COLORGROUPS[["BEAtDKDcell"]][col.split], fontsize=2*font.size)
            )
            draw(hm, 
                 heatmap_legend_side="left",
                 padding = unit(c(1, -8, 1, 0), "inches"),
                 newpage=F)
            dev.off()
            ##-----
            
            return (list(genes.clus=NULL, mds=mds))#genes.clus
        })
        
        names(res) <- gsub(enrfiles, pattern='_withOPLS_stat.tsv', replacement = '')
        #save(res, file='goclustering.RData')
    })
} else if (PROJECT %in% c("Endothelial")) {
    ############################################
    ## GO clustering
    setwd(paste0(DIR, '/limma/output/', PROJECT, '/pathway_combined'))
    #setwd('~/Work/BEAt-DKD/output/Endothelial/pathway_combined/')
    enrfiles <- grep("^GO_.+withOPLS_stat.tsv", list.files(), value=T)
    FUNCTIONASSAY <- TRUE
    
    #' param filtering List of index vectors of the same length
    #' param filtering.lim Vector of limits, with the same length as filtering vectors
    filterPathway <- function(tab, filtering = NULL, filtering.lim = NULL, groups = NULL) {
        if (is.null(filtering)) return (tab)
        colnames.1 <- gsub(colnames(tab)[filtering[[1]]], pattern=".pvalue$|.qvalues$|.NES$", replacement="")
        ## check if filtering is correctly formatted
        if (length(filtering) > 1) {
            correct.order <- all(sapply(2:length(filtering), function(i) {
                colnames.i <- gsub(colnames(tab)[filtering[[i]]], pattern=".pvalue$|.qvalues$|.NES$", replacement="")
                return (all(colnames.i==colnames.1))
            }))
            stopifnot(correct.order)
        }
        ## filtering on each row of tab
        ind.f <- do.call(cbind, mclapply(1:length(filtering[[1]]), mc.cores=length(filtering[[1]]), function(j) {
            apply(do.call(cbind, mclapply(1:length(filtering), mc.cores=length(filtering), function(i) {
                ifelse(is.na(tab[, filtering[[i]][j]]), 
                       NA, 
                       as.numeric(tab[, filtering[[i]][j]]) <= filtering.lim[i])
            })), 1, all)
        }))
        colnames(ind.f) <- colnames.1
        if (!is.null(groups)) {
            ind.group.f <- apply(do.call(cbind, mclapply(groups, mc.cores=length(groups), function(gr) {
                gr.cont <- unique(sub(sub(grep(gr, colnames(ind.f), value=T),
                                          pattern=paste0(gr, "_"), replacement=""),
                                      pattern="_RNAseq|_RNA_loadings|_Prot|_Prot_loadings", replacement=""))
                gr.pass <- t(apply(ind.f[, grep(gr, colnames(ind.f)), drop=F], 1, function(x) {
                    sapply(gr.cont, function(grc) {
                        ## more than 3/4 filtering multiblock/individual satisfactions
                        #print(grep(paste0(gr, "_", grc), colnames(ind.f), value=T))
                        sum(na.omit(x[grep(paste0(gr, "_", grc), colnames(ind.f), value=T)]))/length(grep(grc, names(x))) > 3/4
                    })
                }))
                apply(gr.pass, 1, function(y) any(y[match(c("DiabeticSoup_versus_Basal", 
                                                            "MannitolCtrl_versus_Basal", 
                                                            "DiabeticSoup_versus_MannitolCtrl"), colnames(gr.pass))])) # in at least 1 contrast
                # apply(ind.f[, grep(gr, colnames(ind.f)), drop=F], 1, function(x) {
                #     ## more than 3/4 filtering multiblock/individual satisfactions
                #     sum(na.omit(x[grep(gr, colnames(ind.f), value=T)]))/length(x) > 3/4
                # })
            })), 1, function(x) any(x[match(c("GEnC","ECFC"), groups)])) # in at least 1 cell type of GEnC or ECFC
        } else {
            ind.group.f <- apply(ind.f, 1, any)
        }
        ind.group.f[is.na(ind.group.f)] <- FALSE
        return (tab[ind.group.f, , drop=F])
    }
    # contrs <- c("DiabeticSoup_versus_Basal", "MannitolCtrl_versus_Basal", "DiabeticSoup_versus_MannitolCtrl",
    #             "DiabeticPatient_versus_Basal", "DiabeticPatient_versus_DiabeticSoup", "DiabeticPatient_versus_MannitolCtrl")
    # mclapply(contrs, mc.cores=length(contrs), function(cont) {
        res <- lapply(enrfiles, function(x) {
            ##- clean data ----
            ontology <- strsplit(x, split='_')[[1]][2]
            tab.init <- read.table(x, header=T, sep='\t', quote='')
            
            # tab.init <- tab.init[, c(1:3, 
            #                          grep(cont, colnames(tab.init)))]
            tab.init <- tab.init[, c(1:3,
                                     grep("RNAseq\\.", colnames(tab.init)),
                                     grep("RNA_loadings", colnames(tab.init)),
                                     grep("Prot\\.", colnames(tab.init)),
                                     grep("Prot_loadings", colnames(tab.init)))]
            tab.init <- tab.init[, c(1:3, 
                                     grep("GEnC", colnames(tab.init)), 
                                     grep("ECFC", colnames(tab.init)), 
                                     grep("HMEC", colnames(tab.init)),
                                     grep("RAEC", colnames(tab.init))
                                     )]
            tab.init <- tab.init[, c(1:3, 
                                     grep("NES$", colnames(tab.init)), 
                                     grep("pvalue", colnames(tab.init)), 
                                     grep("qvalue", colnames(tab.init)), 
                                     grep("core_enrichment", colnames(tab.init)))]
            rownames(tab.init) <- gsub(tab.init$ID, pattern="^GO", replacement=paste0("GO_", ontology))
            tab.init$ID <- gmt.id[[ontology]][tab.init$ID,1]
            ##-----
            
            ##- Functional assay ----
            if (FUNCTIONASSAY && ontology=="BP") {
                functionalGO <- c('GO:0001935', 'GO:0043542', 'GO:0001525')
                #Ontology(functionalGO)
                GOoffspring <- as.list(GOBPOFFSPRING)
                GOoffspring <- GOoffspring[!is.na(GOoffspring)]
                #functionalGOoffspring <- lapply(GOoffspring[functionalGO], function(x) {Term((x))})
                tab <- tab.init[tab.init$ID %in% unname(unlist(GOoffspring[functionalGO])),]
                write.table(tab, file="GO_BP_Endothelial_functionalassay_withOPLS_stat.tsv", sep='\t', quote=F, col.names=NA)
            } else {
                tab <- filterPathway(tab.init, filtering = list(p=grep("pvalue$", colnames(tab.init)),
                                                                q=grep("qvalues$", colnames(tab.init))),
                                     filtering.lim = c(p=PLIM, q=QLIM),
                                     groups = cells)
                # rownames(tab) <- gsub(tab$ID, pattern="^GO", replacement=paste0("GO_", ontology))
                # tab$ID <- gmt.id[[ontology]][tab$ID,1]
            }
            ##-----
            
            ##- GO similarity ----
            gos_ss <- sapply(tab$ID, function(y) {
                sapply(tab$ID, function(z) {
                    if (exists(".GOSemSimEnv")) rm(".GOSemSimEnv")
                    goSim(y, z, semData=godatas[[ontology]], measure='Wang')
                })
            })
            colnames(gos_ss) <- rownames(tab)
            rownames(gos_ss) <- rownames(tab)
            gos_dist <- 1 - gos_ss
            gos_hc <- hclust(as.dist(gos_dist), method="ward.D")
            hcut <- 1.6
            clus <- cutree(gos_hc, h=hcut)
            clunum <- length(unique(clus))
            colors <- WEBSAFECOLORS[1:clunum]
            dend <- as.dendrogram(gos_hc) %>% set("branches_k_color", value=colors, k=clunum) #
            ##-----
            
            ##- MDS ----
            mds <- data.frame(cmdscale(gos_dist))
            colnames(mds) <- c("x", "y")
            mds$go <- rownames(mds)
            mds$col <- factor(get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))])
            p <- ggplot(mds, aes(x=x, y=y)) + 
                geom_point(aes(col=col), show.legend = F) + 
                geom_text(aes(label=go, colour=col), size=1.5, check_overlap = T, show.legend = F,
                          vjust=-1) +
                scale_colour_manual(values=levels(mds$col)) +
                xlim(-0.7, 0.7) +
                xlab("Coordinate 1") + ylab("Coordinate 2") +
                theme_bw()
            ggsave(p, filename=gsub(x, pattern=".tsv", replacement=paste0('_', ifelse(FUNCTIONASSAY, "functionalassay", ""), '_mds.pdf')), 
                   width=8, height=8)
            ##-----
            
            ##- clustering ----
            pdf(gsub(x, 
                     pattern=".tsv", 
                     replacement=ifelse(noclusterCell, 
                                        paste0('_', ifelse(FUNCTIONASSAY, "functionalassay", ""), '.pdf'), 
                                        paste0('_', ifelse(FUNCTIONASSAY, "functionalassay", ""), '_clusteredCell.pdf'))), 
                width=31, #25 
                height=4+ifelse(nrow(tab)>200, 
                                nrow(tab)/5, 
                                ifelse(nrow(tab)>100, 
                                       nrow(tab)/4, 
                                       ifelse(nrow(tab)>25, 
                                              nrow(tab)/3.5, 
                                              7))))
            hm.tab <- as.matrix(tab[, grep(".NES$", colnames(tab)), drop=F])
            hm.tab[is.na(hm.tab)] <- 0
            col.split <- sapply(strsplit(colnames(hm.tab), split="_"), function(hmcol) hmcol[1])
            #colnames(hm.tab) <- sub(colnames(hm.tab), pattern=paste(paste0(cells, "_"), collapse="|"), replacement = "")
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern=".NES", replacement = "")
            #colnames(hm.tab) <- sub(colnames(hm.tab), pattern=cont, replacement = "")
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern="^_|_$", replacement = "")
            colnames(hm.tab) <- sub(colnames(hm.tab), pattern="__", replacement = "_")
            hm <- Heatmap(hm.tab,
                          column_title=if (FUNCTIONASSAY) {
                              paste0(titles[ontology],
                                     "\n * p < 0.05, ** p < 0.01, *** p < 0.001")
                          } else {
                              paste0("Selected ", 
                                     titles[ontology], 
                                     "\n enriched (p < ", PLIM, 
                                     ", q < ", QLIM, ") in all omics from both individual and integration analyses\n * p < 0.05, ** p < 0.01, *** p < 0.001")
                              },
                          column_title_gp = gpar(fontsize=15, fontface='bold'),
                          column_split = if (noclusterCell) NULL else col.split,
                          heatmap_legend_param=list(title="Normalized Enrichement Score",
                                                    title_position="leftcenter-rot",
                                                    legend_direction="vertical",
                                                    legend_height=unit(2,'inches'),
                                                    title_gp=gpar(fontsize=font.size),
                                                    labels_gp=gpar(fontsize=font.size)),
                          cell_fun=function(j, i, x, y, width, height, col) {
                              grid.text(if (is.na(as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j])) ""
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 1e-3) "***"
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 1e-2) "**"
                                        else if (as.matrix(tab[, grep(".pvalue$", colnames(tab)), drop=F])[i, j] < 0.05) "*"
                                        else "", 
                                        x, y, gp=gpar(fontsize=font.size))
                          },
                          col=blueWhiteRed(50),
                          cluster_rows = dend,
                          cluster_columns = noclusterCell,
                          row_dend_width = unit(4, "cm"),
                          width = unit(20, "cm"),
                          show_column_names = T,
                          row_names_gp = gpar(col=get_leaves_branches_col(dend)[match(1:nrow(tab), order.dendrogram(dend))],
                                              fontsize=font.size),
                          row_names_max_width = max_text_width(tab$ID, #paste(rep(" ", max.length), collapse=''),
                                                               gp=gpar(fontsize=font.size)),
                          column_names_rot = 45,
                          column_names_gp = gpar(col=COLORGROUPS[["BEAtDKDcell"]][col.split], fontsize=1*font.size)
            )
            draw(hm, 
                 heatmap_legend_side="left",
                 padding = unit(c(1, -8, 1, 0), "inches"),
                 #padding = unit(c(1, -2, 1, 6), "inches"),
                 newpage=F)
            dev.off()
            ##-----
            
            return (list(genes.clus=NULL, mds=mds))#genes.clus
        })
        
        names(res) <- gsub(enrfiles, pattern='_withOPLS_stat.tsv', replacement = '')
        #save(res, file='goclustering.RData')
    # })
}



