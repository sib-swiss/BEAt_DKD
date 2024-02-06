##- settings ----
rm(list=ls())
LIBRARIES <- c("clusterProfiler",
               "ggplot2", "reshape2", "gridExtra", "gtable", "ggpubr", "ggVennDiagram",
               "assertthat", "parallel"
)
invisible(sapply(LIBRARIES, function(lib) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
}))
options(stringsAsFactors=F)

#DIR <- '/db/scratch/ttran/BEAt_DKD/limma'
DIR <- '~/Work/BEAtDKD_analysis'
ADJUST.METHOD <- 'BH'
GENEID <- file.path(DIR, "genes.id.txt")
FDR <- 0.05
LOGFC <- 1
GMTFILES <- c(HALLMARK="h.all.v7.1.symbols.gmt",
              GO_BP="c5.bp.v7.1.symbols.gmt",
              GO_CC="c5.cc.v7.1.symbols.gmt",
              GO_MF="c5.mf.v7.1.symbols.gmt",
              KEGG="c2.cp.kegg.v7.1.symbols.gmt",
              REACTOME="c2.cp.reactome.v7.1.symbols.gmt",
              PID="c2.cp.pid.v7.1.symbols.gmt",
              BIOCARTA="c2.cp.biocarta.v7.1.symbols.gmt",
              CGP="c2.cgp.v7.1.symbols.gmt",
              IMMUNO="c7.all.v7.1.symbols.gmt"
)
# GMTFILES.miRNA <- c(GO_BP="miRNA/GO_BP_validated_miRTarBase_all.gmt",
#                     GO_CC="miRNA/GO_CC_validated_miRTarBase_all.gmt",
#                     GO_MF="miRNA/GO_MF_validated_miRTarBase_all.gmt",
#                     KEGG="miRNA/KEGG_validated_miRTarBase_all.gmt",
#                     REACTOME="miRNA/REACTOME_validated_miRTarBase_all.gmt"
# )
# CLUSTERPROT  <- "9606.clusters.proteins.v11.0.txt"  # STRINGdb protein cluster
# CLUSTERANNO  <- "9606.clusters.info.v11.0.txt"      # STRINGdb protein cluster annotation
# CLUSTERRELT  <- "9606.clusters.tree.v11.0.txt"      # STRINGdb protein cluster relation
# ENTREZSTRING <- "human.entrez_2_string.2018.tsv"    # entrez - string mapping
# SYMBOLSTRING <- "human.name_2_string.tsv"           # symbol - string mapping


PROJECT <- "Bristol"
setwd(file.path(DIR, "output", PROJECT))
COLORGROUPS <- list()
COLORGROUPS[["BEAtDKDcell"]] <- setNames(
    c("coral", "gray40", "turquoise3", "gold2"),
    c("ABN",   "GEnC",   "K29",        "PT34"))
COLORGROUPS[["BEAtDKDcond"]] <- setNames(
    c("green2", "magenta2",         "darkgreen", "darkmagenta"),
    c("Basal",  "InsulinResistant", "IR_Basal",  "IR_InsulinResistant")
)
OMICS <- switch (PROJECT,
                 Bristol = c("RNAseq", "Prot")
)
CELLS <- switch(PROJECT,
                Bristol = c("ABN", "GEnC", "K29", "PT34")
)
CELLS.display <- setNames(c("Pod", "GEC", "MC", "PTC"), CELLS)
CELLS.pch <- setNames(c(16, 15, 18, 17), CELLS)
IR <- switch(PROJECT,
             Bristol = c("IR", "")
)
dirs <- intersect(as.vector(outer(as.vector(outer(OMICS, CELLS, paste0)), IR, paste0)), list.files())
dirs.all <- c("RNAseqallall", "Protallall")
##-----

##- gene id mapping ----
idmap <- read.table(GENEID, header=T, sep='\t', quote="")
##-----

##- multiple-testing correction on celltype-specific limma DE analyses ----
## correct by data type: RNA-seq, proteomics, miRNA-seq
## bypass for Lund with only resequenced ECFC samples
pvalues <- mclapply(dirs, mc.cores=max(length(dirs), detectCores()), function(x) {
    tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="")
    return (tab[, grep("^p_value_InsulinResistant_versus_Basal", colnames(tab))])
})
names(pvalues) <- dirs
padjust <- unlist(lapply(c("^RNAseq", "^miRNAseq", "^Prot"), function(x) {
    p.adjust(unlist(pvalues[grep(x, names(pvalues))]), method=ADJUST.METHOD)
}))
table(sub(names(padjust), pattern='[0-9]+$', replacement=''))
# write obtained padjust to DE.tsv
invisible(lapply(dirs, function(x) {
    tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
    tab[, grep(paste0("^FDR", "_InsulinResistant_versus_Basal"), 
               colnames(tab))] <- matrix(
                   padjust[grep(paste0("^", x, "[0-9]"), names(padjust))], 
                   ncol=length(grep(paste0("^FDR", "_InsulinResistant_versus_Basal"), colnames(tab))), 
                   byrow=F)
    write.table(tab, file=paste0(x, '/DE.tsv'), sep='\t', quote=F, col.names=NA) #TODO uncomment
}))
##-----

##- multiple-testing correction on celltype-specific limma DE analyses with different types of correction method ----
## correct by data type: RNA-seq, proteomics, miRNA-seq
## specific analysis for the check on different types of p-value correction
if (PROJECT=='Bristol.FDR') {
    pvalues <- mclapply(dirs, mc.cores=max(length(dirs), detectCores()), function(x) {
        tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="")
        if (PROJECT %in% c('Lund', 'Endothelial')) return (tab[, grep("^p_value_", colnames(tab))])
        else return (tab[, grep("^p_value_InsulinResistant_versus_Basal", colnames(tab))])
    })
    names(pvalues) <- gsub(dirs, pattern="^Bristol", replacement="RNAseq")
    padjust <- do.call(cbind, mclapply(c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"), mc.cores=6, function(adjmet) {
        unlist(mclapply(c("^RNAseq", "^miRNAseq", "^Prot"), mc.cores=3, function(x) {
            p.adjust(unlist(pvalues[grep(x, names(pvalues))]), method=adjmet)
        }))
    }))
    colnames(padjust) <- paste0(c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"), '_adjusted_allCells')
    
    invisible(mclapply(dirs, mc.cores=length(dirs), function(x) {
        tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
        tab <- cbind(tab, 
                     padjust[grep(
                         paste0("^", gsub(x, pattern="^Bristol|^Lund", replacement="RNAseq"), 
                                ifelse(!(PROJECT %in% c("Lund", "Endothelial")), "[0-9]", "")), rownames(padjust)),])
        write.table(tab, file=paste0(x, '/DEs.tsv'), sep='\t', quote=F, col.names=NA)
    }))
    
    prot.mito <- read.csv(paste0(DIR, '../data/raw/5Sept22_mitoproteins.csv'), header=T)
    pvalues.mito <- mclapply(grep("^Prot", dirs, value=T), mc.cores=max(length(dirs), detectCores()), function(x) {
        tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
        tab <- tab[prot.mito$UniProtKB,]
        if (PROJECT %in% c('Lund', 'Endothelial')) return (tab[, grep("^p_value_", colnames(tab))])
        else return (tab[, grep("^p_value_InsulinResistant_versus_Basal", colnames(tab))])
    })
    names(pvalues.mito) <- grep("^Prot", dirs, value=T)
    
    padjust.mito <- do.call(cbind, mclapply(c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"), mc.cores=6, function(adjmet) {
        p.adjust(unlist(pvalues.mito), method=adjmet)
    }))
    colnames(padjust.mito) <- paste0(c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"), '_adjusted_mitochondria')
    
    invisible(mclapply(grep("^Prot", dirs, value=T), mc.cores=length(dirs), function(x) {
        tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
        tab <- tab[prot.mito$UniProtKB,]
        tab <- cbind(tab, 
                     padjust.mito[grep(
                         paste0("^", x, 
                                ifelse(!(PROJECT %in% c("Lund", "Endothelial")), "[0-9]", "")), rownames(padjust.mito)),])
        write.table(tab, file=paste0(x, '/DEs_mito.tsv'), sep='\t', quote=F, col.names=NA)
    }))
}
##-----

##- PCA ----
omics.pca.plots.all <- lapply(dirs.all, function(x) {
    tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
    counts <- tab[,!grepl("^FDR_|^p_value|^UpSetGroup|^name|^Log2FC_", colnames(tab))]
    
    titlex <- ifelse(x=="Protallall", "Proteomics", "RNA-seq")
    constant.idx <- which(apply(counts, 1, function(xx) length(unique(xx))==1))
    if (length(constant.idx) > 0) counts[constant.idx,] <- counts[constant.idx,]*(1+(1e-3)*rnorm(ncol(counts), mean=0, sd=1))
    # PCA
    pca <- prcomp(t(counts), center=T, scale.=T)
    # PCs for plotting
    pca.dims <- list(c(1,2), c(1,3))
    
    pca.df <- data.frame(pca$x[, sort(unique(unlist(pca.dims)))])
    pca.df$cell <- rownames(pca.df) %>% gsub("_.+", "", .) %>% CELLS.display[.]
    pca.df$group <- rownames(pca.df) %>% 
        gsub("_ABN.+|_K29.+|_GEnC.+|_PT34.+|_GTG.+", "", .) %>%
        gsub("ABN_|K29_|GEnC_|PT34_", "", .) %>%
        gsub("^_", "", .)
    pca.df$color <- COLORGROUPS$BEAtDKDcond[pca.df$group]
        
    # for all cell types
    pca.plot <- ggplot(data=pca.df) +
        theme_bw() +
        theme(
            title            = element_text(size=18, face="bold"),
            axis.text        = element_text(size=18, face="plain"),
            axis.title       = element_text(size=18, face="plain"),
            plot.margin      = margin(0,0.5,0.5,0.5, "cm"),
            legend.title     = element_text(size=18, face="bold"),
            legend.text      = element_text(size=18, face="plain"),
            legend.direction = "horizontal",
            legend.position  = "top") + 
        ylim(-120, 120)
    pca.plots.all <- lapply(1L:length(pca.dims), function(i) {
        pc1 <- pca.df[,paste0('PC', pca.dims[[i]][1])]
        pc2 <- pca.df[,paste0('PC', pca.dims[[i]][2])]
        pc1 <- enquo(pc1)
        pc2 <- enquo(pc2)
        p <- pca.plot +
            ggtitle(titlex) +
            geom_point(aes(x=!!pc1,
                           y=!!pc2,
                           colour=group, shape=cell), size=4) +
            xlab(paste0("PC",
                        pca.dims[[i]][1],
                        " (",
                        round((summary(pca)$importance)[2, pca.dims[[i]][1]], 2), ")")) +
            ylab(paste0("PC",
                        pca.dims[[i]][2],
                        " (",
                        round((summary(pca)$importance)[2, pca.dims[[i]][2]], 2), ")")) +
            scale_colour_manual(name="Condition", 
                                values=COLORGROUPS$BEAtDKDcond[as.vector(unique(pca.df$group))],
                                guide=guide_legend(title.position='top')) +
            scale_shape_manual(name="Cell", values=unname(CELLS.pch[names(sort(CELLS.display))]),
                               guide=guide_legend(title.position='top')) +
            stat_ellipse(aes(x=!!pc1,
                             y=!!pc2,
                             colour=cell))
        return (p)
    })
    
    return (pca.plots.all)
})
ag <- ggarrange(plotlist=c(unlist(omics.pca.plots.all, recursive = F))[c(1,3,2,4)], 
                # labels=c("a", "b", rep("", 2)),
                # font.label=list(size=40, face='bold'),
                # label.y=1.02,
                common.legend=T, ncol=2, nrow=2)
ggsave(ag, filename="Fig2AB.pdf", width=12, height=12)

proteomics.info <- t(read.table("db/proteomics/Scaled.tsv", header=F, sep='\t', quote="", nrows=1)[1, -1])
rownames(proteomics.info) <- proteomics.info %>% gsub("^F._[0-9A-Z]+_", "", .)
rownames(proteomics.info) <- rownames(proteomics.info) %>% paste0(., "_",.) %>% gsub("_[0-9]+_", "_", .)
proteomics.info <- data.frame(sample=proteomics.info[,1], 
                              batch=paste0("batch", as.numeric(gsub("_.+|^F", "", proteomics.info[,1])) %% 2 + 1))

omics.pca.plots.each <- lapply(dirs.all, function(x) {
    tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
    counts <- tab[,!grepl("^FDR_|^p_value|^UpSetGroup|^name|^Log2FC_", colnames(tab))]
    
    titlex <- ifelse(x=="Protallall", "Proteomics", "RNA-seq")
    constant.idx <- which(apply(counts, 1, function(xx) length(unique(xx))==1))
    if (length(constant.idx) > 0) counts[constant.idx,] <- counts[constant.idx,]*(1+(1e-3)*rnorm(ncol(counts), mean=0, sd=1))
    
    # for each cell type
    pca.plots.each <- lapply(names(sort(CELLS.display)), function(cel) {
        # PCA
        pca <- prcomp(t(counts[, grep(cel, colnames(counts))]), center=T, scale.=T)
        # PCs for plotting
        pca.dims <- list(c(1,2))
        
        pca.df <- data.frame(pca$x[, sort(unique(unlist(pca.dims)))])
        pca.df$cell <- rownames(pca.df) %>% gsub("_.+", "", .) %>% CELLS.display[.]
        pca.df$group <- rownames(pca.df) %>% 
            gsub("_ABN.+|_K29.+|_GEnC.+|_PT34.+|_GTG.+", "", .) %>%
            gsub("ABN_|K29_|GEnC_|PT34_", "", .) %>%
            gsub("^_", "", .)
        pca.df$color <- COLORGROUPS$BEAtDKDcond[pca.df$group]
        pca.df$batch <- if (x=="Protallall") proteomics.info[rownames(pca.df), 'batch'] else 16
            
        pca.plot <- ggplot(pca.df) +
            theme_bw() +
            theme(
                title            = element_text(size=18, face="bold"),
                axis.text        = element_text(size=18, face="plain"),
                axis.title       = element_text(size=18, face="plain"),
                plot.margin      = margin(0,0.5,0.5,0.5, "cm"),
                legend.title     = element_text(size=18, face="bold"),
                legend.text      = element_text(size=18, face="plain"),
                legend.direction = "horizontal",
                legend.position  = "top")
        pca.plots <- lapply(1L:length(pca.dims), function(i) {
            p <- pca.plot +
                ggtitle(paste0(titlex, ": ", CELLS.display[cel]))
            if (x=="Protallall")
                p <- p + 
                    geom_point(aes(x=PC1,
                                   y=PC2,
                                   colour=group, shape=batch), size=5)#, shape=CELLS.pch[cel]) +
            else           
                p <- p + 
                    geom_point(aes(x=PC1,
                                   y=PC2,
                                   colour=group), shape=16, size=5)
            p <- p +     
                xlab(paste0("PC",
                            pca.dims[[i]][1],
                            " (",
                            round((summary(pca)$importance)[2, pca.dims[[i]][1]], 2), ")")) +
                ylab(paste0("PC",
                            pca.dims[[i]][2],
                            " (",
                            round((summary(pca)$importance)[2, pca.dims[[i]][2]], 2), ")")) +
                scale_colour_manual(name="Condition", 
                                    values=COLORGROUPS$BEAtDKDcond[as.vector(unique(pca.df$group))],
                                    guide=guide_legend(title.position='top')) +
                stat_ellipse(data=pca.df %>% filter(group=="IR_InsulinResistant"),
                             type='norm',
                             level=0.6,
                             aes(x=PC1,
                                 y=PC2,
                                 colour=group))
            if (x=="Protallall")
                p <- p + scale_shape_manual(name="Batch", values=c(15, 17),
                                            guide=guide_legend(title.position='top'))
            return (p)
        })
        return (pca.plots)
    })
    
    return (pca.plots.each)
})
ag <- ggarrange(plotlist=unlist(unlist(omics.pca.plots.each, recursive = F), recursive = F),
                # labels=c("c", rep("", 3), "d", rep("", 3)),
                # font.label=list(size=40, face='bold'),
                # label.y=1.1,
                common.legend=T, legend.grob=get_legend(omics.pca.plots.each[[2]][[1]]),
                ncol=4, nrow=2)
ggsave(ag, filename="SubFig2CD.pdf", width=20, height=10)
##-----

##- number of DE genes per contrast ----
regulated.df <- do.call(rbind, lapply(dirs, function(x) {
    tab <- read.table(paste0(x, '/DE.tsv'), header=T, sep='\t', quote="")
    const <- sub(grep("^FDR_", colnames(tab), value=T), pattern='^FDR_', replacement='')

    df <- melt(setNames(lapply(c(0, 1), function(lfc) {
        lfc.df <- t(sapply(const, function(cnt) {
            fc.cnt <- paste0("Log2FC_", cnt)
            fdr.cnt <- paste0("FDR_", cnt)
            c(Down=sum(tab[, fc.cnt] < -lfc & tab[, fdr.cnt] < FDR),
              Up  =sum(tab[, fc.cnt] >  lfc & tab[, fdr.cnt] < FDR))
        }))
        rownames(lfc.df) <- sub(rownames(lfc.df), pattern='_IR', replacement='')
        rownames(lfc.df) <- paste0(gsub(x, pattern=paste(CELLS, collapse='|'), replacement=''), '_', rownames(lfc.df))
        return (lfc.df)
    }), c("Log2FC > 0", "Log2FC > 1")))

    return (df)
})) #TO REMOVE

regulated.df <- do.call(rbind, lapply(dirs.all, function(x) {
    tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="")
    const <- sub(grep("^FDR_", colnames(tab), value=T), pattern='^FDR_', replacement='')
    df <- melt(setNames(lapply(c(0, 1), function(lfc) {
        lfc.df <- t(sapply(const, function(cnt) {
            fc.cnt <- paste0("Log2FC_", cnt)
            fdr.cnt <- paste0("FDR_", cnt)
            c(Down=sum(tab[, fc.cnt] < -lfc & tab[, fdr.cnt] < FDR),
              Up  =sum(tab[, fc.cnt] >  lfc & tab[, fdr.cnt] < FDR))
        }))
        rownames(lfc.df) <- paste0(gsub(x, pattern=paste(CELLS, collapse='|'), replacement=''), '_', rownames(lfc.df))
        return (lfc.df)
    }), c("Log2FC > 0", "Log2FC > 1")))
    
    return (df)
}))

colnames(regulated.df)[colnames(regulated.df)=="Var1"] <- "Contrast"
colnames(regulated.df)[colnames(regulated.df)=="Var2"] <- "Regulation"
colnames(regulated.df)[colnames(regulated.df)=="value"] <- "Size"
colnames(regulated.df)[colnames(regulated.df)=="L1"] <- "Cutoff"

for (cel in CELLS) regulated.df$Contrast <- regulated.df$Contrast %>% sub(cel, CELLS.display[cel] , .)
regulated.df$Cutoff <- factor(regulated.df$Cutoff,
                              levels=sort(unique(regulated.df$Cutoff)))
invisible(assert_that(!any(paste0(regulated.df[regulated.df$Cutoff=="Log2FC > 1", ]$Contrast,
                                  regulated.df[regulated.df$Cutoff=="Log2FC > 1", ]$Regulation)!=
                               paste0(regulated.df[regulated.df$Cutoff=="Log2FC > 0", ]$Contrast,
                                      regulated.df[regulated.df$Cutoff=="Log2FC > 0", ]$Regulation)),
                      msg="Error in melt!"))
regulated.df$Stack <- regulated.df$Size * (regulated.df$Cutoff=="Log2FC > 1")
regulated.df[regulated.df$Cutoff=="Log2FC > 0",]$Stack <- regulated.df[regulated.df$Cutoff=="Log2FC > 0",]$Size - 
    regulated.df[regulated.df$Cutoff=="Log2FC > 1",]$Size

regulated.df$Cell <- regulated.df$Contrast %>% gsub(".+_in_", "", .) %>% gsub("_.+", "", .)
regulated.df$IR <- ifelse(grepl("IR_|_IR", regulated.df$Contrast), "IR", "")
regulated.df$Omics <- regulated.df$Contrast %>% gsub("allall_.+|IR_.+|_.+", "", .) %>% factor(., levels=c("RNAseq", "Prot"))
regulated.df$Cond <- paste0(regulated.df$Cell, ifelse(regulated.df$IR=="IR", "_", ""), regulated.df$IR)
regulated.df$Cond <- factor(regulated.df$Cond, levels=rev(sort(unique(regulated.df$Cond))))

colors <- c("Down.Log2FC > 0" = "#92C5DE",
            "Up.Log2FC > 0"   = "#F4A582",
            "Down.Log2FC > 1" = "#0571B0",
            "Up.Log2FC > 1"   = "#CA0020")
colors <- colors[levels(interaction(regulated.df$Regulation, regulated.df$Cutoff))]
regulated.df$Color <- factor(ifelse(
    regulated.df$Regulation=="Down",
    ifelse(regulated.df$Cutoff=="Log2FC > 0", colors["Down.Log2FC > 0"], colors["Down.Log2FC > 1"]),
    ifelse(regulated.df$Cutoff=="Log2FC > 0", colors["Up.Log2FC > 0"], colors["Up.Log2FC > 1"])
), levels=colors)
regulated.df[regulated.df$Regulation=="Down",]$Stack <- -1 * regulated.df[regulated.df$Regulation=="Down",]$Stack
names(colors) <- names(colors) %>% 
    sub("Down.Log2FC > 0", "Log2FC < 0", .) %>%
    sub("Down.Log2FC > 1", "Log2FC < -1", .) %>%
    sub("Up.Log2FC >", "Log2FC >", .)
# plot
p <- ggplot(regulated.df, aes(x=Cond, y=Stack, fill=interaction(Regulation, Cutoff))) +
    facet_wrap(~Omics, ncol=2, scales='free',
               labeller=labeller(Omics=setNames(c('Number of DE proteins', "Number of DE transcripts"), 
                                                c("Prot", "RNAseq")))) +
    xlab("") + ylab("") +
    ylim(-max(regulated.df$Size), max(regulated.df$Size)) +
    coord_flip() +
    #ggtitle(paste0("InsulinResistant versus Basal")) +
    theme_bw() +
    theme(title           = element_text(size=12, face='bold'),
          legend.text     = element_text(size=16),
          legend.position = "top",
          strip.text      = element_text(size=18),
          plot.margin     = margin(0,1,0,0, "cm"),
          axis.text       = element_text(size=16),
          axis.text.y     = element_text(size=16)
    )
p <- p + geom_bar(width=0.6, position="stack", stat="identity") +
    geom_hline(yintercept=0, color="white") +
    scale_fill_manual(name="", values=unname(colors), labels=names(colors))
ggsave(p, filename='Fig2CD.pdf', width=12, height=6.1)
##-----

##- Volcano plots ----
volcano.list <- unlist(setNames(mclapply(dirs, mc.cores=length(dirs), function(x) {
    tab <- read.table(paste0(x, '/DE.tsv'), header=T, sep='\t', quote="", row.names=1)
    const <- sub(grep("^FDR_", colnames(tab), value=T), pattern='^FDR_', replacement='')
    setNames(lapply(const, function(cnt) {
        fc.cnt <- paste0("Log2FC_", cnt)
        fdr.cnt <- paste0("FDR_", cnt)
        cnt <- sub(cnt, pattern="_IR", replacement="")
        mat <- cbind(paste0(gsub(x, pattern=paste(CELLS, collapse='|'), replacement=''), '_', cnt), 
                     tab[, c('name', fc.cnt, fdr.cnt)])
        rownames(mat) <- rownames(tab)
        colnames(mat) <- c("Contrast", 'ID', "Log2FC", "FDR")
        return (mat)
    }), const)
}), dirs), recursive = F)

volcano.list <- unlist(setNames(mclapply(dirs.all, mc.cores=length(dirs), function(x) {
    tab <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="", row.names=1)
    const <- sub(grep("^FDR_", colnames(tab), value=T), pattern='^FDR_', replacement='')
    setNames(lapply(const, function(cnt) {
        fc.cnt <- paste0("Log2FC_", cnt)
        fdr.cnt <- paste0("FDR_", cnt)
        mat <- cbind(paste0(gsub(x, pattern=paste(CELLS, collapse='|'), replacement=''), '_', cnt), 
                     tab[, c('name', fc.cnt, fdr.cnt)])
        rownames(mat) <- rownames(tab)
        colnames(mat) <- c("Contrast", 'ID', "Log2FC", "FDR")
        return (mat)
    }), const)
}), dirs.all), recursive = F)
names(volcano.list) <- NULL
volcano.df <- do.call(rbind, volcano.list)
for (cel in CELLS) volcano.df$Contrast <- volcano.df$Contrast %>% sub(cel, CELLS.display[cel] , .)

volcano.df$Cell <- volcano.df$Contrast %>% gsub(".+_in_", "", .) %>% gsub("_.+", "", .)
volcano.df$IR <- ifelse(grepl("IR_|_IR", volcano.df$Contrast), "IR", "")
volcano.df$Omics <- volcano.df$Contrast %>% 
    gsub("allall_.+|IR_.+|_.+", "", .) %>% 
    factor(., levels=c("RNAseq", "Prot"))
volcano.df$Cond <- paste0(volcano.df$Cell, ifelse(volcano.df$IR=="IR", "_", ""), volcano.df$IR)
volcano.df$Cond <- factor(volcano.df$Cond, levels=sort(unique(volcano.df$Cond)))

volcano.df$sig <- factor(ifelse(volcano.df$FDR < FDR,
                                ifelse(volcano.df$Log2FC > 1, 
                                       colors["Log2FC > 1"],
                                       ifelse(volcano.df$Log2FC < -1, 
                                              colors["Log2FC < -1"],
                                              ifelse(volcano.df$Log2FC < 0, 
                                                     colors["Log2FC < 0"],
                                                     colors["Log2FC > 0"]))),
                                "grey70"),
                         levels=c(colors, "grey70"))

# Fig2F volcano plot for RNA-seq
volcano.df.plot <- volcano.df[volcano.df$IR=="IR" & volcano.df$Omics=="RNAseq",]
p <- ggplot(volcano.df.plot, 
            aes(x=Log2FC, y=-log10(FDR), colour=sig)) +
    geom_point(size=1, alpha=1, na.rm=T) +
    geom_text_repel(data=subset(volcano.df.plot, 
                                sqrt(Log2FC^2 + (-log10(FDR))^2) > quantile(sqrt(Log2FC^2 + (-log10(FDR))^2), 0.9) &
                                    FDR<0.05),
                    aes(x=Log2FC, y=-log10(FDR), label=ID), box.padding = 0.5, max.overlaps = 20, size=6) +
    scale_color_manual(name="", values=levels(volcano.df.plot$sig), guide='none') +
    geom_hline(yintercept = -log10(FDR), colour="#990000", linetype="dashed") +
    geom_vline(xintercept =  LOGFC, colour="#990000", linetype="dashed") +
    geom_vline(xintercept = -LOGFC, colour="#990000", linetype="dashed") +
    facet_wrap(~Cond, ncol=2) +
    xlab("Log2 FC RNA") + ylab("-log10 FDR") +
    xlim(-max(abs(volcano.df.plot$Log2FC)), max(abs(volcano.df.plot$Log2FC))) +
    theme_bw() +
    theme(title           = element_text(size=12, face='bold'),
          legend.text     = element_text(size=16),
          legend.position = "top",
          strip.text      = element_text(size=16),
          axis.title      = element_text(size=16),
          axis.text       = element_text(size=16)
    )
ggsave(p, file='Fig2E.pdf', width=12, height=12)

# Fig2F volcano plot for proteomics
volcano.df.plot <- volcano.df[volcano.df$IR=="IR" & volcano.df$Omics=="Prot",]
p <- ggplot(volcano.df.plot, 
            aes(x=Log2FC, y=-log10(FDR), colour=sig)) +
    geom_point(size=1, alpha=1, na.rm=T) +
    geom_text_repel(data=subset(volcano.df.plot, 
                                sqrt(Log2FC^2 + (-log10(FDR))^2) > quantile(sqrt(Log2FC^2 + (-log10(FDR))^2), 0.9) &
                                    FDR<0.05),
                    aes(x=Log2FC, y=-log10(FDR), label=ID), box.padding = 0.5, max.overlaps = 20, size=6) +
    scale_color_manual(name="", values=levels(volcano.df.plot$sig), guide='none') +
    geom_hline(yintercept = -log10(FDR), colour="#990000", linetype="dashed") +
    geom_vline(xintercept =  LOGFC, colour="#990000", linetype="dashed") +
    geom_vline(xintercept = -LOGFC, colour="#990000", linetype="dashed") +
    facet_wrap(~Cond, ncol=2) +
    xlab("Log2 FC Protein") + ylab("-log10 FDR") +
    xlim(-max(abs(volcano.df.plot$Log2FC)), max(abs(volcano.df.plot$Log2FC))) +
    theme_bw() +
    theme(title           = element_text(size=12, face='bold'),
          legend.text     = element_text(size=16),
          legend.position = "top",
          strip.text      = element_text(size=16),
          axis.title      = element_text(size=16),
          axis.text       = element_text(size=16)
    )
ggsave(p, file='Fig2F.pdf', width=12, height=12)
##-----

##- Venn diagrams ----
venn.plots <- lapply(levels(volcano.df$Omics), function(const) {
    venn.df <- volcano.df[volcano.df$Omics==const & volcano.df$IR=="IR",]
    venn.plot <- lapply(c("Log2FC < 0", "Log2FC > 0"), function(reg) {
        venn.df <- venn.df[venn.df$sig %in% colors[reg] & venn.df$FDR < 0.05,]
        venn.list <- setNames(lapply(paste0(sort(CELLS.display), "_IR"), function(cel) {
            na.omit(venn.df[venn.df$Cond==cel, 'ID', drop=T])
        }), paste0(sort(CELLS.display), "_IR"))
        p <- ggVennDiagram(venn.list, label_alpha = 0, label_size=6, set_size=6) + 
            ggtitle(paste0(ifelse(const=="Prot", "Proteomics", "RNA-seq"), ': FDR < 0.05, ', reg)) + 
            scale_fill_gradient(low='white', high='red') +
            scale_x_continuous(expand = expansion(mult = .1)) +
            theme(title           = element_text(size=16, face='bold'),
                  plot.margin     = margin(0,0,0.5,0, "cm"),
                  legend.text     = element_text(size=11),
                  legend.position = "bottom")
        return (p)
    })
    return (venn.plot)
})
ag <- ggarrange(plotlist=unlist(venn.plots, recursive = F),
                # labels=c("c", rep("", 3), "d", rep("", 3)),
                # font.label=list(size=40, face='bold'),
                # label.y=1.1,
                common.legend=T,
                ncol=2, nrow=2)
ggsave(ag, filename=paste0('Venn_log2FC0.pdf'), width=12, height=10)

venn.plots <- lapply(levels(volcano.df$Omics), function(const) {
    venn.df <- volcano.df[volcano.df$Omics==const & volcano.df$IR=="IR",]
    venn.plot <- lapply(c("Log2FC < -1", "Log2FC > 1"), function(reg) {
        venn.df <- venn.df[venn.df$sig %in% colors[reg] & venn.df$FDR < 0.05,]
        venn.list <- setNames(lapply(paste0(sort(CELLS.display), "_IR"), function(cel) {
            na.omit(venn.df[venn.df$Cond==cel, 'ID', drop=T])
        }), paste0(sort(CELLS.display), "_IR"))
        p <- ggVennDiagram(venn.list, label_alpha = 0, label_size=6, set_size=6) + 
            ggtitle(paste0(ifelse(const=="Prot", "Proteomics", "RNA-seq"), ': FDR < 0.05, ', reg)) + 
            scale_fill_gradient(low='white', high='red') +
            scale_x_continuous(expand = expansion(mult = .1)) +
            theme(title           = element_text(size=16, face='bold'),
                  plot.margin     = margin(0,0,0.5,0, "cm"),
                  legend.text     = element_text(size=11),
                  legend.position = "bottom")
        return (p)
    })
    return (venn.plot)
})
ag <- ggarrange(plotlist=unlist(venn.plots, recursive = F),
                # labels=c("c", rep("", 3), "d", rep("", 3)),
                # font.label=list(size=40, face='bold'),
                # label.y=1.1,
                common.legend=T,
                ncol=2, nrow=2)
ggsave(ag, filename=paste0('Venn_log2FC1.pdf'), width=12, height=10)

##-----

##- pathway enrichment settings ----
dirs.pattern <- "" #"IR"#, 
gmt.types <- switch(PROJECT,
                    Bristol = c(names(GMTFILES), "STRING"),
                    Lund = c(names(GMTFILES), "STRING"),
                    Endothelial = c(names(GMTFILES), "STRING"),
                    Helsinki = c(names(GMTFILES.miRNA))
)
##-----

##- combination of individual pathway enrichment on all cell types  ----
dirs <- grep(paste0(dirs.pattern, "$"), dirs, value=T)
dir.create("pathway_individual", showWarnings = T)
# signal2noise.list <- lapply(dirs, function(x) {
#     ranking <- read.table(paste0(x, "/signal2noise.tsv"), header=T, sep='\t', stringsAsFactors=F)
#     return (ranking)
# })
#signal2noise.df <- Reduce(function(x, y) merge(x, y, by = "X", all = T), signal2noise.list)
invisible(sapply(gmt.types, function(epath) {
	gostat <- setNames(lapply(1:length(dirs), function(i) {
		tab <- read.table(paste0(dirs[i], "/", epath, "_stat.tsv"), header = T, sep = "\t", quote = "")
		tab <- tab[, !grepl("setSize$|enrichmentScore$|p.adjust$", colnames(tab))]
		if (PROJECT %in% c("Lund", "Endothelial")) {
		    colnames(tab) <- sub(colnames(tab), pattern="_in_.*\\.", replacement=".")
		} else {
		    colnames(tab) <- sub(colnames(tab), pattern="^.*\\.", replacement="")
		}
		i.cond <- sub(dirs[i], pattern="IR$", replacement = "_IR")
		if (PROJECT=="Bristol") {
		    i.cond <- paste0(i.cond, ifelse(grepl("^Bristol", i.cond), "_RNAseq", "_Prot"))
		    i.cond <- sub(i.cond, pattern="^Bristol|^Prot", replacement="")
		    colnames(tab) <- paste0(i.cond, "_", colnames(tab))
		} else if (PROJECT=="Endothelial") {
		    colnames(tab) <- paste0(sub(i.cond, pattern="^RNAseq|^miRNAseq|^Prot", replacement=""),
		                            "_",
		                            sub(colnames(tab), 
		                                pattern="\\.", 
		                                replacement=ifelse(grepl("^RNAseq", i.cond), "_RNAseq.", 
		                                                   ifelse(grepl("^miRNAseq", i.cond), "_miRNAseq.", "_Prot."))))
		} else {
		    i.cond <- paste0(i.cond, ifelse(grepl("^RNAseq", i.cond), "_RNAseq", ifelse(grepl("^miRNAseq", i.cond), "_miRNAseq", "_Prot")))
		    i.cond <- sub(i.cond, pattern="^RNAseq|^miRNAseq|^Prot", replacement="")
		    colnames(tab) <- paste0(i.cond, "_", colnames(tab))
		}
	
		colnames(tab)[1] <- "ID"
		# if (grepl("^ProtUniprot|^ExosomeUniprot", dirs[i])) {
		# 	tab[, grep("core_enrichment", colnames(tab))] <- sapply(strsplit(tab[, grep("core_enrichment", colnames(tab))], split=','), 
		# function(x) {
		# 		paste(gene.df[x,'symbol'], collapse=',')
		# 	})
		# }
		return(tab)
	}), dirs)
	gostat.tab <- Reduce(function(x, y) merge(x, y, by = "ID", all = T), gostat)
	gostat.tab[, 2] <- apply(gostat.tab[, grep("Ontology", colnames(gostat.tab))], 1, function(x) unique(na.omit(x)))
	colnames(gostat.tab)[2] <- "Ontology"
	gostat.tab[, 3] <- apply(gostat.tab[, grep("Description", colnames(gostat.tab))], 1, function(x) unique(na.omit(x)))
	colnames(gostat.tab)[3] <- "Description"
	gostat.tab <- gostat.tab[, -c(grep("Ontology.|.Ontology|Description.|.Description", colnames(gostat.tab)))]
	gostat.tab <- gostat.tab[, c(1:3, 
								 grep("NES$", colnames(gostat.tab)), 
								 grep("pvalue", colnames(gostat.tab)), 
								 grep("qvalue", colnames(gostat.tab)), 
								 grep("core_enrichment", colnames(gostat.tab)))]
	write.table(gostat.tab, file = paste0("pathway_individual/", epath, "_", PROJECT, dirs.pattern, "_stat.tsv"), row.names = F, sep = "\t", quote = F)
	if (PROJECT=="All") {
	    gostat.tab.kidney <- gostat.tab[, c(1:3, 
	                                        grep("ABN|GEnC|K29|PT34", colnames(gostat.tab)))]
	    write.table(gostat.tab.kidney, file = paste0("pathway_individual/", epath, "_Kidney", "_stat.tsv"), row.names = F, sep = "\t", quote = F)
	    gostat.tab.endothelial <- gostat.tab[, c(1:3, 
	                                            grep("GEnC|ECFC|HMEC|RAEC", colnames(gostat.tab)))]
	    write.table(gostat.tab.endothelial, file = paste0("pathway_individual/", epath, "_Endothelial", "_stat.tsv"), row.names = F, sep = "\t", quote = F)
	}
	# filt.idx <- do.call(cbind, lapply(CELLS, function(cell) {
	# 	cellp <- grep(cell, grep("pvalue", colnames(gostat.tab), value = T), value = T)
	# 	apply(gostat.tab[, cellp, drop=F], 1, function(x) length(na.omit(x)) > 0 && min(na.omit(x)) < 0.05)
	# }))
	# gostat.tab.filter <- gostat.tab[rowSums(filt.idx) > 0, ]
	# write.table(gostat.tab.filter, file = paste0(epath, PROJECT, "_statfilter.tsv"), row.names = F, sep = "\t", quote = F)
}))
##-----

#### pathway enrichment with loadings from OPLS

##- STRING db ----
## string_id to gene name mapping
symbol.string <- read.table(paste0(DIR, "/../data/download/", SYMBOLSTRING), header=F, stringsAsFactors=F, row.names=3)[, -1, drop=F]
colnames(symbol.string) <- c('symbol')
## string clusters
cluster2prot <- read.delim(paste0(DIR, "/../data/download/", CLUSTERPROT), header=T, sep='\t', stringsAsFactors=F)[, c("cluster_id", "protein_id")]
cluster2gene <- data.frame(string_id=cluster2prot$cluster_id, symbol=symbol.string[cluster2prot$protein_id,], stringsAsFactors=F)
cluster2name <- read.delim(paste0(DIR, "/../data/download/", CLUSTERANNO), header=T, sep='\t', stringsAsFactors=F)[, c("cluster_id", "best_described_by")]
## string cluster relationship
cluster.relation <- read.delim(paste0(DIR, "/../data/download/", CLUSTERRELT), header=T, sep='\t', stringsAsFactors=F, row.names=2)[, -1, drop=F]
##-----

##- load files from multiblock analysis ----
if (PROJECT %in% c("Endothelial")) {
    loadingsfiles <- grep("_loadings.csv", 
                          list.files(file.path(DIR, "../multiblock", PROJECT, dirs.pattern, "results"), full.names = T, recursive = T),
                          value=T)
} else {
    loadingsfiles <- list.files(paste0(DIR, "/../multiblock/", PROJECT, "/", dirs.pattern), full.names = T)
}
loadings <- lapply(loadingsfiles, function(x) {
    tab <- switch(PROJECT,
                  Bristol  = read.csv(x, header=T, col.names=c(ifelse(grepl("_genes", x), "ensembl", "uniprot"), "predictive", "VIP")),
                  Helsinki = read.csv(x, header=T, col.names=c("predictive", "VIP", ifelse(grepl("_RNA_", x), "RNA", "miRNA"))),
                  Lund     = read.csv(x, header=T, col.names=c("predictive", "VIP", ifelse(grepl("_RNA_", x), "ensembl", "uniprot"))),
                  Endothelial = read.csv(x, header=T, col.names=if (grepl("_RNA_", x)) c("predictive", "VIP", "ensembl") else
                      if (grepl("_Prot_", x)) c("predictive", "VIP", "uniprot") else
                          if (grepl("_genes", x)) c("ensembl", "predictive", "VIP") else
                              if (grepl("_proteins", x)) c("uniprot", "predictive", "VIP")
                  )
    )
    if (grepl("_genes|_proteins", x)) tab <- tab[, c(2,3,1)]
    ## NB. This is because the order of conditions in filename of OPLS is not canonical!!! TODO: Florence, please fix!
    ## Here, we consider DiabeticSoup_versus_MannitolCtrl
    # if (PROJECT=="Endothelial" && grepl("MannitolCtrl_vs_InsulinResistant", x)) {
    #     tab$predictive <- -1*tab$predictive
    # }
    return (tab)
})
if (PROJECT %in% c("Endothelial")) {
    names(loadings) <- gsub(sapply(strsplit(loadingsfiles, split='/'), function(x) paste(x[length(x)-c(2,0)], collapse="_")), pattern=".csv", replacement="")
} else {
    names(loadings) <- gsub(sapply(strsplit(loadingsfiles, split='/'), function(x) x[length(x)]), pattern=".csv", replacement="")
}

##- check number of studied genes and proteins ----
checkrow <- sapply(dirs, function(x) {
    tab <- read.table(paste0(x, '/normalized.tsv'), header=T, sep='\t', quote="", row.names=1)
    return (rownames(tab))
})
length(Reduce(intersect, checkrow[grep("^Bristol|^RNAseq", names(checkrow))]))
length(Reduce(intersect, checkrow[grep("^Prot", names(checkrow))]))
sapply(loadings, nrow)
lengths(checkrow)
##-----

##- pathway enrichment with loadings on each cell type ----
#Ntermshowed <- 20
enrpathways <- mclapply(1:length(loadings), mc.cores=1, function(i) {
    loading <- loadings[[i]]
    print(names(loadings)[i])
    loading.idtype <- if (grepl("_genes|_RNA", names(loadings)[i])) "ensembl" else 
        if (grepl("_proteins|_Prot", names(loadings)[i])) "uniprot" else 
            if (grepl("_miRNA", names(loadings)[i])) "miRNA"
    
    ##- mapping loading.idtype and gene symbol ----
    if (loading.idtype != "miRNA") {
        gene.df <- idmap[, c(loading.idtype, "symbol")]
        ## remove any duplicates in first column: keep first occurrence
        gene.df <- gene.df[!is.na(gene.df[,1]) & !duplicated(gene.df[,1]), ]
        rownames(gene.df) <- gene.df[,1]
        gene.df <- gene.df[, -1, drop=F]
    }
    ##-----
    
    ##- read gmt files ----
    if (loading.idtype!="miRNA") {
        gmt.dfs <- mclapply(GMTFILES, mc.cores=length(GMTFILES), function(gmtf) {
            gmt.df      <- read.gmt(paste0(DIR, "/../data/download/", gmtf))
            gmt.df$term <- as.vector(gmt.df$term)
            gmt.df$gene <- as.vector(gmt.df$gene)
            return (gmt.df)
        })
    } else {
        gmt.dfs <- mclapply(GMTFILES.miRNA, mc.cores=length(GMTFILES.miRNA), function(gmtf) {
            gmt.df      <- read.gmt(paste0(DIR, "/../data/download/", gmtf))
            gmt.df$term <- as.vector(gmt.df$term)
            gmt.df$gene <- as.vector(gmt.df$gene)
            gmt.df$term  <- gsub(gmt.df$term, pattern="\\+", replacement="PLUS")
            gmt.df$term  <- gsub(gmt.df$term, pattern="'",   replacement="")
            gmt.df$term  <- gsub(gmt.df$term, pattern="-| |\\.|\\(|\\)|/|,|:|>|=|#|\\[|\\]", replacement="_")
            gmt.df$term  <- gsub(gmt.df$term, pattern="__", replacement="_")
            gmt.df$term  <- paste0("GO_", toupper(gmt.df$term))
            return (gmt.df)
        })
    }
    ##-----
    
    log2FC.conts <- ""
    ## NB. GSEA function cannot be run in parallel with seed = TRUE
    enrichedGMT <- lapply(log2FC.conts, function(cnt) {
        print(cnt)
        if (loading.idtype=="miRNA") {
            cnt.signal2noise <- setNames(loading[, paste0("predictive", cnt), drop=T], loading[, 3]) # should be the same as the following if OPLS is synchronized!
        } else {
            cnt.signal2noise <- setNames(loading[, paste0("predictive", cnt), drop=T], gene.df[loading[, 3],]) # better with predictive*sign(vip)?
        }
        
        ##- max signal2noise is used for repeated gene/uniprot names ----
        cnt.signal2noise <- cnt.signal2noise[!is.na(names(cnt.signal2noise)) & names(cnt.signal2noise)!="" & !is.na(cnt.signal2noise)]
        dup.idx <- which(duplicated(names(cnt.signal2noise)))
        cnt.signal2noise.nodup <- cnt.signal2noise
        if (length(dup.idx) > 0)
            cnt.signal2noise.nodup <- cnt.signal2noise.nodup[-dup.idx]
        for (x in unique(names(cnt.signal2noise)[dup.idx]))
            cnt.signal2noise.nodup[x] <- (cnt.signal2noise[names(cnt.signal2noise)==x])[which.max(abs(cnt.signal2noise[names(cnt.signal2noise)==x]))]
        cnt.signal2noise.nodup <- sort(cnt.signal2noise.nodup, decreasing=T)
        
        ## sort by signal2noise ratio
        geneList <- cnt.signal2noise.nodup
        
        ##- GMT enrichment ----
        egmts <- lapply(gmt.dfs, function(gmt.df) {
            term2name <- replicate(unique(gmt.df$term), n=2)
            colnames(term2name) <- c('term', 'name') #c("ont", "gene")
            egmt <- GSEA(geneList      = geneList,
                         TERM2GENE     = gmt.df,
                         TERM2NAME     = term2name,
                         #nPerm         = 1000,
                         seed          = TRUE,
                         eps           = 0,
                         minGSSize     = 5,
                         maxGSSize     = 500,
                         pAdjustMethod = ADJUST.METHOD,
                         pvalueCutoff  = 2,
                         by            = "fgsea",
                         verbose       = FALSE)

            return (egmt)
        })
        names(egmts) <- names(gmt.dfs)
        ##-----
        
        ##- STRINGdb enrichment ----
        if (PROJECT!="Helsinki") {
            estring <- GSEA(geneList      = geneList,
                            TERM2GENE     = cluster2gene,
                            TERM2NAME     = cluster2name,
                            #nPerm         = 1000,
                            seed          = TRUE,
                            eps           = 0,
                            minGSSize     = 5,
                            maxGSSize     = 500,
                            pAdjustMethod = ADJUST.METHOD,
                            pvalueCutoff  = 2,
                            by            = "fgsea",
                            verbose       = FALSE)
            egmts <- c(egmts, STRING=estring)
        }
        ##-----

        return (list(egmts=egmts, signal2noise=cnt.signal2noise))
    })
    names(enrichedGMT) <- log2FC.conts
    return (enrichedGMT)
})
names(enrpathways) <- names(loadings)
dir.create("pathway_multiblock", showWarnings = T)
invisible(mclapply(1:length(enrpathways), mc.cores=length(enrpathways), function(i) {
    enrichedGMT <- enrpathways[[i]]
    sapply(names(enrichedGMT[[1]][[1]]), function(epath) {
        enrichedGMTresult <- mclapply(enrichedGMT, mc.cores = length(enrichedGMT), function(x) {
            attr(x$egmts[[epath]], "result")
        })
        resultgmts <- Reduce(union, mclapply(enrichedGMTresult, mc.cores = length(enrichedGMTresult), function(x) {
            rownames(x)
        }))
        enrichedGMTresult.short <- do.call(cbind, mclapply(enrichedGMTresult, mc.cores = length(enrichedGMTresult), function(x) {
            if (epath=='ekegg') {
                # x[, "core_enrichment"] <- sapply(strsplit(x[, "core_enrichment"], split = "/"), function(y) {
                #     paste(rownames(gene.df[gene.df[, 1] %in% y, , drop = F]), collapse = ",")
                # })
            } else {
                x[, "core_enrichment"] <- gsub(x[, "core_enrichment"], pattern = "/", replacement = ",")
            }
            return (data.frame(x[match(resultgmts, rownames(x)), c("setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "core_enrichment")],
                       row.names=resultgmts,
                       stringsAsFactors=F))
        }))
        enrichedGMTresult.order <- enrichedGMTresult.short[, unlist(lapply(c("setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "core_enrichment"), function(x) {
            grep(x, colnames(enrichedGMTresult.short), value = T)
        }))]
        enrichedGMTresult.order <- cbind(Ontology=epath,
                                         Description=apply(do.call(cbind, 
                                                                   lapply(enrichedGMTresult, function(x) x[match(resultgmts, rownames(x)), "Description"])),
                                                           1, function(xx) unique(na.omit(xx))),
                                         enrichedGMTresult.order)
        write.table(enrichedGMTresult.order, file = paste0("pathway_multiblock/", names(enrpathways)[i], "_", epath, "_stat.tsv"), sep = "\t", quote = F, col.names = NA)
    })
}))
##-----

##- combination of multiblock pathway enrichment on all cell types ----
invisible(sapply(gmt.types, function(epath) {
    pwfiles <- list.files("pathway_multiblock", pattern = epath, full.names = T)
    pwfiles <- pwfiles[!grepl("all_", pwfiles)]
    gostat <- setNames(lapply(pwfiles, function(x) {
        tab <- read.table(x, header = T, sep = "\t", quote = "")
        tab <- tab[, !grepl("setSize$|enrichmentScore$|p.adjust$", colnames(tab))]
        colnames(tab) <- paste0(gsub(strsplit(x, split='/')[[1]][2], pattern=paste0("_", epath, "_stat.tsv"), replacement=""),
                                ".", colnames(tab))
        colnames(tab)[1] <- "ID"
        
        return(tab)
    }), 
    gsub(lapply(strsplit(pwfiles, split='/'), function(x) x[length(x)]), pattern=paste0("_", epath, "_stat.tsv"), replacement=""))
    
    gostat.tab <- Reduce(function(x, y) merge(x, y, by = "ID", all = T), gostat)
    gostat.tab[, 2] <- apply(gostat.tab[, grep("Ontology", colnames(gostat.tab))], 1, function(x) unique(na.omit(x)))
    colnames(gostat.tab)[2] <- "Ontology"
    gostat.tab[, 3] <- apply(gostat.tab[, grep("Description", colnames(gostat.tab))], 1, function(x) unique(na.omit(x)))
    colnames(gostat.tab)[3] <- "Description"
    gostat.tab <- gostat.tab[, -c(grep("Ontology.|.Ontology|Description.|.Description", colnames(gostat.tab)))]
    gostat.tab <- gostat.tab[, c(1:3, 
                                 grep("NES$", colnames(gostat.tab)), 
                                 grep("pvalue", colnames(gostat.tab)), 
                                 grep("qvalue", colnames(gostat.tab)), 
                                 grep("core_enrichment", colnames(gostat.tab)))]
    write.table(gostat.tab, file = paste0("pathway_multiblock/", epath, "_", PROJECT, 
                                          dirs.pattern, "_stat.tsv"), row.names = F, sep = "\t", quote = F)
    if (PROJECT=="All") {
        gostat.tab.kidney <- gostat.tab[, c(1:3, 
                                            grep("ABN|GEnC|K29|PT34", colnames(gostat.tab)))]
        write.table(gostat.tab.kidney, file = paste0("pathway_multiblock/", epath, "_Kidney", "_stat.tsv"), row.names = F, sep = "\t", quote = F)
        gostat.tab.endothelial <- gostat.tab[, c(1:3, 
                                                grep("GEnC|ECFC|HMEC|RAEC", colnames(gostat.tab)))]
        write.table(gostat.tab.endothelial, file = paste0("pathway_multiblock/", epath, "_Endothelial", "_stat.tsv"), row.names = F, sep = "\t", quote = F)
    }
}))
##-----


##- combine pathways from loadings with those from individual analysis ----
dir.create("pathway_combined", showWarnings = T)
invisible(sapply(gmt.types, function(epath) {
    pwfiles <- c(list.files("pathway_individual", pattern=paste0("^", epath), full.names = T),
                 list.files("pathway_multiblock", pattern=paste0("^", epath), full.names = T))
    gostat <- lapply(pwfiles, function(x) {
        tab <- read.table(x, header = T, sep = "\t", quote = "")
    })
    gostat.tab <- Reduce(function(x, y) merge(x, y, by = "ID", all = T), gostat)
    gostat.tab[, 2] <- apply(gostat.tab[, grep("Ontology", colnames(gostat.tab))], 1, function(x) unique(na.omit(x)))
    colnames(gostat.tab)[2] <- "Ontology"
    gostat.tab[, 3] <- apply(gostat.tab[, grep("Description", colnames(gostat.tab))], 1, function(x) unique(na.omit(x)))
    colnames(gostat.tab)[3] <- "Description"
    gostat.tab <- gostat.tab[, !grepl("Ontology.|.Ontology|Description.|.Description", colnames(gostat.tab))]
    # if (PROJECT=="Endothelial")
    #     colnames(gostat.tab) <- sub(colnames(gostat.tab), pattern='GEnC_IR_ConsensusOPLS', replacement="GEnC_DiabeticSoup_versus_Basal")
    if (PROJECT=="Bristol") {
        gostat.tab <- gostat.tab[, c(1:3,
                                     grep("RNAseq", colnames(gostat.tab)),
                                     grep("Prot", colnames(gostat.tab)),
                                     grep("loadings_genes", colnames(gostat.tab)),
                                     grep("loadings_proteins", colnames(gostat.tab)))]
    } else if (PROJECT=="Lund") {
        gostat.tab <- gostat.tab[, c(1:3,
                                     grep("ECFC_RNAseq", colnames(gostat.tab)),
                                     grep("RNA_loadings", colnames(gostat.tab)),
                                     grep("ECFC_Prot", colnames(gostat.tab)),
                                     grep("Prot_loadings", colnames(gostat.tab)))]
    } else if (PROJECT=="Endothelial") {
        gostat.tab <- gostat.tab[, c(1:3,
                                     grep("RNAseq\\.", colnames(gostat.tab)),
                                     grep("RNA_loadings", colnames(gostat.tab)),
                                     grep("Prot\\.", colnames(gostat.tab)),
                                     grep("Prot_loadings", colnames(gostat.tab)))]
    } else if (PROJECT=="Helsinki") {
        gostat.tab <- gostat.tab[, c(1:3,
                                     grep("_RNAseq", colnames(gostat.tab)),
                                     grep("_miRNAseq", colnames(gostat.tab)),
                                     grep("_Prot", colnames(gostat.tab)))]
    }
    if (PROJECT=="Lund") {
        gostat.tab <- gostat.tab[, c(1:3,
                                     grep("InsulinResistant_versus_Basal", colnames(gostat.tab)),
                                     grep("MannitolCtrl_versus_Basal", colnames(gostat.tab)),
                                     grep("MannitolCtrl_versus_InsulinResistant", colnames(gostat.tab)))]
    } else if (PROJECT=="Endothelial") {
        gostat.tab <- gostat.tab[, c(1:3, 
                                     grep("GEnC", colnames(gostat.tab)), 
                                     grep("ECFC", colnames(gostat.tab)), 
                                     grep("HMEC", colnames(gostat.tab)),
                                     grep("RAEC", colnames(gostat.tab)))]
        gostat.tab <- gostat.tab[, c(1:3,
                                     grep("DiabeticSoup_versus_Basal", colnames(gostat.tab)),
                                     grep("MannitolCtrl_versus_Basal", colnames(gostat.tab)),
                                     grep("DiabeticSoup_versus_MannitolCtrl", colnames(gostat.tab)),
                                     grep("DiabeticPatient_versus_Basal", colnames(gostat.tab)),
                                     grep("DiabeticPatient_versus_DiabeticSoup", colnames(gostat.tab)),
                                     grep("DiabeticPatient_versus_MannitolCtrl", colnames(gostat.tab)))]
    } else {
        gostat.tab <- gostat.tab[, c(1:3, 
                                     grep("ABN", colnames(gostat.tab)), 
                                     grep("GEnC", colnames(gostat.tab)), 
                                     grep("K29", colnames(gostat.tab)), 
                                     grep("PT34", colnames(gostat.tab)))]
    }
    gostat.tab <- gostat.tab[, c(1:3, 
                                 grep("NES$", colnames(gostat.tab)), 
                                 grep("pvalue", colnames(gostat.tab)), 
                                 grep("qvalue", colnames(gostat.tab)), 
                                 grep("core_enrichment", colnames(gostat.tab)))]
    # gostat.tab[, grep("DiabeticSoup_versus_MannitolCtrl.+NES", colnames(gostat.tab), value=T)] <- -1 * gostat.tab[, grep("DiabeticSoup_versus_MannitolCtrl.+NES", colnames(gostat.tab), value=T)]
    # NB. loadings in this contrast does not follow the case_versus_control order
    gostat.tab[, grep("DiabeticSoup_versus_MannitolCtrl.+loadings.NES", colnames(gostat.tab), value=T)] <- -1 * gostat.tab[, grep("DiabeticSoup_versus_MannitolCtrl.+loadings.NES", colnames(gostat.tab), value=T)]
    # colnames(gostat.tab) <- gsub(colnames(gostat.tab), 
    #                              pattern="DiabeticSoup_versus_MannitolCtrl", 
    #                              replacement="MannitolCtrl_versus_DiabeticSoup")
    write.table(gostat.tab, 
                file = paste0("pathway_combined/", epath, "_", PROJECT, dirs.pattern, "_withOPLS_stat.tsv"), 
                row.names = F, sep = "\t", quote = F)
    # filt.idx <- do.call(cbind, lapply(c("ABN","GEnC","K29","PT34"), function(cell) {
    # 	cellp <- grep(cell, grep("qvalue", colnames(gostat.tab), value = T), value = T)
    # 	apply(gostat.tab[, cellp, drop=F], 1, function(x) length(na.omit(x)) > 0 && min(na.omit(x)) < 0.05)
    # }))
    # gostat.tab.filter <- gostat.tab[rowSums(filt.idx) > 0, ]
    # write.table(gostat.tab.filter, file = paste0(epath, "_", PROJECT, dirs.pattern, "_withOPLS_statfilter.tsv"), row.names = F, sep = "\t", quote = F)
}))
##-----
#save(enrpathways, file="enrpathways.RData")


##- IR Bristol ----
irgene <- "INSR"
irdirs <- c("BristolallIR", "ProtallIR")
irdirs <- dirs

## INSR related genes
irstring <- rownames(symbol.string[symbol.string$symbol %in% irgene,,drop=F]) # i.e. c("9606.ENSP00000303830")
ircluster <- cluster2prot[cluster2prot$protein_id %in% irstring, 1]
ircluster.plot <- cluster.relation[setdiff(ircluster, cluster.relation[ircluster,]),]
ircluster.proteins <- cluster2prot[cluster2prot$cluster_id %in% ircluster.plot, 2]
ircluster.genes <- symbol.string[ircluster.proteins,]
paste(ircluster.proteins, collapse="','")


ir.dfs <- as.data.frame(t(sapply(irdirs, function(x) {
    DE.df <- read.table(paste0(x, '/1_DE.tsv'), header=T, sep='\t', quote="")
    ir.df <- DE.df[DE.df$name %in% ircluster.genes, , drop=F]
    rownames(ir.df) <- ir.df$name
    ir.df <- ir.df[, -c(1:2), drop=F]
    ## heatmap
    invisible(sapply(c('row', 'column', 'none'), function(scaletype) {
        pdf(paste0(x, '_irgenes_', scaletype, '.pdf'))
        par(mfrow=c(1,2))
        heatmap.2(as.matrix(ir.df[, -c(grep("versus", colnames(ir.df)), ncol(ir.df))]),
                  scale=scaletype,
                  main="STRINGdb cluster of INSR", symbreaks=F, symkey=F, trace='none', margin=c(20,10))
        heatmap.2(as.matrix(ir.df[, grep("Log2FC", colnames(ir.df))]),
                  scale=scaletype, cexCol = 0.8,
                  main="STRINGdb cluster of INSR", symbreaks=F, symkey=F, trace='none', margin=c(20,10))
        dev.off()
    }))
    ## boxplot
    lapply(rownames(ir.df), function(y) {
        irbox.df <- data.frame(count=ir.df[y, ], 
                            exp.design[colnames(counts.norm.log), c("sample", "cond", "cell", "color", "group")], 
                            stringsAsFactors=F)
    })
})))


## regression: DE gene number ~ INSR level 
ir.df <- as.data.frame(t(sapply(dirs[!grepl("Lund", dirs)], function(x) {
    tab <- read.table(paste0(x, '/DE.tsv'), header=T, sep='\t', quote="")
    
    return (c(log2FC=tab[tab$name %in% irgene, grep("^Log2FC", colnames(tab))],
              corIB=mean(cor(tab[, grepl("Basal", colnames(tab)) & !grepl("versus", colnames(tab))], tab[, grepl("InsulinResistant", colnames(tab)) & !grepl("versus", colnames(tab))])),
              Basal=mean(2^as.numeric(tab[tab$name %in% irgene, grepl("Basal", colnames(tab)) & !grepl("versus", colnames(tab))])),
              InsulinResistant=mean(2^as.numeric(tab[tab$name %in% irgene, grepl("InsulinResistant", colnames(tab)) & !grepl("versus", colnames(tab))]))))
})))

ir.df$mean <- (ir.df$Basal + ir.df$InsulinResistant)/2
ir.df$Basal <- log2(ir.df$Basal)
ir.df$InsulinResistant <- log2(ir.df$InsulinResistant)
ir.df$mean <- log2(ir.df$mean)
ir.df$IR <- factor(ifelse(grepl("IR$", rownames(ir.df)), "IR", "noIR"))
rownames(ir.df) <- gsub(rownames(ir.df), pattern="^Bristol|^Lund", replacement="RNAseq")
ir.df$numDEs <- log2(numDEs[rownames(ir.df)] + 1)
ir.df$cell <- factor(gsub(rownames(ir.df), pattern="RNAseq|Exosome|Prot|IR$", replacement=""))
ir.df$type <- factor(ifelse(grepl("RNAseq", rownames(ir.df)),
                            "RNAseq",
                            ifelse(grepl("Exosome", rownames(ir.df)), "Exosome", "Prot")),
                     levels=c("RNAseq", "Prot", "Exosome"))

ps <- list(p1=ggplot(ir.df, aes(x=Basal, y=InsulinResistant)) + 
               xlab("log2(mean counts) INSR in Basal") + 
               ylab("log2(mean counts) INSR in InsulinResistant") + 
               geom_smooth(method='lm'),
           p2=ggplot(ir.df, aes(x=Basal, y=numDEs)) + 
               xlab("log2(mean counts) INSR in Basal") + 
               ylab("log2(1 + # DE genes)"),
           p3=ggplot(ir.df, aes(x=InsulinResistant, y=numDEs)) +
               xlab("log2(mean counts) INSR in InsulinResistant") + 
               ylab("log2(1 + # DE genes)"),
           p4=ggplot(ir.df, aes(x=mean, y=numDEs)) + 
               xlab("log2(mean counts) INSR") + 
               ylab("log2(1 + # DE genes)"),
           # p5=ggplot(ir.df, aes(x=log2FC, y=numDEs)) +
           #     xlab("log2FC INSR InsulinResistant vs Basal") +
           #     ylab("log2(1 + # DE genes)"),
           p6=ggplot(ir.df, aes(x=Basal, y=corIB)) + 
               xlab("log2(mean counts) INSR in Basal") + 
               ylab("cor(InsulinResistant, Basal)"),
           p7=ggplot(ir.df, aes(x=InsulinResistant, y=corIB)) +
               xlab("log2(mean counts) INSR in InsulinResistant") + 
               ylab("cor(InsulinResistant, Basal)"),
           p8=ggplot(ir.df, aes(x=mean, y=corIB)) + 
               xlab("log2(mean counts) INSR") + 
               ylab("cor(InsulinResistant, Basal)")
           # p9=ggplot(ir.df, aes(x=log2FC, y=corIB)) +
           #     xlab("log2FC INSR InsulinResistant vs Basal") +
           #     ylab("cor(InsulinResistant, Basal)")
)
ps <- lapply(ps, function(p) {
    p <- p + geom_point(aes(colour=cell, shape=type, alpha=IR, size=IR)) +
        scale_color_manual(values = COLORGROUPS[["BEAtDKDcell"]][levels(ir.df$cell)], guide=guide_legend(title.position='top', nrow=2)) +
        scale_size_manual(values=c(4,6)) + 
        scale_alpha_manual(values=c(0.4, 1)) +
        theme_bw() +
        theme(title           = element_text(size=12, face='bold'),
              legend.text     = element_text(size=10),
              legend.position = "none",
              axis.text       = element_text(size=12))
    return (p)
})
pl <- ggplot(ir.df, aes(x=Basal, y=InsulinResistant, colour=cell, shape=type, alpha=IR, size=IR)) +
    geom_point() + 
    scale_color_manual(values = COLORGROUPS[["BEAtDKDcell"]][levels(ir.df$cell)]) +
    scale_size_manual(values=c(4,6)) + 
    scale_alpha_manual(values=c(0.4, 1))

legend <- gtable_filter(ggplot_gtable(ggplot_build(pl + theme(legend.position = "right", legend.text=element_text(size=16)))), "guide-box")
#ps[8] <- list(legend)
#ag <- grid.arrange(grobs=ps, layout_matrix=rbind(c(1:4), c(9, 5:7)))
ag <- arrangeGrob(ps[[1]], ps[[2]], ps[[3]], ps[[4]], legend, ps[[5]], ps[[6]], ps[[7]], layout_matrix=rbind(c(1:4), c(5:8)))
ggsave(ag, filename="IR.pdf", width=22, height=11)

## correlation on replicate pairwise
ircor.df <- as.data.frame(do.call(rbind, lapply(dirs[!grepl("Lund", dirs)], function(x) {
    tab <- read.table(paste0(x, '/DE.tsv'), header=T, sep='\t', quote="")
    eg <- expand.grid(colnames(tab)[grepl("Basal", colnames(tab)) & !grepl("versus", colnames(tab))], colnames(tab)[grepl("InsulinResistant", colnames(tab)) & !grepl("versus", colnames(tab))])
    res <- apply(eg, 1, function(y) {
        c(dir=x, 
          meanIR=log2(mean(2^as.numeric(tab[tab$name %in% irgene, y]))),
          minIR=log2(min(2^as.numeric(tab[tab$name %in% irgene, y]))), 
          maxIR=log2(max(2^as.numeric(tab[tab$name %in% irgene, y]))),
          cor=cor(tab[, y[1]], tab[, y[2]]))
    })
    return (t(res))
})))
ircor.df$dir <- gsub(ircor.df$dir, pattern="^Bristol|^Lund", replacement="RNAseq")
ircor.df$IR <- factor(ifelse(grepl("IR$", ircor.df$dir), "IR", "noIR"))
ircor.df$cell <- factor(gsub(ircor.df$dir, pattern="RNAseq|Exosome|Prot|IR$", replacement=""))
ircor.df$type <- factor(ifelse(grepl("RNAseq", ircor.df$dir),
                               "RNAseq",
                               ifelse(grepl("Exosome", ircor.df$dir), "Exosome", "Prot")),
                        levels=c("RNAseq", "Prot", "Exosome"))
ircor.df[, c("meanIR", "minIR", "maxIR", "cor")] <- apply(ircor.df[, c("meanIR", "minIR", "maxIR", "cor")], c(1,2), as.numeric)
p <- ggplot(ircor.df, aes(x=meanIR, y=cor)) + 
    geom_point(aes(color=cell, shape=type, alpha=IR, size=IR)) + 
    scale_color_manual(values = COLORGROUPS[["BEAtDKDcell"]][levels(ir.df$cell)]) +
    scale_size_manual(values=c(2.1,3.2)) + 
    scale_alpha_manual(values=c(0.4, 1)) +
    xlab("log2(mean counts) INSR") +
    ylab("cor(InsulinResistant, Basal)") +
    theme_bw() +
    theme(title           = element_text(size=12, face='bold'),
          legend.text     = element_text(size=10),
          legend.position = "right",
          axis.text       = element_text(size=10))
ggsave(p, filename="IRcor.pdf", width=12, height=8)

fit.RNAseq <- lm(cor~meanIR, data=ircor.df, subset = grep("RNAseq", ircor.df$type))
fit.Prot <- lm(cor~meanIR, data=ircor.df, subset = grep("Prot", ircor.df$type))
fit.Exosome <- lm(cor~meanIR, data=ircor.df, subset = grep("Exosome", ircor.df$type))
fit.RNAseq
fit.Prot
fit.Exosome

opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(fit.RNAseq, las = 1)      # Residuals, Fitted, ...
par(opar)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(fit.Prot, las = 1)      # Residuals, Fitted, ...
par(opar)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(fit.Exosome, las = 1)      # Residuals, Fitted, ...
par(opar)

sapply(levels(ircor.df$type), function (x) {
    t.test(ircor.df[ircor.df$type==x & ircor.df$IR=="IR", "cor"], ircor.df[ircor.df$type==x & ircor.df$IR=="noIR", "cor"],
           alternative = "less")
})
##-----
