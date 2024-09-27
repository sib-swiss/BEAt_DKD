#!/usr/bin/sh

## ABN = Pod, GEnC = GEC, K29 = MC, PT34 = PTC
for out in Bristol
do
    for omics in RNAseq Prot
    do
        for cell in all ABN GEnC K29 PT34
        do
            for ir in "" "IR" "all" 
            do
                nohup Rscript -e 'library(rmarkdown); library(knitr); 
                PROJECT <- "BEAt_DKD";
                OUT <- "'$out'";
                OMICS <- "'$omics'";
                CELL <- "'$cell'";
                IR <- "'$ir'";
                DESIGNFILE <- paste0(OUT, "_design.cfg");
                COUNTFILE <- ifelse(OMICS=="Prot", "data/clean/proteomics/Scaled.tsv", NA);
                ADJUSTALL <- TRUE; FDR <- 0.05; LOGFC <- 1;
                GENEANNOTFILE <- NA; GENESUB <- "all";
                render("~/src/Bristol.Rmd", output_dir=paste0("output/", OUT, "/", OMICS, CELL, IR, "/"), output_file="report.pdf", clean=T);' &
            done
        done
    done
done
