### Integrative transcriptome and proteome profiling of insulin-resistant kidney cell models and human biopsies reveals common and cell-type-specific mechanisms underpinning Diabetic Kidney Disease


#### Cell-specific differential expression analysis

Execute `src/DE.sh` to produce DE tables for all cell types as well as for each individual cell type with and with and without insulin receptor over-expression. The input and output can be obtained from the NCBI BioProject PRJNA905899 or upon request.

#### Downstream analysis

The script `src/Bristol_postprocess.R` provides a post-processing and produced some figures in the article. `src/Bristol_GOcluster.R` provides the pathway enrichment analysis and clustering of the pathways. `src/ConsensusOPLS.m` and `src/CCSWA_+IR_irvoom.Rmd` provides the ConsensusOPLS analysis and post-processing.
