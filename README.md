# Pipelines + statistical and visualization scripts

This repository presents analyses for the manuscript
[Single-cell transcriptomics of resected human traumatic brain injury tissues reveals acute activation of endogenous retroviruses in oligodendroglia](https://www.biorxiv.org/content/10.1101/2022.09.07.506982v1).

## Organization
- `src` - Directory containing all pipelines and scripts
- `src/config_files` - Json files required to run the [**Snakefile_snRNAseq_QC**](./src/Snakefile_snRNAseq_QC) pipeline.
- [**Snakefile_bulkRNAseq**](./src/Snakefile_bulkRNAseq) - Snakemake pipeline to preprocess bulk RNAseq data.
- [**Snakefile_snRNAseq_QC**](./src/Snakefile_snRNAseq_QC) - Snakemake file which reproduces the quality filtering steps.
- [**main_combinedUMAP_perCluster.py**](./src/main_combinedUMAP_perCluster.py) - trusTEr script to quantify TEs per cluster, per condition.
- [**main_combinedUMAP_perCluster_perSample.py**](./src/main_combinedUMAP_perCluster_perSample.py) - trusTEr script to quantify TEs per cluster, per sample.
- [**Snakefile_snRNAseq_celltype**](./src/Snakefile_snRNAseq_celltype) - Snakemake pipeline to index and convert to bigWig files the cluster BAM files output from trusTEr.
- `src/sc_pipeline` - R scripts used by [**Snakefile_snRNAseq_QC**](./src/Snakefile_snRNAseq_QC)
- `src/rcode` - Directory containing helper functions (used in `sc_pipeline/` scripts) and R markdowns for downstream analyses and visualization for the figures. Main scripts to look for:
	+ [**TBI_bulk.Rmd**](./src/rcode/TBI_bulk.Rmd) - R markdown for the statistical analyses and visualization of bulk RNAseq data
		* Normalization
		* Differential expression analysis of ERVs
		* Characterization of putative proteins and LTRs
	+ [**TBI_snRNAseq_TE_expression.Rmd**](./src/rcode/TBI_snRNAseq_TE_expression.Rmd) - Visualization of trusTEr's output (TE quantification per cluster)
	+ [**TBI_snRNAseq_gene_expression.Rmd**](./src/rcode/TBI_snRNAseq_TE_expression.Rmd) - Cell characterization and differential gene expression analysis of snRNAseq
		* Cell characterization and visualization of the data
		* Cell cycle scoring
		* Differential gene expression analysis per cell type
		* Gene set enrichment analysis per cell type
		* Enrichment scores for innate immunity related genes
		* Visualization of gene expression
	+ [**TBI_snRNAseq_gene_expression_QC_preprocessing.Rmd**] - Visualization of quality control metrics for snRNAseq data
- `data/` - Directory containing input matrix counts (saved using `git lfs`)
- `output/` - Directory where output artifacts will be generated upon rerunning the pipeline
- `tables/` - Directory containing tables with differential gene expression analysis per cell type and gene set enrichment analysis results.

## Installation
#### On Mac M1 
- Install using Homebrew.
```shell
brew install r openjdk@11 fribidi python@3.9 git-lfs
pip3 install snakemake  # for running the pipeline 
git lfs install  # for getting the input matrix files

# see https://krishanv.com/posts-output/install-gmp/ and https://pat-s.me/transitioning-from-x86-to-arm64-on-macos-experiences-of-an-r-user/#r-packages---source-installations 
sudo ln -sfn /opt/homebrew/opt/openjdk@11/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk-11.jdk
```


#### Install from `renv` lock file
```R
require("renv")
renv::restore()
```
If this does not work, please do a fresh install using the steps below. 
#### Fresh install using `renv`
```R 
require("renv")
renv::init(bioconductor = T)
renv::install(c("getopt", "futile.logger", "bioc::harmony", "Seurat", "bioc::monocle", "bioc::scater", "bioc::clusterProfiler", "bioc::simpleSingleCell", "ellipse"))
renv::snapshot()
```

## Running the QC steps

To reproduce the QC analysis, please do the following steps: 

- Run `git lfs fetch --all` to download the raw data files into `data/` directory.
- Run the following command: 

```shell 
snakemake -f tbi_qc.smk merge_samples 
``` 
This will create run the QC steps in `output/` sub-directory. The final output file is `output/merged/ALL_Harmony/results-merged.rds` which is used in downstream analysis.

