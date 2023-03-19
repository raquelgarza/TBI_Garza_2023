# Quality Filtering for TBI Manuscript

## Organization 
- `sc_pipeline/` - R scripts used in the snakemake pipeline
- `rcode/` - Directory containing helper functions (used in `sc_pipeline/` scripts)
- `tbi_qc.smk` - Main snakemake file which reproduces the quality filtering steps. 


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
renv::init(bioconductor = T)
renv::install()
```
If this does not work, please do a fresh install using the steps below. 
#### Fresh install using `renv`
```R 
require("renv")
renv::init(bioconductor = T, bare=T)
renv::install(c("getopt", "futile.logger", "Seurat", "bioc::monocle", "bioc::scater", "bioc::clusterProfiler", "bioc::simpleSingleCell", "ellipse"))
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
