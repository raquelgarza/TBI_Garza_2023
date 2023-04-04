""" File for reproducing TBI results """

import datetime
import glob
import logging
import pandas as pd
import os
import re
# from snakemake.utils import validate
import json
import multiprocessing
import subprocess
import shlex

from multiprocessing.pool import ThreadPool

import fnmatch

def get_files(base_dir: str, file_pattern: str) -> list:
  result = []
  for root, dir_names, filenames in os.walk(base_dir):
    for file_name in fnmatch.filter(filenames, file_pattern):
      c_file = os.path.join(root, file_name)
      print("found file: {}".format(c_file))
      result.append(c_file)
  return result

def call_proc(cmd):
    """ This runs in a separate thread. """
    #subprocess.call(shlex.split(cmd))  # This will block until cmd finishes
    p = subprocess.Popen(shlex.split(cmd))
    p.communicate()


configfile: "tbi_qc.json"

DATA_DIR=config["DATA_DIR"]
PATTERN=config["PATTERN"]
OUT_DIR=config["OUT_DIR"]
HUMAN_GENE_IDS_FILE=config["HUMAN_GENE_IDS_FILE"]

# print("DATA_DIR: {} OUT_DIR: {}".format(DATA_DIR, OUT_DIR))

## hard-coded here - can be moved to config if changed
CELLRANGER_FILENAMES = {}
CELLRANGER_FILENAMES["matrix"]="matrix.mtx.gz"
CELLRANGER_FILENAMES["features"]="features.tsv.gz"
CELLRANGER_FILENAMES["barcodes"]="barcodes.tsv.gz"

## this is directly generated from reading the data_directory 
## see steps `find_samples` and `create_sampledf` which can be used to update sampleDf 
SAMPLEDF_FILE_PRE=os.path.join(OUT_DIR, "sampleDf_pre.tsv")
SAMPLEDF_FILE=os.path.join(OUT_DIR, "sampleDf.tsv")
UNPROCESSED_RDS_FILE=os.path.join(OUT_DIR, "unprocessed.rds")

DECONTAMINATED_DIR=os.path.join(OUT_DIR, "decontaminated")
DECONTAMINATED_RESULT_FILE=os.path.join(DECONTAMINATED_DIR, "decontaminated.rds")
DECONTAMINATED_SAMPLE_FILE=os.path.join(DECONTAMINATED_DIR, "decontaminated_sampleDf.tsv")

QF_DIR=os.path.join(OUT_DIR, "quality_filtering")
QF_RESULT_FILE=os.path.join(QF_DIR, "results-quality-filtered.rds")
QF_SAMPLE_FILE=os.path.join(QF_DIR, "sampleDf-quality-filtered.tsv")
QF_LOG_FILE=os.path.join(QF_DIR,  "do_quality_filtering.log")

PROCESSED_DIR=os.path.join(OUT_DIR, "processed")


MERGED_SAMPLE_DIR=os.path.join(OUT_DIR, "merged")
# MERGE_SAMPLES=[key for key in config["MERGE_SAMPLES"]]
# MERGE_STRATEGIES=[key for key in config["MERGE_STRATEGIES"]]

MERGE_SAMPLES=["ALL"]
MERGE_STRATEGIES=["SeuratV2", "Harmony"] # "SeuratV3" is crashing out

MERGED_SUBDIR = [ x + "_" + y for x in MERGE_SAMPLES for y in MERGE_STRATEGIES] 

UNPROCESSED_READ10X_RDS = os.path.join(OUT_DIR, "unprocessed_Read10X.rds")
UNPROCESSED_READ10X_SAMPLEDF = os.path.join(OUT_DIR, "sampleDf_Read10X.tsv")

rule all:
    input:
      PROCESSED_DIR

rule create_pre_sampledf:
    output:
      SAMPLEDF_FILE_PRE
    params:
      script="sc_pipeline/create_sampledf_from_cellranger_files.R",
      pattern=config["PATTERN"]
    shell:"Rscript {params.script} -i {DATA_DIR} -o {output} -p {params.pattern}; echo '\n\n'; cat {output};"

rule create_sampledf: 
   input: 
       SAMPLEDF_FILE_PRE
   output: 
       SAMPLEDF_FILE
   shell: "cp {input} {output}"

rule Read10X:
    input:
        SAMPLEDF_FILE
    output:
        rds=UNPROCESSED_READ10X_RDS,
        out_sample_df=UNPROCESSED_READ10X_SAMPLEDF
    params:
        script="sc_pipeline/read_10x_using_seurat.R"
    shell: "Rscript {params.script} {input} {output.rds} --out_sample_df {output.out_sample_df}"

rule read_cellranger_files:
    input:
        SAMPLEDF_FILE
    output:
        UNPROCESSED_RDS_FILE
    log:
        os.path.join(OUT_DIR, "read_cellranger_files.log")
    params: 
        script="sc_pipeline/read_cellranger_files.R"
    shell: "Rscript {params.script} -i {input} -f {output} -l {log}"


rule do_quality_filtering:
    input: 
      SAMPLEDF_FILE=SAMPLEDF_FILE,
      RDS_FILE=UNPROCESSED_RDS_FILE
    output:
      SAMPLEDF=QF_SAMPLE_FILE,
      RDS=QF_RESULT_FILE
    log: QF_LOG_FILE
    params:
        QF_SCRIPT="sc_pipeline/do_quality_filtering.R",
        OUT_DIR=QF_DIR,
        FILTER_MITO_NMADS=0, 
        FILTER_MIN_GENES=1000
    shell:
        "Rscript {params.QF_SCRIPT} --SAMPLEDF_FILE {input.SAMPLEDF_FILE} --RDS_FILE {input.RDS_FILE} --OUT_DIR {params.OUT_DIR} --FILTERED_SAMPLEDF_FILE {output.SAMPLEDF} --RESULTS_FILE {output.RDS} --LOG_FILE {log} --FILTER_MITO_NMADS {params.FILTER_MITO_NMADS} --FILTER_MIN_NUM_GENES_EXPRESSED {params.FILTER_MIN_GENES}"


rule merge_samples:
  input: 
    RDS=QF_RESULT_FILE,
    SAMPLEDF=QF_SAMPLE_FILE
  output:
    sample_file=expand("{out_dir}/{merge_sample}_{merge_strategy}/merge_options.json",
                           out_dir=MERGED_SAMPLE_DIR,
                           merge_sample=MERGE_SAMPLES,
                           merge_strategy=MERGE_STRATEGIES), 
    RDS=expand("{out_dir}/{merge_sample}_{merge_strategy}/results-merged.rds",
                   out_dir=MERGED_SAMPLE_DIR,
                   merge_sample=MERGE_SAMPLES,
                   merge_strategy=MERGE_STRATEGIES)
  log: 
    expand("{out_dir}/{merge_sample}_{merge_strategy}/merge_seurat_samples.log",
               out_dir=MERGED_SAMPLE_DIR,
               merge_sample=MERGE_SAMPLES,
               merge_strategy=MERGE_STRATEGIES)
  params:
    script="sc_pipeline/merge_seurat_samples.R",
  run: 
    os.makedirs(MERGED_SAMPLE_DIR, exist_ok=True)
    for merge_sample in MERGE_SAMPLES: 
      for merge_strategy in MERGE_STRATEGIES: 
        print("sample: {} strategy: {}".format(merge_sample, merge_strategy))
        curr_out_dir = MERGED_SAMPLE_DIR + "/" + merge_sample + "_" + merge_strategy + "/"
        sample_file = curr_out_dir +"merge_options.json"
        rds_file = curr_out_dir + "results-merged.rds"
        log_file = curr_out_dir + "merge_seurat_samples.log"
        
        ## write json options
        with open(sample_file, "w") as fid: 
          result = dict()
          result["sample"]=config["MERGE_SAMPLES"][merge_sample]
          result["strategy"]=config["MERGE_STRATEGIES"][merge_strategy]
          json.dump(result, fid)
        
        cmd = "Rscript " + params.script 
        cmd = cmd + " --IN_SAMPLEDF_FILE " + input.SAMPLEDF
        cmd = cmd + " --IN_RDS_FILE "  + input.RDS 
        cmd = cmd + " --OUT_RDS_FILE " + rds_file
        cmd = cmd + " --MERGE_JSON_FILE " + sample_file
        cmd = cmd + " --RESOLUTION " + str(config["SEURAT_RESOLUTION"])
        cmd = cmd + " --LOG_FILE " + log_file
        # run the script
        print(cmd)
        os.system(cmd)  

rule do_panglodb_enrichment: 
  input: get_files(PROCESSED_DIR, "processed.rds")
  output: [z.replace("processed.rds", "enrichmentResult.rds") for z in get_files(PROCESSED_DIR, "processed.rds")]
  params: 
    script="sc_pipeline/do_panglodb_enrichment.R"
  run: 
    for filename in input:
      curr_out_dir = os.path.dirname(filename)
      curr_log_file = filename.replace("processed.rds", "do_panglodb_enrichment.log")
      print("get_enrichment_results: {} -> {}".format(filename, curr_out_dir))
      cmd = "Rscript " + params.script
      cmd = cmd + " -i " + filename + " -o " + curr_out_dir 
      cmd = cmd + " -p " + config["PANGLODB_DIR"] 
      cmd = cmd + " -l " + curr_log_file
      os.system(cmd)
      
