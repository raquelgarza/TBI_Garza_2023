#!/bin/python
import truster
import os

path_parent = os.path.dirname(os.getcwd())
raw_path = {"healthy":"CTG_Seq109_110/210202_A00681_0301_BHWFL5DMXX/HWFL5DMXX/outs/fastq_path/10Xsc/", "TBI" : ["CTG_Seq093_98_SingleNuclei/FastQFiles/", "Seq073_Glioma/Seq073_75/HNWHHDMXX/outs/fastq_path/10Xsc/", "CTG_Seq093_Seq104/MJ/"]}
lunarc = "config_files/lunarc_config.json"
modules = "config_files/software_modules.json"

tbi = truster.Experiment("tbi_per_condition", lunarc, modules)

tbi.register_sample(sample_id = "TBI_HuBrainCTL_Nuclei501F_Hg38", sample_name = "TBI_HuBrainCTL_Nuclei501F_Hg38", raw_path = os.path.join(raw_path["healthy"], "Seq109_11"))
tbi.register_sample(sample_id = "TBI_HuBrainCTL_Nuclei501T_Hg38", sample_name = "TBI_HuBrainCTL_Nuclei501T_Hg38", raw_path = os.path.join(raw_path["healthy"], "Seq109_12"))
tbi.register_sample(sample_id = "TBI_HuBrainCTL_Nuclei529F_Hg38", sample_name = "TBI_HuBrainCTL_Nuclei529F_Hg38", raw_path = os.path.join(raw_path["healthy"], "Seq109_13"))
tbi.register_sample(sample_id = "TBI_HuBrainCTL_Nuclei529T_Hg38", sample_name = "TBI_HuBrainCTL_Nuclei529T_Hg38", raw_path = os.path.join(raw_path["healthy"], "Seq109_14"))
tbi.register_sample(sample_id = "TBI_HuBrainCTL_Nuclei502T_Hg38", sample_name = "TBI_HuBrainCTL_Nuclei502T_Hg38", raw_path = os.path.join(raw_path["healthy"], "Seq109_6"))

tbi.register_sample(sample_id = "HuBrain_TBI_no6", sample_name = "HuBrain_TBI_no6", raw_path = os.path.join(raw_path["TBI"][1], "Seq073_1"))
tbi.register_sample(sample_id = "HuBrain_TBI_no7", sample_name = "HuBrain_TBI_no7", raw_path = os.path.join(raw_path["TBI"][1], "Seq073_2"))
tbi.register_sample(sample_id = "HuBrain_TBI_no8", sample_name = "HuBrain_TBI_no8", raw_path = os.path.join(raw_path["TBI"][1], "Seq073_3"))
tbi.register_sample(sample_id = "MJ_TBI_nr10", sample_name = "MJ_TBI_nr10", raw_path = os.path.join(raw_path["TBI"][2], "MJ_TBI_nr10"))
tbi.register_sample(sample_id = "MJ_TBI_nr11", sample_name = "MJ_TBI_nr11", raw_path = os.path.join(raw_path["TBI"][2], "MJ_TBI_nr11"))
tbi.register_sample(sample_id = "MJ_TBI_nr16", sample_name = "MJ_TBI_nr16", raw_path = os.path.join(raw_path["TBI"][2], "MJ_TBI_nr16"))
tbi.register_sample(sample_id = "MJ_TBI_nr19", sample_name = "MJ_TBI_nr19", raw_path = os.path.join(raw_path["TBI"][2], "MJ_TBI_nr19"))
tbi.register_sample(sample_id = "MJ_TBI_nr20", sample_name = "MJ_TBI_nr20", raw_path = os.path.join(raw_path["TBI"][2], "MJ_TBI_nr20"))
tbi.register_sample(sample_id = "MJ_TBI_nr21", sample_name = "MJ_TBI_nr21", raw_path = os.path.join(raw_path["TBI"][2], "MJ_TBI_nr21"))
tbi.register_sample(sample_id = "MJ_TBI_nr1", sample_name = "MJ_TBI_nr1", raw_path = os.path.join(raw_path["TBI"][0], "Seq093_1"))
tbi.register_sample(sample_id = "MJ_TBI_nr2", sample_name = "MJ_TBI_nr2", raw_path = os.path.join(raw_path["TBI"][0], "Seq093_2"))
tbi.register_sample(sample_id = "MJ_TBI_nr3", sample_name = "MJ_TBI_nr3", raw_path = os.path.join(raw_path["TBI"][0], "Seq093_7"))

print([sample.sample_name for sample in tbi.samples.values()])
quantification_dir = os.path.join(path_parent, "1_count_premRNA/")
index_dir = "premRNAREF_SingleCells/GRCh38_premRNA/" 


for sample_id in list(tbi.samples.keys()):
        tbi.set_quantification_outdir(sample_id = sample_id, cellranger_outdir = os.path.join(quantification_dir, sample_id))

merged_dir = os.path.join(path_parent, "3_combinedUMAP_perCluster") 
gene_gtf = "annotations/hg38/gencode/v38/gencode.v38.annotation.gtf"
te_gtf = "annotations/hg38/repeatmasker/hg38_rmsk_TEtranscripts.gtf"
star_index = "GRCh38.p13_gencode.v38_STAR/" 

tbi.set_merge_samples_outdir(merged_dir)
merged_pipeline_dir = os.path.join(merged_dir, "clusterPipeline_per_condition")

tbi.process_clusters(mode = "merged", outdir = merged_pipeline_dir, gene_gtf = gene_gtf, te_gtf = te_gtf, star_index = star_index, RAM = 48725506423, jobs=17, groups = {"control" : ["TBI_HuBrainCTL_Nuclei501F_Hg38", "TBI_HuBrainCTL_Nuclei501T_Hg38", "TBI_HuBrainCTL_Nuclei502T_Hg38", "TBI_HuBrainCTL_Nuclei529F_Hg38", "TBI_HuBrainCTL_Nuclei529T_Hg38"], "tbi" : ["HuBrain_TBI_no6", "HuBrain_TBI_no7", "HuBrain_TBI_no8", "MJ_TBI_nr1", "MJ_TBI_nr2", "MJ_TBI_nr3", "MJ_TBI_nr10", "MJ_TBI_nr11", "MJ_TBI_nr16", "MJ_TBI_nr19", "MJ_TBI_nr20", "MJ_TBI_nr21"]}, unique = True, tsv_to_bam = False, filter_UMIs = False, bam_to_fastq = False, concatenate_lanes = False, merge_clusters = False, map_cluster = True, TE_counts = False, normalize_TE_counts = False) 

