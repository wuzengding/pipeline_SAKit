import os
import sys
print(sys.executable)
import pandas as pd
import tempfile
#import pprint
import ntpath
from Bio import SeqIO
import gzip
from collections import Counter
import itertools
import logging

if config["species"]=="mouse":
    configfile: "/mnt/user/wzd/05.pipeline_dev/pipeline_SAKit_v1.0/config/config_mouse.yml"
elif config["species"]=="human":
    configfile: "/mnt/user/wzd/05.pipeline_dev/pipeline_SAKit_v1.0/config/config.yml"

configfile: config["SampleYml"]
#print(config["species"])
#print(config["SampleYml"])
#print(config["PACBIOREADS"]["Subreadbam"])
#print(config["NGSShortReads"]["fastq_R1"])

#print(config)
SampleID = config["SampleID"]
outpath = config["OUTPATH"]
rule all:
    input:
        # Isoform  outputs
        sqanti3_classification_final_results = os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_classification_final_results.tsv".format(SampleID)),
        sqanti3_classification_filtered_out = os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_classification_filtered_out_isoforms_results.tsv".format(SampleID)),
        
        ## Fusion outputs
        fusion_class_fusionhub = os.path.join(outpath,"results/06.Fusion_calling/{0}.fusion_classification_final_results_fusionhub.tsv".format(SampleID)),
        
        ## SVSNV outputs
        pbsv_vcf_annot_tsv = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.snpeff.tsv".format(SampleID)),
        pbsv_snv_vcf_tsv = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.pacbio.single_variant.vcf".format(SampleID)),
        snpeff_res = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.single_variant.snpeff.vcf".format(SampleID)),

        ## QC
        qcsum = os.path.join(outpath,"results/00.QC/{0}.runqc.txt".format(SampleID)),
        
        # Make cmd
        cmdlog = os.path.join(outpath,"cmd/{0}.all_command.sh".format(SampleID))

include: "workflow/rules/01.ccs.smk.py"
include: "workflow/rules/02.demultiplex.smk.py"
include: "workflow/rules/03.refine.smk.py"
include: "workflow/rules/04.cluster.smk.py"
include: "workflow/rules/05.align.smk.py"
include: "workflow/rules/06.1.Isoform_calling.smk.py"
include: "workflow/rules/06.2.Isoform_filter_bycounts.smk.py"
include: "workflow/rules/06.3.Isoform_filter_subset.smk.py"
include: "workflow/rules/12.Align_NGS.smk.py"
include: "workflow/rules/07.1.Isoform_quantfs.smk.py"
include: "workflow/rules/07.2.Isoform_classification.smk.py"
include: "workflow/rules/07.3.Isoform_classification_filter.smk.py"
include: "workflow/rules/08.Isoform_summary.smk.py"
include: "workflow/rules/09.1.Fusion_calling.smk.py"
include: "workflow/rules/09.2.Fusion_classification.smk.py"
include: "workflow/rules/10.Fusion_annotation.py"
include: "workflow/rules/11.SV_calling.smk.py"
include: "workflow/rules/13.SnvIndelCalling.smk.py"
include: "workflow/rules/15.QC.py"
include: "workflow/rules/16.Make_cmd.py"
