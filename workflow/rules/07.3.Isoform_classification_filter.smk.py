
rule Isoform_classification_filter:
    input:
        sqanti3_classification = rules.Isoform_classification.output[0],
        sqanti3_classification_corrected_fa = rules.Isoform_classification.output[2],
        sqanti3_classification_corrected_gtf = rules.Isoform_classification.output[3],
        sqanti3_classification_corrected_faa = rules.Isoform_classification.output[4]

    
    output:
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.filtered_lite_reasons.txt".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.filtered_lite_SQANTI3_report.pdf".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.filtered_lite.gtf".format(sample=SampleID)), 
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.filtered_lite_junctions.txt".format(sample=SampleID)), 
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.filtered_lite_classification.txt".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.filtered_lite.fasta".format(sample=SampleID))
    
    params:
        pythonpath = config["Env"]["python_lib"],
        filterModel = "ML",
        intrapriming = 0.6,
        outprefix = SampleID,
        min_cov = 3,
        output_dir = os.path.join(outpath, "results/05.Isoform_Novel/"),
        input_dir = os.path.join(config["OUTPATH"], "results/IsoformClassification/")
    
    conda:
        "SQANTI3.env"

    log:
        os.path.join(outpath,"log/{0}.Isoform_classification_filter.log".format(SampleID))

    threads:
        1
    
    shell:
        """
        export PYTHONPATH=/mnt/user/wzd/03.biotools/software/cDNA_Cupcake/sequence:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake

            {config[ScriptTools][sqanti_filter]}  \
            {input.sqanti3_classification} {input.sqanti3_classification_corrected_fa}  {input.sqanti3_classification_corrected_gtf} \
            --faa {input.sqanti3_classification_corrected_faa} \
            --intrapriming {params.intrapriming} \
            --min_cov {params.min_cov} \
            --report pdf \
            --saturation \
             >{log} 2>&1
             
        unset PYTHONPATH
         """