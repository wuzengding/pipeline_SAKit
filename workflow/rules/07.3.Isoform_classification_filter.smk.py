
rule Isoform_classification_filter:
    input:
        sqanti3_classification = rules.Isoform_classification.output[0],
        sqanti3_classification_corrected_fa = rules.Isoform_classification.output[2],
        sqanti3_classification_corrected_gtf = rules.Isoform_classification.output[3],
        sqanti3_classification_corrected_faa = rules.Isoform_classification.output[4]

    
    output:
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}_isoform_filtered_lite_reasons.txt".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}_isoform_filtered_lite_SQANTI3_report.pdf".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}_isoform_filtered_lite.gtf".format(sample=SampleID)), 
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}_isoform_filtered_lite_junctions.txt".format(sample=SampleID)), 
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}_isoform_filtered_lite_classification.txt".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}_isoform_filtered_lite.fasta".format(sample=SampleID))
    
    params:
        pythonpath = config["Env"]["python_lib"],
        filterModel = "ML",
        intrapriming = 0.6,
        outprefix = SampleID,
        min_cov = 3,
        output_dir = os.path.join(outpath, "results/05.Isoform_Novel/"),
        input_dir = os.path.join(config["OUTPATH"], "results/IsoformClassification/")


    log:
        os.path.join(outpath,"log/Isoform_classification_filter.log")

    threads:
        1
    
    shell:
        """
            {config[Env][python3]} {config[ScriptTools][sqanti_filter]}  \
            {params.filterModel}  \
            --isoforms {input.sqanti3_classification_corrected_fa} \
            --gtf {input.sqanti3_classification_corrected_gtf} \
            --faa {input.sqanti3_classification_corrected_faa} \
            --intrapriming {params.intrapriming} \
            --min_cov {params.min_cov} \
            --report pdf \
            --dir {params.output_dir} \
            --output {params.outprefix} \
            {input.sqanti3_classification} 2>{log}
         """