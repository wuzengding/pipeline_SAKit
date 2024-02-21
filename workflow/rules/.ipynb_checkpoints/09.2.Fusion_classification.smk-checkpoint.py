rule Fusion_classification:
    input:
        fusion_gff = rules.Fusion_calling.output.fusion_gff,
        fusion_rep = rules.Fusion_calling.output.fusion_rename_rep

    output:
        fusion_class = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_classification.txt".format(sample=SampleID)),
        fusion_junc = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_junctions.txt".format(sample=SampleID)),
        refAnnotation_lq_isoforms = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_corrected.genePred".format(sample=SampleID)),
        fusion_gtf_corrected = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_corrected.gtf".format(sample=SampleID)),
        fusion_rep_corrected = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_corrected.fasta".format(sample=SampleID))

    params:
        genome = config["Reference"]["genome"],
        genome_annotation = config["Reference"]["annotgtf"],
        classification_dir = os.path.join(outpath,"results/06.Fusion_calling/"),
        classification_prefix = "{sample}.fusion".format(sample=SampleID)

    threads:
        1
    log:
        os.path.join(outpath,"log/{0}.Fusion_classification.log".format(SampleID))
    conda:
        "SQANTI3.env"
        
    shell:
        """
        export PYTHONPATH=/mnt/user/wzd/03.biotools/software/cDNA_Cupcake/sequence:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake
        
        {config[ScriptTools][sqanti_qc]} \
                --is_fusion \
                -d {params.classification_dir} \
                --orf_input {input.fusion_rep} \
                -o {params.classification_prefix} \
                --report pdf \
                {input.fusion_gff} \
                {params.genome_annotation} \
                {params.genome}   > {log} 2>&1
               
        unset PYTHONPATH
        """
