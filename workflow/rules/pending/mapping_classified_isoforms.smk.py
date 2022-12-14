rule run_classified_tx_mapping:
    input:
        sqanti3_output_filtered  = rules.Isoform_classification_filter.output[5]

    output:
        sqanti3_output_filtered_bam_tx = os.path.join(outpath,"results/FinalResults/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.bam".format(SampleID))

    params:
        genome = config["Reference"]["genome"],
        ax = "splice", 
        secondary = "no",  
        o = "-O6,24", 
        b4 = "-B4",
        uf = "-uf",
        temp_out_sam = os.path.join(outpath,"results/FinalResults/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.sam".format(SampleID)),
        picard = config["Softwaretools"]["picard"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/data/tmp",
        sort_order = "coordinate"

    log:
        os.path.join(outpath,"log/run_classified_tx_mapping.log")

    threads:
        18

    shell:
        """
            {config[Softwaretools][minimap2]} \
                -t {threads} \
                -ax {params.ax} \
                --secondary={params.secondary} \
                {params.o} \
                {params.b4} \
                {params.uf} \
                {params.genome} \
                {input.sqanti3_output_filtered} > {params.temp_out_sam} 2> {log}
            
            java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                -jar  {params.picard} \
                SortSam \
                I={params.temp_out_sam} \
                O={output.sqanti3_output_filtered_bam_tx} \
                SORT_ORDER={params.sort_order}  2>> {log}
            
            samtools index {output.sqanti3_output_filtered_bam_tx} 2>> {log}

            rm {params.temp_out_sam}
        """
