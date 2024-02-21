def filter_fusion_reads( final_pbfusion_fa, fusion_class_cor_fasta, fusion_final_results):

    fusion_final_results_pd =  pd.read_csv( fusion_final_results, sep ='\t')
    pbfusion_name_list  = fusion_final_results_pd ["UniqueID"]

    with open( final_pbfusion_fa, "w" ) as output_handle:
        for seq_record in SeqIO.parse( fusion_class_cor_fasta, "fasta"):
            for pb in pbfusion_name_list:
                if pb.lower() in seq_record.description.lower():
                    #print(seq_record.format("fasta"))
                    sequences =seq_record
                    SeqIO.write(sequences, output_handle, "fasta")
                    break


rule run_filter_pbfusion_to_bam:
    input:
        fusion_rep = rules.Fusion_calling.output.fusion_rep,
        fusion_final_results = rules.run_fusion_events.output.final_fusion_class_fusionhub

    output:
        final_pbfusion_fa = os.path.join(outpath,"results/FinalResults/final_pbfusion.fasta"),
        final_pbfusion_mapped_bam = os.path.join(outpath,"results/FinalResults/{0}_final_pbfusion_mapped_hg38_sorted.bam".format(SampleID)),
        final_pbfusion_mapped_bam_bai = os.path.join(outpath,"results/FinalResults/{0}_final_pbfusion_mapped_hg38_sorted.bam.bai".format(SampleID))

    params:
        genome = config["Reference"]["genome"],
        ax = "splice", 
        secondary = "no",  
        o = "-O6,24", 
        b4 = "-B4",
        uf = "-uf",
        temp_out_sam = os.path.join(outpath,"results/FinalResults/final_pbfusion.sam"),
        picard = config["Softwaretools"]["picard"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/data/tmp",
        sort_order = "coordinate"
    threads:
        18
    log:
        os.path.join(outpath,"log/run_filter_pbfusion_to_bam.log")
    run:
        filter_fusion_reads( \
                                output.final_pbfusion_fa, 
                                input.fusion_rep, 
                                input.fusion_final_results \
                            )

        shell(\
                """
                    {config[Softwaretools][minimap2]} \
                        -t {threads} \
                        -ax {params.ax} \
                        --secondary={params.secondary} \
                        {params.o} \
                        {params.b4} \
                        {params.uf} \
                        {params.genome} \
                        {output.final_pbfusion_fa} > {params.temp_out_sam} 2> {log}
            
                    java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                        -jar  {params.picard} \
                        SortSam \
                        I={params.temp_out_sam} \
                        O={output.final_pbfusion_mapped_bam} \
                        SORT_ORDER={params.sort_order}  2>> {log}
            
                    samtools index {output.final_pbfusion_mapped_bam} 2>> {log}

                    rm {params.temp_out_sam}
                """  
            )