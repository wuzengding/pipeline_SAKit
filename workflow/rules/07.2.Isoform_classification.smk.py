rule Isoform_classification:
    input:
        collapsed_isoform_gtf = rules.Isoform_filter_subset.output[0],
        fl_count_file = rules.Isoform_filter_subset.output[3],
        salmon_short_reads_align = rules.Isoform_quantfs.output[0],
        #junction_shortreads = os.path.join(outpath,"results/03.Align_NGS/{0}_SJ.out.tab".format(SampleID))

    output:
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_classification.txt".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_junctions.txt".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_corrected.fasta".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_corrected.gtf".format(sample=SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{sample}.isoform_corrected.faa".format(sample=SampleID)),
    
    params:
        pythonpath = config["Env"]["python_lib"],
        aligner = "minimap2",
        chunks = 5,
        isoannot = "--isoAnnotLite",
        outprefix = ".".join([SampleID,"isoform"]),
        output_dir = os.path.join(outpath, "results/05.Isoform_Novel/"),
        isoannotlitegff3 = config["Reference"]["annotgff3"],
        genome = config["Reference"]["genome"],
        genome_annotation = config["Reference"]["annotgtf"],
        orf_input_path = os.path.join(outpath, rules.Isoform_filter_subset.output[2])

    conda:
        "SQANTI3.env"
        
    log:
        os.path.join(outpath, "log/{0}.Isoform_classification.log".format(SampleID))

    threads:
        18
    
    shell:
        """
        export PYTHONPATH=/mnt/user/wzd/03.biotools/software/cDNA_Cupcake/sequence:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake
        
        if [ -z {config[NGSShortReads][fastq_R1]} ]; then
            {config[ScriptTools][sqanti_qc]} \
                --dir {params.output_dir} \
                -t {threads} \
                --aligner_choice {params.aligner} \
                --chunks {params.chunks} \
                -o {params.outprefix} \
                {params.isoannot} \
                --gff3 {params.isoannotlitegff3} \
                --fl_count {input.fl_count_file} \
                --orf_input {params.orf_input_path} \
                {input.collapsed_isoform_gtf} {params.genome_annotation} {params.genome} >{log} 2>&1
        
        else
            {config[ScriptTools][sqanti_qc]} \
                --dir {params.output_dir} \
                -t {threads} \
                --aligner_choice {params.aligner} \
                --chunks {params.chunks} \
                -o {params.outprefix} \
                {params.isoannot} \
                --gff3 {params.isoannotlitegff3} \
                --fl_count {input.fl_count_file} \
                --orf_input {params.orf_input_path} \
                --expression {rules.Isoform_quantfs.params.quant_file} \
                --coverage "{rules.Align_NGS.params.output_prefix}_SJ.out.tab" \
                --report pdf \
                {input.collapsed_isoform_gtf} {params.genome_annotation} {params.genome} >{log} 2>&1
        
            unset PYTHONPATH
        fi
        """
    

    
    #    if not config["NGSShortReads"]["fastq_R1"]:
    #        shell(\
    #            """
    #            export PYTHONPATH=$PYTHONPATH:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake/sequence
    #            export PYTHONPATH=$PYTHONPATH:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake
    #            
    #             {config[ScriptTools][sqanti_qc]} \
    #                    --dir {params.output_dir} \
    #                    -t {threads} \
    #                    --aligner_choice {params.aligner} \
    #                    --chunks {params.chunks}  \
    #                    -o {params.outprefix} \
    #                    {params.isoannot} \
    #                    --gff3 {params.isoannotlitegff3}  \
    #                    --fl_count {input.fl_count_file}  \
    #                    --orf_input {params.orf_input_path} \
    #                    {input.collapsed_isoform_gtf} {params.genome_annotation}  {params.genome} >{log} 2>&1
    #
    #            unset PYTHONPATH
    #            """
    #            )
    #    else:
    #        junction_shortreads = rules.Align_NGS.params.output_prefix + "_SJ.out.tab"
    #        expression = rules.Isoform_quantfs.params.quant_file
    #        shell(\
    #            """
    #            export PYTHONPATH=$PYTHONPATH:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake/sequence
    #            export PYTHONPATH=$PYTHONPATH:/mnt/user/wzd/03.biotools/software/cDNA_Cupcake
    #            
    #            {config[ScriptTools][sqanti_qc]} \
    #               --dir {params.output_dir} \
    #               -t {threads} \
    #               --aligner_choice {params.aligner} \
    #               --chunks {params.chunks}  \
    #               -o {params.outprefix} \
    #               {params.isoannot} \
    #               --gff3 {params.isoannotlitegff3}  \
    #               --fl_count {input.fl_count_file}  \
    #               --orf_input {params.orf_input_path} \
    #               --expression {expression} \
    #               --coverage {junction_shortreads} \
    #               --report pdf \
    #               {input.collapsed_isoform_gtf} {params.genome_annotation}  {params.genome}  >{log} 2>&1
    #               
    #            unset PYTHONPATH
    #            """
    #             )