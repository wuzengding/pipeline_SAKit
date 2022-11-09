rule Align:
    input:
        hq_tx = rules.Cluster.output.hq_transcripts

    output:
        sam_tx = os.path.join(outpath,"results/02.Cluster_Align/hq_isoforms_mapped_hg38.sam".format(SampleID)),
        samsort_tx = os.path.join(outpath,"results/02.Cluster_Align/hq_isoforms_mapped_hg38.sorted.sam".format(SampleID)),
        bamsort_tx = os.path.join(outpath,"results/02.Cluster_Align/{0}_hq_isoforms_mapped_hg38_sorted.bam".format(SampleID)),
        hq_isoforms_mapped_bai = os.path.join(outpath,"results/02.Cluster_Align/{0}_hq_isoforms_mapped_hg38_sorted.bam.bai".format(SampleID))
    params:
        ## hg_isoforms.fasta alignment
        genome = config["Reference"]["genome"],
        ax = config["align"]["parameter"]["ax"], 
        secondary = config["align"]["parameter"]["secondary"],  
        o = config["align"]["parameter"]["o"], 
        b4 = config["align"]["parameter"]["b4"],
        uf = config["align"]["parameter"]["uf"],
        read_groups = "'@RG\\tID:HQ\\tSM:{sample}'".format(sample = SampleID),
        hard_clip_off = config["align"]["parameter"]["hard_clip_off"],
        picard = config["SoftwareTools"]["picard"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/igm/temp",
        sort_order = "coordinate"


    log:
        os.path.join(outpath,"log/align.log")

    threads:
        18

    shell:
        """
            {config[SoftwareTools][minimap2]} \
                -t {threads} \
                -ax {params.ax} \
                --secondary={params.secondary} \
                -R {params.read_groups} \
                {params.hard_clip_off} \
                {params.o} \
                {params.b4} \
                {params.uf} \
                {params.genome} \
                {input} > {output.sam_tx} 2> {log}
            
            java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                -jar  {params.picard} \
                SortSam \
                I={output.sam_tx} \
                O={output.samsort_tx} \
                SORT_ORDER={params.sort_order}  2> {log}
                
            samtools view  \
                -@ {threads} \
                -S -b -h \
                -o {output.bamsort_tx} {output.samsort_tx} > {log}
                
            samtools index {output.bamsort_tx} {output.hq_isoforms_mapped_bai} >> {log}
        """
