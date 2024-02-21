rule Cluster:
    input:
        flnc_xml = rules.Refine.output.flnc_xml
    output:
        clusteredbam = os.path.join(outpath,"results/02.Cluster_Align/{0}.flnc.clustered.bam".format(SampleID)),
        cluster_report = os.path.join(outpath,"results/02.Cluster_Align/{0}.flnc.clustered.cluster_report.csv".format(SampleID)),
        hq_transcripts = os.path.join(outpath,"results/02.Cluster_Align/{0}.flnc.clustered.hq.fasta".format(SampleID)),
    log:
        os.path.join(outpath,"log/{0}.Clustered.log".format(SampleID))
    params:
        clustertool=config["cluster"]["software"]["clustertool"],
        poa = config["cluster"]["parameter"]["poa"],
        qvs = config["cluster"]["parameter"]["qvs"],
        outputpattern =config["cluster"]["parameter"]["outputpattern"],
        thread = config["cluster"]["parameter"]["thread"],
    shell:
        """
        {params.clustertool} cluster {input.flnc_xml} {output.clusteredbam} \
        --poa-cov {params.poa} \
        {params.qvs}  \
        {params.outputpattern} \
        -j {params.thread} >> {log} 2>&1 
        gunzip {outpath}/results/02.Cluster_Align/{SampleID}.flnc.clustered.hq.fasta.gz  >> {log}
        """