rule Refine:
    input:
        demulxml = os.path.join(outpath,"results/01.CCS_Trim/{0}.demultiplex.consensusreadset.xml".format(SampleID)),
        primerfa = config["refine"]["configfile"]["primerfa"]
    output:
        flnc_xml = os.path.join(outpath,"results/01.CCS_Trim/{0}.flnc.consensusreadset.xml".format(SampleID)),
        flnc_bam = os.path.join(outpath,"results/01.CCS_Trim/{0}.flnc.bam".format(SampleID)),
        flnc_fq = os.path.join(outpath,"results/01.CCS_Trim/{0}.flnc.fastq".format(SampleID))
    log:
        os.path.join(outpath,"log/{0}.Refine.log".format(SampleID))
    params:
        refinetool=config["refine"]["software"]["refinetool"],
        minpolyalength=config["refine"]["parameter"]["minpolyalength"],
        requirepolya=config["refine"]["parameter"]["requirepolya"],
        minrq = config["refine"]["parameter"]["minrq"],
        thread = config["refine"]["parameter"]["thread"],
        prefix_output = "results/01.CCS_Trim/{0}.flnc".format(SampleID),
        dont_compress = "-u"
    shell:
        """
        {params.refinetool} refine {input.demulxml} {input.primerfa} {output.flnc_xml} \
        --min-polya-length {params.minpolyalength} \
        --min-rq {params.minrq} {params.requirepolya} \
        -j {params.thread}  >> {log} 2>&1
                
        {config[SoftwareTools][index]} {output.flnc_bam} >> {log}
        {config[SoftwareTools][bam2fastq]} {params.dont_compress} \
        -o {params.prefix_output} {output.flnc_bam} >> {log}
        """