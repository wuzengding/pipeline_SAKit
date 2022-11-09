rule CCS:
    input:
        subreadbam = config["PACBIOREADS"]["Subreadbam"]
    output:
        ccsxml = os.path.join(outpath,"results/01.CCS_Trim/{0}.ccs.xml".format(SampleID)),
        ccsbam = os.path.join(outpath,"results/01.CCS_Trim/{0}.ccs.bam".format(SampleID)),
        ccsreport = os.path.join(outpath,"results/01.CCS_Trim/{0}.ccs.report.txt".format(SampleID)),
        ccsmetrics = os.path.join(outpath,"results/01.CCS_Trim/{0}.ccs.zmw_metrics.json.gz".format(SampleID))
    log:
        os.path.join(outpath,"log/{0}.ccs.log".format(SampleID))
    params:
        ccstool = config["ccs"]["software"]["ccstool"],
        minpass = config["ccs"]["parameter"]["minpass"],
        toppass = config["ccs"]["parameter"]["toppass"],
        minsnr  = config["ccs"]["parameter"]["minsnr"],
        minlength = config["ccs"]["parameter"]["minlength"],
        maxlength = config["ccs"]["parameter"]["maxlength"],
        chunk  = config["ccs"]["parameter"]["chunk"],
        minrq = config["ccs"]["parameter"]["minrq"],
        thread = config["ccs"]["parameter"]["thread"]
    shell:
        """
        {params.ccstool} {input.subreadbam} {output.ccsxml} \
        --report-file {output.ccsreport} \
        --metrics-json {output.ccsmetrics} \
        --min-passes {params.minpass} \
        --top-passes {params.toppass} \
        --min-snr {params.minsnr} \
        --min-length {params.minlength} \
        --max-length {params.maxlength} \
        --min-rq {params.minrq} \
        -j {params.thread}  >> {log} 2>&1
        """