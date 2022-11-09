rule DeMultiplex:
    input:
        ccsxml = os.path.join(outpath,"results/01.CCS_Trim/{0}.ccs.xml".format(SampleID)),
        barcodefa = config["demultiplex"]["configfile"]["barcodefa"]
    output:
        demulxml = os.path.join(outpath,"results/01.CCS_Trim/{0}.demultiplex.consensusreadset.xml".format(SampleID)),
    log:
        os.path.join(outpath,"log/{0}.demultiplex.log".format(SampleID))
    params:
        #prefix = os.path.join(outpath,"{batchid}/{SampleID}/03.DEMULTIPLEX/{SampleID}.demultiplex"),
        limatool = config["demultiplex"]["software"]["limatool"],
        libraryDesign = config["demultiplex"]["parameter"]["libraryDesign"],
        readsflankedbyadapter=config["demultiplex"]["parameter"]["readsflankedbyadapter"],
        maxscoredbarcode=config["demultiplex"]["parameter"]["maxscoredbarcode"],
        maxscoredadapter=config["demultiplex"]["parameter"]["maxscoredadapter"],
        minpass=config["demultiplex"]["parameter"]["minpass"],
        minlength=config["demultiplex"]["parameter"]["minlength"],
        maxinputlength=config["demultiplex"]["parameter"]["maxinputlength"],
        badadapterratio=config["demultiplex"]["parameter"]["badadapterratio"],
        minscore=config["demultiplex"]["parameter"]["minscore"],
        seqmode=config["demultiplex"]["parameter"]["seqmode"],
        thread=config["demultiplex"]["parameter"]["thread"],
    shell:
        """
        {params.limatool} {input.ccsxml} {input.barcodefa} {output.demulxml} --max-scored-barcodes {params.maxscoredbarcode} \
        --max-scored-adapters {params.maxscoredadapter} \
        --min-passes {params.minpass}  \
        --min-length {params.minlength} \
        --max-input-length {params.maxinputlength} \
        --bad-adapter-ratio {params.badadapterratio} \
        --min-score {params.minscore} \
        {params.seqmode} \
        -j {params.thread} >> {log} 2>&1
        
        """