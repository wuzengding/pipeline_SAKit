import ntpath
import pandas as pd

rule Isoform_quantfs:
    input:
        isofiltered = rules.Isoform_filter_subset.output[2],
        r1 = config["NGSShortReads"]["fastq_R1"],
        r2 = config["NGSShortReads"]["fastq_R2"]

    output:
        quant_file = os.path.join(outpath,"results/05.Isoform_Novel/Salmon_Isoformfiltered/quant.sf")
    
    params:
        output_dir = os.path.join(outpath,"results/05.Isoform_Novel/Salmon_Isoformfiltered"),
        indexoutput = os.path.join(outpath,"results/05.Isoform_Novel/Salmon_Index"),
        numBootstraps = 50,
        kmer = 31
    log:
        os.path.join(outpath,"log/{0}.Isoform_quantfs.log".format(SampleID))

    threads:
        18
    
    shell:
        """
        salmon index \
            -t {input.isofiltered} \
            -i {params.indexoutput}  \
            -k {params.kmer} \
            -p {threads} > {log} 2>&1

        salmon --no-version-check quant \
            --validateMappings \
            -i {params.indexoutput} \
            --numBootstraps {params.numBootstraps}  \
            -p {threads} \
            -l A  \
            -1 {input.r1} \
            -2 {input.r2} \
            -o  {params.output_dir} > {log} 2>&1
            
        sed -i '1s/Name/target_id/; 1s/Length/length/; 1s/EffectiveLength/eff_length/; 1s/NumReads/est_counts/; 1s/TPM/tpm/' {output.quant_file} 
        """
