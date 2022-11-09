import ntpath
import pandas as pd

def clean_salmon_output( quant_file):
    read_quant_sf = pd.read_table(quant_file , sep ="\t")
    col_name = ["target_id", "length",  "eff_length",  "est_counts",  "tpm"]
    temp = read_quant_sf[["Name", "Length", "EffectiveLength", "NumReads", "TPM"]]
    temp.columns = col_name 
    temp.to_csv(quant_file, sep ="\t")

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
        quant_file = os.path.join(outpath,"results/05.Isoform_Novel/Salmon_Isoformfiltered/quant.sf"),
        kmer = 31
    log:
        os.path.join(outpath,"log/Isoform_quantfs.log")

    threads:
        18
    
    run:
    
        shell(\
            """
                salmon index \
                -t {input.isofiltered} \
                -i {params.indexoutput}  \
                -k {params.kmer} \
                -p {threads} > {log} 2>&1
            """)
            
        shell(\
            """
                salmon --no-version-check quant \
                        --validateMappings \
                        -i {params.indexoutput} \
                        --numBootstraps {params.numBootstraps}  \
                        -p {threads} \
                        -l A  \
                        -1 {input.r1} \
                        -2 {input.r2}  \
                        -o  {params.output_dir} > {log} 2>&1
                
                #touch {output}
            """)
        
        clean_salmon_output(params.quant_file)
