
rule SV_calling:
    input:
        hq_mapped_hg38_bam = rules.Align.output.bamsort_tx

    output:
        pbsv_svsig = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.svsig.gz".format(SampleID)),
        pbsv_vcf = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.var.vcf".format(SampleID)),
        pbsv_vcf_annot = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.annot.var.vcf".format(SampleID)),
        report = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.html".format(SampleID)),
        pbsv_vcf_annot_tsv = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.annot.var.tsv".format(SampleID)),
    
    params:
        genome = config["Reference"]["genome"],
        min_ref_span= config["PBsvParam"]["min_ref_span"],
        call_min_read_perc_one_sample = config["PBsvParam"]["call_min_read_perc_one_sample"],
        ccs_flag = "--ccs",
        snpeff_java_mem = "-Xmx8g",
        snpeff = config["SoftwareTools"]["snpeff"],
        snpsift = config["SoftwareTools"]["snpsift"],
        species = config["Reference"]["species"],
        sep = "';'",
        empty = "'.'"



    threads:
        16

    log:
         os.path.join(outpath,"log/SV_calling.log")

    run:
        if params.species == 'hs':
            species = "hg38"
        elif  params.species == 'mm':
            species = "mm10"

        shell("""
            pbsv discover \
                -m {params.min_ref_span} \
                {input.hq_mapped_hg38_bam} {output.pbsv_svsig} > {log} 2>&1

            pbsv call  \
            -P {params.call_min_read_perc_one_sample} \
            {params.ccs_flag} \
            -j {threads} \
            {params.genome}  {output.pbsv_svsig}  {output.pbsv_vcf} >> {log} 2>&1

            java {params.snpeff_java_mem} \
                 -jar {params.snpeff} \
                 -v -stats {output.report} {species} {output.pbsv_vcf} > {output.pbsv_vcf_annot} 2> {log}
            
            java {params.snpeff_java_mem} \
                 -jar {params.snpsift} extractFields \
                 -s {params.sep} -e {params.empty} {output.pbsv_vcf_annot}  CHROM POS END SVLEN SVTYPE FILTER ANN[*].GENE ANN[*].IMPACT ANN[*].EFFECT ANN[*].RANK > {output.pbsv_vcf_annot_tsv} 2> {log}
            

        """)
