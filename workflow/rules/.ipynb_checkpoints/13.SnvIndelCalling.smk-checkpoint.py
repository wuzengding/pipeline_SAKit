rule SnvIndelCalling:
    input:
        rules.Align_NGS.output[0]
    output:
        snv_vcf = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.single_variant.vcf".format(SampleID)),
        snpeff_res = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.single_variant.snpeff.vcf".format(SampleID)),
        report = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.single_variant.snpeff.html".format(SampleID)),
    params:
        varcallpath=config["SoftwareTools"]["vardict"],
        snpeff = config["SoftwareTools"]["snpeff"],
        snpsift = config["SoftwareTools"]["snpsift"],
        species = config["species"],
        bedfile = config["Reference"]["bedfile"],
        genome = config["Reference"]["genome"],
        af_cutoff = 0.001,
        prefix = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.single_variant".format(SampleID))
    threads:
        40
    log:
         os.path.join(outpath,"log/{0}.SnvIndelCalling.log".format(SampleID))
    shell:
        """
        echo {params.species}
        
        if [ {params.species} == "human" ]; then
            species="hg38"
        elif [ {params.species} == "mouse" ]; then
            species="mm10"
        fi

        /usr/bin/java -jar -Xmx20g {params.varcallpath}/lib/VarDict-1.8.1.jar  \
            -G {params.genome} -f {params.af_cutoff} -N {SampleID} -b {input[0]} -U -c 1 -S 2 -E 3 -g 4 {params.bedfile} -th 40| \
            {params.varcallpath}/bin/teststrandbias.R | \
            {params.varcallpath}/bin/var2vcf_valid.pl -N {SampleID} -E -f {params.af_cutoff} > {output.snv_vcf}
            
        /usr/bin/java -jar -Xmx20g {params.snpeff}/snpEff.jar  \
            -c {params.snpeff}/snpEff.config -noLog -v -noInteraction  \
            -nextProt $species -stats {output.report} \
            -fastaProt {params.prefix}.prot.fa  \
            {output.snv_vcf}  \
            > {output.snpeff_res}  \
            2> {log}
        """
