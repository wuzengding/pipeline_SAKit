cutoff = config["isoform_filter_bycounts"]["parameter"]["cutoff"]
rule Isoform_filter_subset:
    input:
        rules.Isoform_filter_bycounts.output[2],
        rules.Isoform_filter_bycounts.output[3],
        rules.Isoform_filter_bycounts.output[4]
    
    output:
        os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.filtered.gff".format(sample=SampleID, cutoff = cutoff)),
        os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.filtered.rep.fa".format(sample=SampleID, cutoff = cutoff)),
        os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.filtered.rep.renamed.fasta".format(sample=SampleID, cutoff = cutoff)),
        os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.filtered.abundance.txt".format(sample=SampleID, cutoff = cutoff))
    
    params:
        pythonpath = config["Env"]["python_lib"],
        filter_away_subset = config["ScriptTools"]["filter_away_subset"],
        prefix = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{filter_count_cutoff}".format(sample=SampleID, filter_count_cutoff = cutoff )),
        sed = "'s/|.*//g'"

    log:
        os.path.join(outpath,"log/isoform_filter_subset.log")

    threads:
        1
    
    shell:
        """
            {config[Env][python3]} {params.filter_away_subset} {params.prefix}> {log} 2>&1
            sed {params.sed} {output[1]}> {output[2]}
        """