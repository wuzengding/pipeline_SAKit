import ntpath
import pandas as pd


rule Align_NGS:
    input:
        r1 = config["NGSShortReads"]["fastq_R1"],
        r2 = config["NGSShortReads"]["fastq_R2"]

    output:
        os.path.join(outpath,"results/03.Align_NGS/{0}_short_reads_genome_alignment.status".format(SampleID)),
        os.path.join(outpath,"results/03.Align_NGS/{0}_Aligned.sortedByCoord.out.bam".format(SampleID)),
        os.path.join(outpath,"results/03.Align_NGS/{0}_SJ.out.tab".format(SampleID))
    
    params:
        star_aligner_parameters = "--genomeLoad NoSharedMemory \
                                        --outSAMtype BAM SortedByCoordinate \
                                        --limitBAMsortRAM 100000000000 \
                                        --outSAMunmapped Within KeepPairs \
                                        --twopassMode Basic \
                                        --readFilesCommand zcat",
        star_index = config["Reference"]["star_index"],
        output_prefix = os.path.join(outpath,"results/03.Align_NGS/{0}_".format(SampleID))

    log:
        os.path.join(outpath,"log/{0}.short_reads_genome_alignment.log".format(SampleID))

    threads:
        18

    shell:
        """
            {config[SoftwareTools][short_reads_align]} {params.star_aligner_parameters} \
                --genomeDir {params.star_index} \
                --runThreadN {threads} \
                --outFileNamePrefix {params.output_prefix} \
                --readFilesIn {input.r1} {input.r2} > {log} 2>&1
            {config[SoftwareTools][samtools]} index {output[1]}
            touch {output}
        """
