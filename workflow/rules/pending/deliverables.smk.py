rule Novel_protein:
    input:
        
    
    
    
        collapsed_gff = rules.Isoform_calling.output[0],
        collapsed_rep_fa = rules.Isoform_calling.output[1],
        collapsed_group = rules.Isoform_calling.output[2],
        collapsed_abundance = rules.Isoform_filter_bycounts.output[1],
        cupcake_ignored_ids = os.path.join(outpath,"results/03.Collapse/cupcake.ignored_ids.txt"),
        collapsed_gff_fl = rules.Isoform_filter_bycounts.output[2],
        collapsed_rep_fa_fl =  rules.Isoform_filter_bycounts.output[4],
        collapsed_abundance_fl =  rules.Isoform_filter_bycounts.output[3],
        collapsed_gff_fl_away = rules.Isoform_filter_subset.output[0],
        collapsed_rep_fa_fl_away =  rules.Isoform_filter_subset.output[1],
        collapsed_abundance_fl_away =  rules.Isoform_filter_subset.output[3],
        #collapsed_filtered_hq_bam_tx  = rules.Isoform_calling._filtered_hq_mapping.output.collapsed_filtered_hq_bam_tx
        sqanti3_class = rules.Isoform_classification.output[0],
        sqanti3_junc = rules.Isoform_classification.output[1],
        sqanti3_corrected_fa = rules.Isoform_classification.output[2],
        sqanti3_corrected_gtf = rules.Isoform_classification.output[3],
        sqanti3_class_flt = rules.Isoform_classification_filter.output[4],
        sqanti3_junc_flt = rules.Isoform_classification_filter.output[3],
        sqanti3_corrected_fa_flt = rules.Isoform_classification_filter.output[5],
        sqanti3_corrected_gtf_flt = rules.Isoform_classification_filter.output[2],
        sqanti3_corrected_reason_flt = rules.Isoform_classification_filter.output[0],
        sqanti3_corrected_report_flt = rules.Isoform_classification_filter.output[1],
        pbsv_svsig = rules.SV_calling.output[0],
        pbsv_vcf =  rules.SV_calling.output[1],
        pbsv_vcf_annot = rules.SV_calling.output[2],
        toi_list = rules.run_toi_list.output[0],
        pbsv_result = rules.run_toi_list.output[1],
        sqanti3_classification_final_results = rules.Isoform_summary.output[0],
        sqanti3_classification_filtered_out = rules.Isoform_summary.output[1],
        dis_gene_asso = rules.run_dis_gene_asso.output[0],
        # Fusion Pipeline Outputs
        fusion_gff = rules.Fusion_calling.output[0],
        fusion_rep = rules.Fusion_calling.output[1],
        fusion_abundance = rules.Fusion_calling.output[2],
        hq_isoforms_mapped = rules.Align.output.samsort_tx,
        hq_isoforms_mapped_bai = os.path.join(outpath,"results/FinalResults/{0}_hq_isoforms_mapped_hg38_sorted.bam.bai".format(SampleID)),
        fusion_class = rules.Fusion_classification.output[0],
        fusion_rep_corrected = os.path.join(outpath,"results/FusionClassification/lq_isoforms.fasta.fusion_corrected.fasta"),
        fusion_gtf_corrected = os.path.join(outpath,"results/FusionClassification/lq_isoforms.fasta.fusion_corrected.gtf"),
        refAnnotation_lq_isoforms = os.path.join(outpath,"results/FusionClassification/refAnnotation_lq_isoforms.fasta.fusion.genePred"),
        final_pbfusion_mapped_bam =  rules.run_filter_pbfusion_to_bam.output[1],
        final_pbfusion_mapped_bam_bai = os.path.join(outpath,"results/FinalResults/{0}_final_pbfusion_mapped_hg38_sorted.bam.bai".format(SampleID)),
        fusion_annot =  rules.run_fusion_collect_info.output[0],
        fusion_annot_fltrd = os.path.join(outpath,"results/06.Fusion_calling/lq_isoforms.fasta.fusion.annotated_ignored.txt"),
        fusion_classification_final_results =  rules.run_fusion_events.output[0],
        singlevariantvcf = rules.singlevariant.output[0],
        ccs_fl_bam = rules.run_mapping_css_fastq.output[1]


    output:
        collapsed_gff = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.gff"),
        collapsed_rep_fa = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.rep.fa"),
        collapsed_group = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.group.txt"),
        collapsed_abundance = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.abundance.txt"),
        cupcake_ignored_ids = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.ignored_ids.txt"),
        collapsed_gff_fl = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.gff".format(filter_count_cutoff = cutoff)),
        collapsed_rep_fa_fl = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.rep.fa".format(filter_count_cutoff = cutoff)),
        collapsed_abundance_fl = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.abundance.txt".format(filter_count_cutoff = cutoff)),
        collapsed_gff_fl_away = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.gff".format(filter_count_cutoff = cutoff)),
        collapsed_rep_fa_fl_away = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.rep.fa".format(filter_count_cutoff = cutoff)),
        collapsed_abundance_fl_away = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.abundance.txt".format(filter_count_cutoff = cutoff)),
        sqanti3_class = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.txt"),
        sqanti3_junc = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_junctions.txt"),
        sqanti3_corrected_fa = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_corrected.fasta"),
        sqanti3_corrected_gtf = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_corrected.gtf"),
        sqanti3_class_flt = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_classification.txt"),
        sqanti3_junc_flt = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_junctions.txt"),
        sqanti3_corrected_fa_flt = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite.fasta"),
        sqanti3_corrected_gtf_flt = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite.gtf"),
        sqanti3_corrected_reason_flt = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_reasons.txt"),
        sqanti3_corrected_report_flt = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_SQANTI3_report.pdf"),
        pbsv_svsig = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step3/{0}_hq_isoforms.svsig.gz".format(SampleID)),
        pbsv_vcf = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step3/{0}_hq_isoforms.var.vcf".format(SampleID)),
        pbsv_vcf_annot =  os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step3/{0}_hq_isoforms.annot.var.vcf".format(SampleID)),
        toi_list = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step4/{0}_toi_list.tsv".format(SampleID)),
        pbsv_result = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step4/{0}_pbsv_results.tsv".format(SampleID)),
        sqanti3_classification_final_results = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step4/{0}_sqanti3_canonical_results.tsv".format(SampleID)),
        sqanti3_classification_filtered_out = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step4/{0}_sqanti3_noncanonical_results.tsv".format(SampleID)),
        dis_gene_asso = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step5/gene_dis_ass_results.tsv"),
        fusion_gff = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step1/lq_isoforms.fasta.fusion.gff"),
        fusion_rep = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step1/lq_isoforms.fasta.fusion.rep.fa"),
        fusion_abundance =  os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step1/lq_isoforms.fasta.fusion.abundance.txt"),
        hq_isoforms_mapped = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step2/{0}_hq_isoforms_mapped_hg38_sorted.bam".format(SampleID)),
        hq_isoforms_mapped_bai = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step2/{0}_hq_isoforms_mapped_hg38_sorted.bam.bai".format(SampleID)),
        fusion_class = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step2/lq_isoforms.fasta.fusion_classification.txt"),
        fusion_rep_corrected = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step2/lq_isoforms.fasta.fusion_corrected.fasta"),
        fusion_gtf_corrected = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step2/lq_isoforms.fasta.fusion_corrected.gtf"),
        refAnnotation_lq_isoforms = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step2/refAnnotation_lq_isoforms.fasta.fusion.genePred"),
        final_pbfusion_mapped_bam = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step3/{0}_final_pbfusion_mapped_hg38_sorted.bam".format(SampleID)),
        final_pbfusion_mapped_bam_bai = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step3/{0}_final_pbfusion_mapped_hg38_sorted.bam.bai".format(SampleID)),
        fusion_annot =  os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step3/lq_isoforms.fasta.fusion.annotated.txt"),
        fusion_annot_fltrd = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step3/lq_isoforms.fasta.fusion.annotated_ignored.txt"),
        fusion_classification_final_results = os.path.join(outpath,"results/Deliverables/Fusion/Fusion_Step4/{0}_fusion_classification_final_results_fusionhub.tsv".format(SampleID)),
        singlevariantvcf = os.path.join(outpath,"results/Deliverables/SNV_INDEL/{0}.singlevar.vcf".format(SampleID))
    params:
        junction_shortreads_out = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2.1/{0}SJ.out.tab".format(SampleID)),
        expression_out = os.path.join(outpath,"results/Deliverables/Isoform/Isoform_Step2.1/cupcake.collapsed.filtered.rep.renamed.fasta.salmon")


    threads:
        1

    log:
        os.path.join(outpath,"log/run_deliverables.log")

    run:
        for no in range(len(output) -1):
            orig_file = input[no]
            dest_file = output[no]
            shell(\
                    """
                    ln -s {orig_file} {dest_file}
                    """ \
            )
        
        # Check if short reads are processed, then link STAR Junc and Salmon results to Isoform_Step2.1
        step2_1_dir = os.path.dirname(params.junction_shortreads_out)
        shell("mkdir {step2_1_dir}")
        if  config["ILLUMINASHORTREADS"]["ill_fastq_R1"]:
            junction_shortreads = rules.Align_NGS.params.output_prefix + "SJ.out.tab"
            expression = rules.Isoform_quantfs.params.output_dir
            shell(\
                    """  
                        ln -s {junction_shortreads} {params.junction_shortreads_out} 
                        ln -s {expression} {params.expression_out} 
                    """ \
                )
        else:
            shell(\
                    """
                        touch {params.junction_shortreads_out} 
                        mkdir {params.expression_out} 
                    """ \
                )


