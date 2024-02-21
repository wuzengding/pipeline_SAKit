rule Make_cmd:
    input:
        run_cmd = os.path.join(outpath,'Run.sh')
    output:
        cmdfile = os.path.join(outpath,"cmd/{0}.all_command.sh".format(SampleID)),
        cmdlog = os.path.join(outpath,"cmd/{0}.dry_command.log".format(SampleID))
    log:
        os.path.join(outpath,'log/{0}.Make_cmd.log'.format(SampleID))
    params:
        getcmdtool = config["InhouseScript"]["getcmdtool"]
        
    shell:
        '''
        CMD="$(cat {input.run_cmd}|grep snakemake) -n"
        ${{CMD}} > {output.cmdlog} 2>{log}
        
        python {params.getcmdtool} -f {output.cmdlog} -o {outpath} -l {output.cmdfile} >> {log}
        '''
