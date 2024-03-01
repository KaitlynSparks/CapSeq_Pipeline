if __name__ ==  '__main__':
    import capseq_pipeline
    import glob
    import pybedtools
    import pandas as pd
    import os
    from ucsc_genomes_downloader import Genome

    filter_for_results = True
    filter_for_genes = True
    filter_for_transcripts = False
    clear_temp = False

    if filter_for_results:
        with open('input/CaptureSEQ_results_of_interest.txt') as f:
            lines = f.read().splitlines()
        results_bedfile_name_list = lines
    else:
        results_bedfile_name_list = glob.glob('data/*.bed').split('/')[-1]

    for results_bedfile_name in results_bedfile_name_list:
        results_bedtool = pybedtools.BedTool('data/CaptureSEQ_results/' + results_bedfile_name)
        gene_masked_bedfile_name = 'temp_files/gene_masked_' + results_bedfile_name
        
        results_bedtool = results_bedtool.sort(genome="hg38")
        if filter_for_genes:
            gene_masking_bedtool = capseq_pipeline.gene_filter()
            results_bedtool = results_bedtool.intersect(gene_masking_bedtool, f=0.9, u=True, ) #s=True
        results_bedtool.saveas(gene_masked_bedfile_name)    
        
        gene_masked_bed = pd.read_csv(gene_masked_bedfile_name, sep = '\t')
        transcript_and_gene_masked_bed = gene_masked_bed
        if filter_for_transcripts:
            transcript_and_gene_masked_bed = capseq_pipeline.transcript_filter(transcript_and_gene_masked_bed)
        transcript_and_gene_masked_bed_file_name = 'temp_files/transcript_and_' + gene_masked_bedfile_name.split('/')[-1]
        transcript_and_gene_masked_bed.to_csv(transcript_and_gene_masked_bed_file_name,sep='\t', index=False)

        final_bedtool = pybedtools.BedTool(transcript_and_gene_masked_bed_file_name)
        fasta_file_name = 'temp_files/sequences/sequences_' + results_bedfile_name[0:-4] + '.fa'
        sequences = final_bedtool.sequence(fi= 'data/hg38.fa',split=True, fo=fasta_file_name)


if clear_temp:
    files = glob.glob('temp_files/*')
    for f in files:
        os.remove(f)

    
    

    
    




