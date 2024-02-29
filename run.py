if __name__ ==  '__main__':
    import capseq_pipeline
    import glob
    import pybedtools
    import pandas as pd
    import os

    filter_for_results = True
    filter_for_genes = True
    filter_for_transcripts = True
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
            results_bedtool = gene_masking_bedtool.intersect(results_bedtool, F = 1.0, wb=True)
        results_bedtool.saveas(gene_masked_bedfile_name)    
        
        gene_masked_bed = pd.read_csv(gene_masked_bedfile_name, sep = '\t')
        transcript_and_gene_masked_bed = gene_masked_bed
        if filter_for_transcripts:
            transcript_and_gene_masked_bed = capseq_pipeline.transcript_filter(transcript_and_gene_masked_bed)
        
        transcript_and_gene_masked_bed.to_csv('temp_files/transcript_and_' + gene_masked_bedfile_name.split('/')[-1],sep='\t')


if clear_temp:
    files = glob.glob('temp_files/*')
    for f in files:
        os.remove(f)

    
    

    
    




