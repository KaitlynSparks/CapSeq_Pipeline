if __name__ ==  '__main__':
    
    #PREAMBLE
    import pandas as pd
    from ucsc_genomes_downloader import Genome
    import pybedtools
    from pybiomart import Server, Dataset
    import requests, sys
    import glob

    '''
    #GET GENE INTERVALS
            
        #create list of gene ensembl IDs from genes_of_interest.txt file
    with open('data/input/genes_of_interest.txt') as f:
        lines = f.read().splitlines()
    gene_ensembl_ID_list = lines
    
        #specify info to get from query
    query_attributes=['chromosome_name', 
                      'start_position',
                      'end_position']
    
        #specify what to filter query on
    query_filters={'transcript_is_canonical': True,
            'link_ensembl_gene_id': gene_ensembl_ID_list}
    
        #set what database to use
    ensembl_dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    
        #actually submit query
    gene_region_df = ensembl_dataset.query(attributes=query_attributes, filters = query_filters)
        
        #add chr before chromosome number
    def chr(x):
        return 'chr' + str(x)
    gene_region_df['Chromosome/scaffold name'] = gene_region_df['Chromosome/scaffold name'].apply(chr) 

        #create mask bedfile from gene regions
    gene_region_list = gene_region_df.values.tolist()
    joint_bedtool = pybedtools.BedTool(' '.join(map(str,gene_region_list[0])), from_string=True)
    for i in range(1,len(gene_region_list)):
        single_bedtool = pybedtools.BedTool(' '.join(map(str,gene_region_list[i])), from_string=True)
        joint_bedtool = joint_bedtool.cat(single_bedtool, force_truncate=True)
        
        #sort and save mask bed file
    joint_bedtool = joint_bedtool.sort(genome="hg38")    
    joint_bedtool.saveas('data/temp_files/masking_bedfile.bed')
    '''

    '''
    #MASK RESULTS FOR GENES OF INTEREST

        #if masking file saved and above code commented out, run this:
    
    masking_bedtool = pybedtools.BedTool('data/temp_files/masking_bedfile.bed')

        #otherwise, run:
    
    #masking_bedfile = joint_bedtool

        #then run:
    with open('data/input/CaptureSEQ_results_of_interest.txt') as f:
        lines = f.read().splitlines()
    capseq_results_bedfile_list = lines

    headers = ['chrom', 
               'chromStart',
               'chromEnd',
               'name',
               'score',
               'strand',
               'thickStart',
               'thickEnd',
               'itemRgb',
               'blockCount',
               'blockSizes',
               'blockStarts']
    
    for results_bedfile in capseq_results_bedfile_list:
        results_bedtool = pybedtools.BedTool('data/input/CaptureSEQ_results/' + results_bedfile)
        results_bedtool = results_bedtool.sort(genome="hg38")
        masked_bedtool = masking_bedtool.intersect(results_bedtool, F = 1.0, wb=True)
        masked_bedtool.saveas('data/temp_files/masked_' + results_bedfile)
    '''

    #FILTER MASKED RESULTS FOR TRANSCRIPTS OF INTEREST
    
        #get transcript list
    with open('data/input/transcripts_of_interest.txt') as f:
        lines = f.read().splitlines()
    transcript_ensemble_ID_list = lines







    #for gene_id in gene_list:

"""     server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/" + gene_id +"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    hg38 = Genome(assembly="hg38")
    sequences = hg38.bed_to_sequence(transcripts) """
