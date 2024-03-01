
import pandas as pd
from ucsc_genomes_downloader import Genome
import pybedtools
from pybiomart import Server, Dataset
import requests, sys

def gene_filter():
    #remove "test" to run proper analysis
    with open('input/genes_of_interest.txt') as f:
        lines = f.read().splitlines()
    gene_ensembl_ID_list = lines
    
    #specify info to get from query
    query_attributes=['chromosome_name', 
                    'start_position',
                    'end_position',
                    'strand']
    
    #specify what to filter query on
    query_filters={'transcript_is_canonical': True,
                   'link_ensembl_gene_id': gene_ensembl_ID_list}
    
    #set what database to use
    ensembl_dataset = Dataset(name='hsapiens_gene_ensembl',
                              host='http://www.ensembl.org')
    
    #actually submit query
    gene_region_df = ensembl_dataset.query(attributes=query_attributes,
                                           filters = query_filters)
        
    #add chr before chromosome number
    def chr(x):
        return 'chr' + str(x)
    gene_region_df['Chromosome/scaffold name'] = gene_region_df['Chromosome/scaffold name'].apply(chr) 
    
    #make strandedness a string so that it doesnt mess everything up later lol
    def strandedness_translator(x):
        return '+' if x==1 else '-'
        
    gene_region_df['Strand'] = gene_region_df['Strand'].apply(strandedness_translator)

    #create mask bedfile from gene regions
    gene_region_list = gene_region_df.values.tolist()
    joint_bedtool = pybedtools.BedTool(' '.join(map(str,gene_region_list)), from_string=True)
    
    #below can be more efficient
    for i in range(1,len(gene_region_list)):
        single_bedtool = pybedtools.BedTool(' '.join(map(str,gene_region_list[i])), from_string=True)
        joint_bedtool = joint_bedtool.cat(single_bedtool, force_truncate=True)
        
    #sort and save mask bed file
    gene_masking_bedtool = joint_bedtool.sort(genome="hg38")

    gene_masking_bedtool.saveas('temp_files/gene_masking_bedfile.bed')
    return gene_masking_bedtool

def transcript_filter(bed_file):
    with open('input/transcripts_of_interest.txt') as f:
        lines = f.read().splitlines()
    transcript_ensemble_ID_list = lines
    transcript_and_gene_masked_bed = bed_file[bed_file.iloc[:,3].str.contains('|'.join(transcript_ensemble_ID_list), na=False)]
    return transcript_and_gene_masked_bed










"""     server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/" + gene_id +"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    hg38 = Genome(assembly="hg38")
    sequences = hg38.bed_to_sequence(transcripts) """


#perhaps useful in future
'''
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
'''