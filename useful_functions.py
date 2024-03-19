from pybiomart import Dataset

def gene_filter(id_list):
    gene_ensembl_id_list = id_list
    
    #specify info to get from query
    query_attributes=['chromosome_name', 
                    'start_position',
                    'end_position',
                    'ensembl_gene_id',
                    'external_gene_name',
                    'strand']
    
    #specify what to filter query on
    query_filters={'transcript_is_canonical': True,
                   'link_ensembl_gene_id': gene_ensembl_id_list}
    
    #set what database to use
    ensembl_dataset = Dataset(name='hsapiens_gene_ensembl',
                              host='http://www.ensembl.org')
    
    #actually submit query
    gene_region_df = ensembl_dataset.query(attributes=query_attributes,
                                           filters = query_filters)
    
    #sort by chrommosome number
    gene_region_df.sort_values('Chromosome/scaffold name', inplace=True)
    
    #add chr before chromosome number
    def add_chr(x):
        return 'chr' + str(x)
    gene_region_df['Chromosome/scaffold name'] = gene_region_df['Chromosome/scaffold name'].apply(add_chr) 
    
    #make strandedness a string so that it doesnt mess everything up later lol
    def strandedness_translator(x):
        return '+' if x==1 else '-'
    gene_region_df['Strand'] = gene_region_df['Strand'].apply(strandedness_translator)
    
    gene_region_df.to_csv('temp_files/get_sequences/gene_masking_bedfile.bed',sep='\t', header=False ,index=False)
    return

def transcript_filter(bed_file, id_list):
    transcript_ensemble_ID_list = id_list
    transcript_and_gene_masked_bed = bed_file[bed_file.iloc[:,3].str.contains('|'.join(transcript_ensemble_ID_list), na=False)]
    return transcript_and_gene_masked_bed


def chop_sequence(sequence):
    sequence = sequence[sequence.find("M"):-1]
    return sequence

