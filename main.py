if __name__ ==  '__main__':
    import pandas as pd
    from ucsc_genomes_downloader import Genome
    import pybedtools
    from pybiomart import Server, Dataset
    import requests, sys
    
    #GET GENE INTERVALS

    gene_list = ['ENSG00000165995',
                 'ENSG00000141837',
                 'ENSG00000148408',
                 'ENSG00000151067',
                 'ENSG00000151067']
    
    query_attributes=['chromosome_name', 
                      'start_position',
                      'end_position',
                      'ensembl_gene_id']
    
    query_filters={'transcript_is_canonical': True,
            'link_ensembl_gene_id': gene_list}

    ensembl_dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    
    gene_region_df = ensembl_dataset.query(attributes=query_attributes, filters = query_filters)
    gene_region_list = gene_region_df.values.tolist()
    print(gene_region_list)
    #pybedtools.create_interval_from_list




    #for gene_id in gene_list:

"""     server = "https://rest.ensembl.org"
    ext = "/xrefs/symbol/homo_sapiens/" + gene_id +"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    hg38 = Genome(assembly="hg38")
    adult_DLPFC_bedtool = pybedtools.BedTool("../data/adult_DLPFC_coloured_by_expr_str.bed")
    gene_bed = 
    sequence_adult_DLPFC = sequence(fi=fasta)
    adult_DLPFC = pd.read_csv("data/adult_DLPFC_coloured_by_expr_str.bed", sep="\t", names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
    fetal_DLPFC = pd.read_csv("data/fetal_DLPFC_coloured_by_expr_str.bed", sep="\t", names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
    Caud = pd.read_csv("data/CAUD_coloured_by_expr_str.bed", sep="\t", names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
    Hip = pd.read_csv("data/HIP_coloured_by_expr_str.bed", sep="\t", names = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])

    #df_list = ['adult_DLPFC_filtered', 'fetal_DLPFC_filtered', 'Caud_filtered', 'Hip_filtered']
    transcripts = adult_DLPFC[adult_DLPFC.name.str.contains("ENST00000327702.12|ENST00000347598.9|ENST00000399641.6")]
    sequences = hg38.bed_to_sequence(transcripts) """
