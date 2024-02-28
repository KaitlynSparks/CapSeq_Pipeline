if __name__ ==  '__main__':
    import pandas as pd
    from pybiomart import Server, Dataset

    dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    #print(dataset.default_attributes)
    attribute = dataset.list_attributes()
    attribute.to_csv('attributes.csv')
    #filters = dataset.list_filters()
    #print(filters)
    #print(dataset.query(filters={'transcript_is_canonical': True,'link_ensembl_gene_id': 'ENSG00000165995'}))
