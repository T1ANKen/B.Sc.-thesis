import pandas as pd
from goatools.obo_parser import GODag


#打包处理数据的函数
def load_basic(path):
    print('引入go-basic数据...')
    return GODag(path)

def load_goa(path):
    print('引入annotation数据...')
    columns = [
        'DB', 'UniProt_ID', 'Gene_Symbol', 'Qualifier', 'GO_ID',
        'DB_Ref', 'Evidence', 'With_From', 'Aspect', 'Name',
        'Synonym', 'Type', 'Taxon', 'Date', 'Assigned_By', 'Ext', 'Product_ID'
    ]
    return pd.read_csv(path, sep='\t', comment='!', names=columns, dtype=str)

def load_alias(path):
    print('引入alias数据...')
    return pd.read_csv(path, sep='\t')

def load_string_links(path):
    print('引入links数据...')
    return pd.read_csv(path, sep=' ')

def load_query_lookup(path):
    print('引入lookup数据...')
    return pd.read_csv(path, sep="\t", header=None, 
                   names=['seqid', 'header', 'setid'])

def load_feature_table(path):
    print('引入NCBI feature table数据...')
    return pd.read_csv(path, sep="\t", low_memory=False)

