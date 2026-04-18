import os

#定义根目录与所有的数据路径
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROCESSED_DIR = os.path.join(BASE_DIR, 'data', 'process')
RAW_DATA_DIR = os.path.join(BASE_DIR, 'data', 'raw')
SPACEDUST_OUTPUT_DATA_DIR = os.path.join(RAW_DATA_DIR, 'spacedust_batch_output')

ECOLI_GOA_PATH = os.path.join(RAW_DATA_DIR, 'GO', '18.E_coli_MG1655.goa.txt')
ECOLI_STRING_ALIAS_PATH = os.path.join(RAW_DATA_DIR, 'STRING', 'E.coli.protein.aliases.v12.0.txt')
ECOLI_STRING_LINKS_PATH = os.path.join(RAW_DATA_DIR, 'STRING', '511145.protein.links.detailed.v12.0.txt')
BSUBTILIS_GOA_PATH = os.path.join(RAW_DATA_DIR, 'GO', '6.B_subtilis_168.goa.txt')
BSUBTILIS_STRING_ALIAS_PATH = os.path.join(RAW_DATA_DIR, 'STRING', 'B.subtilis.protein.aliases.v12.0.txt')
BSUBTILIS_STRING_LINKS_PATH = os.path.join(RAW_DATA_DIR, 'STRING', 'B.subtilis.protein.links.v12.0.txt')
OBO_PATH = os.path.join(RAW_DATA_DIR, 'GO', 'go-basic.obo')
FEATURE_TABLE_PATH = os.path.join(RAW_DATA_DIR, 'Ecoli_feature_table.txt')


MATCH_TABLE_PATH = os.path.join(PROCESSED_DIR, 'ecoli_match.csv')
GO_PAIRS_PATH = os.path.join(PROCESSED_DIR, 'go_pairs.csv')
GO_STRING_PAIRS_PATH = os.path.join(PROCESSED_DIR, 'go-string_match_pairs.csv')
LINKS_THRESHOLD = 700


