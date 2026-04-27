"""
多物种分析配置文件
定义 8 个物种的所有数据路径和输出路径
"""

import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_SPECIES_DIR = os.path.join(BASE_DIR, 'data', 'raw', 'species')
OUTPUT_DIR = os.path.join(BASE_DIR, 'data', 'output')

# 物种配置字典
# 每个物种包含：name(显示名), taxon_id(STRING 数据库 ID), 以及所有输入数据路径和输出目录
SPECIES_CONFIG = {
    'ecoli': {
        'name': 'E. coli MG1655',
        'taxon_id': '511145',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Ecoli', 'aliases', 'E.coli.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Ecoli', 'links', '511145.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Ecoli', 'goa', '18.E_coli_MG1655.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Ecoli', 'feature', 'GCF_000005845.2_ASM584v2_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Ecoli', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Ecoli', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'ecoli'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'ecoli', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'ecoli', 'figures'),
    },
    'bsubtilis': {
        'name': 'B. subtilis 168',
        'taxon_id': '224308',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Bsubtilis', 'aliases', '224308.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Bsubtilis', 'links', '224308.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Bsubtilis', 'goa', '6.B_subtilis_168.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Bsubtilis', 'feature', 'GCF_000009045.1_ASM904v1_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Bsubtilis', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Bsubtilis', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'bsubtilis'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'bsubtilis', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'bsubtilis', 'figures'),
    },
    'cdifficile': {
        'name': 'C. difficile 630',
        'taxon_id': '272563',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Cdifficile', 'aliases', '272563.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Cdifficile', 'links', '272563.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Cdifficile', 'goa', '25374.C_difficile_630.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Cdifficile', 'feature', 'GCF_000009205.2_ASM920v2_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Cdifficile', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Cdifficile', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'cdifficile'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'cdifficile', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'cdifficile', 'figures'),
    },
    'cjejuni': {
        'name': 'C. jejuni O-2',
        'taxon_id': '192222',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Cjejuni', 'aliases', '192222.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Cjejuni', 'links', '192222.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Cjejuni', 'goa', '10.C_jejuni_O-2.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Cjejuni', 'feature', 'GCF_000011865.1_ASM1186v1_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Cjejuni', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Cjejuni', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'cjejuni'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'cjejuni', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'cjejuni', 'figures'),
    },
    'elenta': {
        'name': 'E. lenta',
        'taxon_id': '479437',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Elenta', 'aliases', '479437.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Elenta', 'links', '479437.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Elenta', 'goa', '33183.E_lenta.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Elenta', 'feature', 'GCF_000024265.1_ASM2426v1_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Elenta', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Elenta', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'elenta'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'elenta', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'elenta', 'figures'),
    },
    'mtuberculosis': {
        'name': 'M. tuberculosis ATCC 25618',
        'taxon_id': '83332',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Mtuberculosis', 'aliases', '83332.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Mtuberculosis', 'links', '83332.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Mtuberculosis', 'goa', '30.M_tuberculosis_ATCC_25618.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Mtuberculosis', 'feature', 'GCF_000195955.2_ASM19595v2_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Mtuberculosis', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Mtuberculosis', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'mtuberculosis'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'mtuberculosis', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'mtuberculosis', 'figures'),
    },
    'selongatus': {
        'name': 'S. elongatus',
        'taxon_id': '1140',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Selongatus', 'aliases', '1140.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Selongatus', 'links', '1140.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Selongatus', 'goa', '4861838.S_elongatus.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Selongatus', 'feature', 'GCF_000012525.1_ASM1252v1_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Selongatus', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Selongatus', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'selongatus'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'selongatus', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'selongatus', 'figures'),
    },
    'smeliloti': {
        'name': 'S. meliloti',
        'taxon_id': '382',
        'aliases_path': os.path.join(RAW_SPECIES_DIR, 'Smeliloti', 'aliases', '382.protein.aliases.v12.0.txt'),
        'links_path': os.path.join(RAW_SPECIES_DIR, 'Smeliloti', 'links', '382.protein.links.detailed.v12.0.txt'),
        'goa_path': os.path.join(RAW_SPECIES_DIR, 'Smeliloti', 'goa', '58.R_meliloti.goa.txt'),
        'feature_table_path': os.path.join(RAW_SPECIES_DIR, 'Smeliloti', 'feature', 'GCF_000346065.1_ASM34606v1_feature_table.txt'),
        'matrix_path': os.path.join(RAW_SPECIES_DIR, 'Smeliloti', 'conservation', 'original_conservation_matrix.csv'),
        'querydb_lookup_path': os.path.join(RAW_SPECIES_DIR, 'Smeliloti', 'conservation', 'querydb.lookup'),
        'output_dir': os.path.join(OUTPUT_DIR, 'smeliloti'),
        'output_data_dir': os.path.join(OUTPUT_DIR, 'smeliloti', 'data'),
        'output_figure_dir': os.path.join(OUTPUT_DIR, 'smeliloti', 'figures'),
    },
}

# GO-STRING 匹配对文件路径（共用）
GO_STRING_PAIRS_PATH = os.path.join(BASE_DIR, 'data', 'process', 'go-string_match_pairs.csv')

# STRING 分数阈值
LINKS_THRESHOLD = 700