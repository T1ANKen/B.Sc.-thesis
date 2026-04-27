import pandas as pd

import data_loader


#定义数据预处理函数

evidence_classify = {
        # Experimental (High)
        'EXP': 'High', 'IDA': 'High', 'IPI': 'High', 'IMP': 'High', 'IGI': 'High', 'IEP': 'High', 'TAS': 'Low', 'IC': 'Low',
        # Author Statement / Curated (Medium)
        'ISS': 'Low', 'ISO': 'Low', 'ISA': 'Low',
        'ISM': 'Low', 'IGC': 'Low', 'IBA': 'Low', 'IBD': 'Low', 'IKR': 'Low',
        'IRD': 'Low', 'RCA': 'Low',
        # Electronic (Low)
        'IEA': 'Low'
    }


def recalculate_string_score(row, channels):
    """按照 STRING 官方贝叶斯公式重新合并分数"""
    E_product = 1.0
    PRIOR = 0.041  # STRING 官方随机基线概率
    for ch in channels:
        # 将分数还原为 0-1 的概率
        S_i = row[ch] / 1000.0
        # 只有当分数大于随机基线时，才认为有真正的证据力
        if S_i > PRIOR:
            # 去除先验背景音
            E_i = (S_i - PRIOR) / (1.0 - PRIOR)
            # 累乘"没有发生互作"的概率
            E_product *= (1.0 - E_i)
            
    # 计算综合证据概率 = 1 - 都不发生的概率
    E_combined = 1.0 - E_product
    
    # 最后把随机基线概率加回来
    S_final = E_combined + PRIOR * (1.0 - E_combined)
    
    # 重新放大 1000 倍，转化为整数得分
    return int(S_final * 1000)





def match_table(df_goa, df_alias):
    print('建立映射表...')
    df_alias = df_alias[df_alias['source']=='UniProt_AC']
    alias_dict = dict(zip(df_alias['alias'], df_alias['#string_protein_id']))

    mask= ~df_goa['Qualifier'].str.contains('NOT', na=False)
    df_goa = df_goa[mask][['Gene_Symbol', 'Name', 'UniProt_ID', 'GO_ID']].copy()
    df_goa.rename(columns={'Name': 'Protein_Name'}, inplace=True)




    df_goa['String_ID'] = df_goa['UniProt_ID'].map(alias_dict)


    table = df_goa.dropna(subset = ['String_ID'])
    table = table.drop_duplicates().reset_index(drop=True)

    return table

#这个函数只处理 BP 类的 GOterm，方便后续核心分析
def match_bp_table(df_goa, df_alias, df_NCBI, taxon_id, df_links=None):
    """
    建立只有 BP 的映射表，支持多物种
    
    Args:
        df_goa: GOA 数据框
        df_alias: STRING aliases 数据框
        df_NCBI: NCBI Feature Table 数据框
        taxon_id: 物种分类 ID（用于生成 STRING ID 前缀）
        df_links: STRING links 数据框（可选，用于检测 ID 格式）
    """
    print('建立只有 BP 的映射表...')
    df_alias = df_alias[df_alias['source']=='UniProt_AC']
    alias_dict = dict(zip(df_alias['alias'], df_alias['#string_protein_id']))

    mask = (df_goa['Aspect'] == 'P') & (~df_goa['Qualifier'].str.contains('NOT', na=False))
    df_goa = df_goa[mask][['Gene_Symbol', 'Name', 'UniProt_ID', 'GO_ID', 'Evidence']].copy()
    df_goa.rename(columns={'Name': 'Protein_Name'}, inplace=True)

    df_sorted = df_goa.sort_values(
        by=['Gene_Symbol', 'GO_ID', 'Evidence'],
        ascending=[True, True, False]
    )
    df_goa = df_sorted.drop_duplicates(subset=['Gene_Symbol', 'GO_ID'], keep='first')
    df_goa['String_ID'] = df_goa['UniProt_ID'].map(alias_dict)

    table = df_goa
    
    gene_rows = df_NCBI[df_NCBI['# feature'] == 'gene']
    locus_to_symbol = dict(zip(gene_rows['symbol'], gene_rows['locus_tag']))
    
    # 对于特定物种（C. difficile 272563, E. lenta 479437, C. jejuni 192222），需要从 attributes 列提取 old_locus_tag
    old_locus_to_tag = {}
    if taxon_id in (272563, 479437, 192222, 1140):
        import re
        for _, row in gene_rows.iterrows():
            attrs = row.get('attributes', '')
            if pd.notna(attrs):
                match = re.search(r'old_locus_tag=([^;]+)', str(attrs))
                if match:
                    old_locus = match.group(1).strip().strip("'\"")
                    old_locus_to_tag[old_locus] = row['locus_tag']
        print(f"  [match_bp_table] taxon_id={taxon_id}: 从 attributes 提取了 {len(old_locus_to_tag)} 个 old_locus_tag")
    
    # 检测 STRING links 的 ID 格式，决定是否需要去掉 locus_tag 中的下划线
    remove_underscore = False
    if df_links is not None and len(df_links) > 0:
        # 检查 STRING links 中的 protein1 是否包含下划线
        has_underscore_in_links = df_links['protein1'].str.contains('_').any()
        # 检查 NCBI locus_tag 是否包含下划线
        has_underscore_in_ncbi = gene_rows['locus_tag'].astype(str).str.contains('_').any()
        # 如果 NCBI 有下划线但 STRING 没有，则需要去掉下划线
        remove_underscore = has_underscore_in_ncbi and not has_underscore_in_links
        print(f"  NCBI locus_tag 含下划线：{has_underscore_in_ncbi}")
        print(f"  STRING links 含下划线：{has_underscore_in_links}")
        print(f"  需要去掉下划线：{remove_underscore}")
    
    def get_locus_tag(gene_name):
        # 对于特定物种，优先尝试使用 old_locus_tag 映射
        if taxon_id in (272563, 479437, 192222, 1140) and gene_name in old_locus_to_tag:
            name = old_locus_to_tag[gene_name]
        else:
            name = locus_to_symbol.get(gene_name, None)
        # 如果需要，去掉 locus_tag 中的下划线以匹配 STRING 格式
        if name and remove_underscore:
            name = name.replace('_', '')
        return f'{taxon_id}.{name}' if name else None

    table['String_ID_new'] = table['Gene_Symbol'].apply(get_locus_tag)
    table['String_ID'] = table['String_ID'].fillna(table['String_ID_new'])
    table = table.drop(columns=['String_ID_new'])
    table = table.dropna(subset = ['String_ID'])

    print(table.isnull().sum())
    print("\n=== 证据等级分布 (去重后) ===")
    print(table['Evidence'].value_counts())
    return table



