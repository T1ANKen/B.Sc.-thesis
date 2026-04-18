import pandas as pd
import numpy as np
import config
import os

#从特定物种的query.lookup中，添加每个基因的ID，用来连接conservation_matrix和match_pairs

def query_mapping(df_query, df_feature_table, taxon_id):
    #从右切 3 刀，和spacedust测试逻辑相同
    temp = df_query['header'].str.rsplit("_", n=3, expand=True)
    df_query['idx'] = temp[1].astype(int)
    raw_start = temp[2].astype(int)
    raw_end = temp[3].astype(int)
    df_query['start'] = np.minimum(raw_start, raw_end)
    df_query['end'] = np.maximum(raw_start, raw_end)
    df_query['pred_strand'] = np.where(raw_start < raw_end, '+', '-')

    #只关心蛋白质编码基因 (CDS)
    cds_df = df_feature_table[df_feature_table['# feature'] == 'CDS'].copy()

    #提取核心列并清除空值
    cds_mapping = cds_df[['start', 'end', 'locus_tag', 'symbol', 'product_accession']].copy()
    cds_mapping = cds_mapping.dropna(subset=['locus_tag'])
    cds_mapping['start'] = cds_mapping['start'].astype(int)
    cds_mapping['end'] = cds_mapping['end'].astype(int)


    df_q_sorted = df_query.sort_values('end')
    cds_sorted = cds_mapping.sort_values('end')

    #核心:允许 end 坐标存在最多 30 个碱基的漂移 (tolerance=30)
    #direction='nearest' 意味着找物理位置最接近的那个官方基因
    mapped_q = pd.merge_asof(
        df_q_sorted, 
        cds_sorted, 
        on='end',
        direction='nearest',
        tolerance=30
    )
    mapped_q = mapped_q.sort_values('seqid').reset_index(drop=True)
    matched_count = mapped_q['locus_tag'].notna().sum()
    total_count = len(mapped_q)
    print(f"匹配完成！共计 {total_count} 个基因，成功匹配上 NCBI 官方 ID 的有 {matched_count} 个。")
    print(f"匹配率: {(matched_count/total_count)*100:.2f}%\n")

    #生成匹配的locus tag：
    #只给成功匹配的基因加前缀，失败的保持空值
    #mapped_q.loc[mapped_q['locus_tag'].notna(), 'STRING_ID'] = "511145." + mapped_q['locus_tag'].astype(str)
    mapped_q.loc[mapped_q['locus_tag'].notna(), 'STRING_ID'] = \
        f"{taxon_id}." + mapped_q['locus_tag'].astype(str)

    print("预览前 5 行合并后的lookup：")
    print(mapped_q[['seqid', 'header', 'locus_tag', 'symbol', 'STRING_ID', 'start_x', 'start_y']].head())
    return mapped_q





def query_mapping(df_query, df_feature_table, taxon_id, df_links):
    #从右切 3 刀，和spacedust测试逻辑相同
    temp = df_query['header'].str.rsplit("_", n=3, expand=True)
    df_query['idx'] = temp[1].astype(int)
    raw_start = temp[2].astype(int)
    raw_end = temp[3].astype(int)
    df_query['start'] = np.minimum(raw_start, raw_end)
    df_query['end'] = np.maximum(raw_start, raw_end)
    df_query['pred_strand'] = np.where(raw_start < raw_end, '+', '-')
 
    #只关心蛋白质编码基因 (CDS)
    cds_df = df_feature_table[df_feature_table['# feature'] == 'CDS'].copy()
 
    #提取核心列并清除空值
    cds_mapping = cds_df[['start', 'end', 'locus_tag', 'symbol', 'product_accession']].copy()
    cds_mapping = cds_mapping.dropna(subset=['locus_tag'])
    cds_mapping['start'] = cds_mapping['start'].astype(int)
    cds_mapping['end'] = cds_mapping['end'].astype(int)
    
    # 对于特定物种（C. difficile 272563, E. lenta 479437, C. jejuni 192222, S. elongatus 1140），需要从 gene 行的 attributes 列获取 old_locus_tag
    # 注意：taxon_id 是字符串类型，所以要用字符串比较
    if str(taxon_id) in ('479437', '272563', '192222', '1140'):
        import re
        # 从 gene 行提取 old_locus_tag 映射
        gene_rows = df_feature_table[df_feature_table['# feature'] == 'gene'].copy()
        old_locus_map = {}
        for _, row in gene_rows.iterrows():
            attrs = row.get('attributes', '')
            if pd.notna(attrs):
                match = re.search(r'old_locus_tag=([^;]+)', str(attrs))
                if match:
                    old_locus = match.group(1).strip().strip("'\"")
                    # 用 locus_tag 作为 key，因为 CDS 和 gene 共享相同的 locus_tag
                    old_locus_map[row['locus_tag']] = old_locus
        
        # 将 old_locus_tag 映射到 cds_mapping
        cds_mapping['old_locus_tag'] = cds_mapping['locus_tag'].map(old_locus_map)
        print(f"  [querymapper] taxon_id={taxon_id}: 从 gene 行提取了 {cds_mapping['old_locus_tag'].notna().sum()} 个 old_locus_tag")
 
    df_q_sorted = df_query.sort_values('end')
    cds_sorted = cds_mapping.sort_values('end')
 
    #核心:允许 end 坐标存在最多 30 个碱基的漂移 (tolerance=30)
    #direction='nearest' 意味着找物理位置最接近的那个官方基因
    mapped_q = pd.merge_asof(
        df_q_sorted, 
        cds_sorted, 
        on='end',
        direction='nearest',
        tolerance=30
    )
    mapped_q = mapped_q.sort_values('seqid').reset_index(drop=True)
    matched_count = mapped_q['locus_tag'].notna().sum()
    total_count = len(mapped_q)
    print(f"匹配完成！共计 {total_count} 个基因，成功匹配上 NCBI 官方 ID 的有 {matched_count} 个。")
    print(f"匹配率: {(matched_count/total_count)*100:.2f}%\n")
 
    # 检测是否需要去掉 locus_tag 中的下划线，与 data_process.py 的逻辑保持一致：
    # NCBI 新版 locus_tag 用下划线（如 BSU_00010），
    # 但 STRING 数据库沿用旧版格式（如 BSU00010，无下划线）。
    # 通过检查 STRING links 里的 protein ID 格式来自动判断。
    remove_underscore = False
    if df_links is not None and len(df_links) > 0:
        has_underscore_in_links = df_links['protein1'].astype(str).str.contains('_').any()
        has_underscore_in_ncbi = cds_mapping['locus_tag'].astype(str).str.contains('_').any()
        remove_underscore = has_underscore_in_ncbi and not has_underscore_in_links
        print(f"  [querymapper] NCBI locus_tag 含下划线：{has_underscore_in_ncbi}")
        print(f"  [querymapper] STRING links 含下划线：{has_underscore_in_links}")
        print(f"  [querymapper] 需要去掉下划线：{remove_underscore}")
 
    # 生成匹配的 STRING_ID：用传入的 taxon_id，不再硬编码 511145
    # 对于 E. lenta (taxon_id=479437)，使用 old_locus_tag 而不是 locus_tag
    def make_string_id(row):
        if str(taxon_id) in ('479437', '272563', '1140') and 'old_locus_tag' in row and pd.notna(row.get('old_locus_tag')):
            tag = str(row['old_locus_tag'])
        else:
            tag = str(row['locus_tag'])
            if remove_underscore:
                tag = tag.replace('_', '')
        return f"{taxon_id}.{tag}"
 
    mapped_q.loc[mapped_q['locus_tag'].notna(), 'STRING_ID'] = \
        mapped_q.loc[mapped_q['locus_tag'].notna()].apply(make_string_id, axis=1)
 

    return mapped_q




