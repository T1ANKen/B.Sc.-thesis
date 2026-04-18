from itertools import combinations
import pandas as pd
import config
from collections import Counter

#核心分析：找出三类关联情况不同的蛋白对

#evidence 分类函数
evidence_classify = {
        # Experimental (High)
        'EXP': 'High', 'IDA': 'High', 'IPI': 'High', 'IMP': 'High', 'IGI': 'High', 'IEP': 'High', 'TAS': 'High', 'IC': 'High',
        # Author Statement / Curated (Medium)
        'ISS': 'Low', 'ISO': 'Low', 'ISA': 'Low',
        'ISM': 'Low', 'IGC': 'Low', 'IBA': 'Low', 'IBD': 'Low', 'IKR': 'Low',
        'IRD': 'Low', 'RCA': 'Low',
        # Electronic (Low)
        'IEA': 'Low'
    }

evidence_rank = {'High': 3, 'Medium': 2, 'Low': 1}


def analyze_links(df_match, df_links):
    print(df_links)

    groups = df_match.groupby('GO_ID')[['String_ID', 'Evidence']].apply(
        lambda x: list(zip(x['String_ID'], x['Evidence']))
    )
    print(f"共有 {len(groups)} 个 GO Term。")

    pair_data = []
    pair_counter = Counter()

    # 定义强度计算规则
        
    def get_strength(e1, e2):
        score = min(e1, e2)
        if score == 3:
            return 'High'
        elif score == 2:
            return 'Medium'
        elif score == 1:
            return 'Low'
        
    def get_strength(e1, e2):
        score = e1 + e2

        # High(3) + High(3) = 6 -> High
        # High(3) + Medium(2) = 5 -> High
        # High(3) + Low(1) = 4 -> Medium
        # Medium(2) + Medium(2) = 4 -> Medium
        # Medium(2) + Low(1) = 3 -> Low
        # Low(1) + Low(1) = 2 -> Low

        if score >= 5: return 'High'
        if score >= 4: return 'Medium'
        return 'Low'
    
    def get_strength(e1, e2):
        score = e1 + e2  #若两个都是 High，则为 High，都是 Low 则为 Low，High-Low 则为 Medium
        if score == 6:
            return 'High'
        elif score == 2:
            return 'Low'
        elif score == 4:
            return 'Medium'

    for go_id, items in groups.items():
        n = len(items)
        if n < 2: continue

        # items 是 [('id1', 'High'), ('id2', 'Low')]
        # 组合两两配对
        for (p1_id, p1_ev), (p2_id, p2_ev) in combinations(items, 2):

            p1_class = evidence_classify.get(p1_ev, 'Low')
            p1_rank = evidence_rank.get(p1_class, 1)
            p2_class = evidence_classify.get(p2_ev, 'Low')
            p2_rank = evidence_rank.get(p2_class, 1)
            strength = get_strength(p1_rank, p2_rank)

            # 排序 ID，保证一致性
            if p1_id < p2_id:
                pair_key = (p1_id, p2_id)
                # 记录：ID1, ID2, GO_ID, Ev1, Ev2, Strength
                pair_data.append([p1_id, p2_id, go_id, p1_ev, p2_ev, p1_class, p2_class, strength])
            else:
                pair_key = (p2_id, p1_id)
                pair_data.append([p2_id, p1_id, go_id, p2_ev, p1_ev, p2_class, p1_class,  strength])

            pair_counter[pair_key] += 1

    print(f"初步生成了 {len(pair_data)} 行关系数据。")
    df_go_pairs = pd.DataFrame(pair_data, columns=[
        'Protein1', 'Protein2', 'GO_ID', 'Evidence1', 'Evidence2',
        'Ev1_Level', 'Ev2_Level', 'GO_Strength'
    ])


    df_go_pairs['temp_key'] = list(zip(df_go_pairs['Protein1'], df_go_pairs['Protein2']))
    df_go_pairs['Support_by_GOTerm'] = df_go_pairs['temp_key'].map(pair_counter)
    df_go_pairs = df_go_pairs.drop(columns=['temp_key'])

    print("go_pairs 的前 5 行预览：")
    print(df_go_pairs.head())

    print(f"\n已生成 go_pairs 文件")



    #生成 pairs 后，和 string 数据库里映射上
    df_links['p1_sorted'] = df_links[['protein1', 'protein2']].min(axis=1)
    df_links['p2_sorted'] = df_links[['protein1', 'protein2']].max(axis=1)

    # string-links 里 A-B 和 B-A 的蛋白对，只保留一对，且遵守左小右大
    df_links = df_links.drop_duplicates(subset=['p1_sorted', 'p2_sorted'])
    def classify_string_score(score):
        if score >= 700:
            return 'High'
        elif score >= 400:
            return 'Medium'
        elif score > 0:
            return 'Low'
        else:
            return 'None'

    df_merged = pd.merge(
        df_go_pairs,
        df_links,
        left_on=['Protein1', 'Protein2'],
        right_on=['p1_sorted', 'p2_sorted'],
        how='outer'
    )
    df_merged['Protein1'] = df_merged['Protein1'].fillna(df_merged['p1_sorted'])
    df_merged['Protein2'] = df_merged['Protein2'].fillna(df_merged['p2_sorted'])
    df_merged['original_combined_score'] = df_merged['original_combined_score'].fillna(0).astype(int)
    df_merged['combined_score'] = df_merged['combined_score'].fillna(0).astype(int)
    df_merged['GO_Strength'] = df_merged['GO_Strength'].fillna('None')
    df_merged['Support_by_GOTerm'] = df_merged['Support_by_GOTerm'].fillna(0).astype(int)
    df_merged['STRING_Strength'] = df_merged['combined_score'].apply(classify_string_score)
    df_merged['GO_ID'] = df_merged['GO_ID'].fillna('None')

    #df_merged = df_merged[['Protein1', 'Protein2', 'GO_ID', 'Evidence1', 'Evidence2',
    #                       'Ev1_Level', 'Ev2_Level', 'combined_score', 'GO_Strength', 
    #                       'STRING_Strength', 'Support_by_GOTerm']]
    




    has_score = df_merged[df_merged['combined_score'] > 0]
    no_score = df_merged[df_merged['combined_score'] == 0]

    print(f"总共有 {len(df_merged)} 对 GO 共注释蛋白对。")
    print(f"1. 有 STRING 关联的 (Category 1): {len(has_score)} 对 ({(len(has_score) / len(df_merged) * 100):.1f}%)")
    print(f"2. 无 STRING 关联的 (Category 2): {len(no_score)} 对 ({(len(no_score) / len(df_merged) * 100):.1f}%)")
    return df_merged

#从 NCBI 数据里添加基因名与方向信息
def add_gene_names(df_merged, df_NCBI, taxon_id):
    len1 = len(df_merged)
    gene_rows = df_NCBI[df_NCBI['# feature'] == 'gene']
    
    # 创建两个字典：一个用原始 locus_tag（可能有下划线），一个用去掉下划线的 locus_tag
    locus_to_symbol = dict(zip(gene_rows['locus_tag'], gene_rows['symbol']))
    locus_to_strand = dict(zip(gene_rows['locus_tag'], gene_rows['strand']))
    
    # 创建去掉下划线的版本，用于匹配 STRING 格式
    no_us_to_symbol = dict(zip(gene_rows['locus_tag'].str.replace('_', ''), gene_rows['symbol']))
    no_us_to_strand = dict(zip(gene_rows['locus_tag'].str.replace('_', ''), gene_rows['strand']))
    
    # 对于特定物种（如 E. lenta taxon_id=479437），需要从 attributes 列提取 old_locus_tag
    old_locus_to_symbol = {}
    old_locus_to_strand = {}
    if str(taxon_id) in ('479437', '272563', '1140'):
        import re
        for _, row in gene_rows.iterrows():
            attrs = row.get('attributes', '')
            if pd.notna(attrs):
                match = re.search(r'old_locus_tag=([^;]+)', str(attrs))
                if match:
                    old_locus = match.group(1)
                    old_locus_to_symbol[old_locus] = row['symbol']
                    old_locus_to_strand[old_locus] = row['strand']
        print(f"  E. lenta: 从 attributes 提取了 {len(old_locus_to_symbol)} 个 old_locus_tag")

    def get_name(string_id):
        # 提取："224308.BSU00010" -> "BSU00010"
        locus_tag = str(string_id).split('.')[1]
        
        # 对于 E. lenta，优先使用 old_locus_tag 匹配
        if str(taxon_id) in ('479437', '272563', '1140') and locus_tag in old_locus_to_symbol:
            return old_locus_to_symbol[locus_tag]
        
        # 先尝试用无下划线版本查找，再用原始版本查找
        if locus_tag in no_us_to_symbol:
            return no_us_to_symbol[locus_tag]
        return locus_to_symbol.get(locus_tag, locus_tag)
    
    def get_strand(string_id):
        # 提取： "224308.BSU00010" -> "BSU00010"
        locus_tag = str(string_id).split('.')[1]
        
        # 对于 E. lenta，优先使用 old_locus_tag 匹配
        if str(taxon_id) in ('479437', '272563', '1140') and locus_tag in old_locus_to_strand:
            return old_locus_to_strand[locus_tag]
        
        # 先尝试用无下划线版本查找，再用原始版本查找
        if locus_tag in no_us_to_strand:
            return no_us_to_strand[locus_tag]
        return locus_to_strand.get(locus_tag, None)


    df_merged['Gene1'] = df_merged['Protein1'].apply(get_name)
    df_merged['Gene2'] = df_merged['Protein2'].apply(get_name)
    df_merged['Strand1'] = df_merged['Protein1'].apply(get_strand)
    df_merged['Strand2'] = df_merged['Protein2'].apply(get_strand)
    df_merged = df_merged.dropna(subset = ['Strand2'])
    len2 = len(df_merged)
    print(f"  -> 原来有{len1}对蛋白对，添加基因名和方向信息后，剩余 {len2} 对蛋白对")
    return df_merged
