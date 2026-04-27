import pandas as pd
import plotly.express as px
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.decomposition import TruncatedSVD
import species_config
import os
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network



def plot_interactive_bipartite(match_table, top_n=100):
    # === 1. 数据检查 (这是最容易出错的地方) ===

    df = match_table.head(5000)
    print(f"数据读取成功，共 {len(df)} 行。")

    # 筛选 Top N 功能
    top_go_ids = df['GO_ID'].value_counts().head(top_n).index.tolist()
    subset = df[df['GO_ID'].isin(top_go_ids)]

    print(f"筛选 Top {top_n} GO Term 后，剩余 {len(subset)} 行数据。")
    if len(subset) == 0:
        print("警告：筛选后数据为空！请检查你的 GO_ID 列是否正确。")
        return

    # === 2. 初始化网络 (关键修改) ===
    # cdn_resources='in_line': 这一步非常重要！
    # 它可以把所有 JS 脚本打包进 HTML，防止因为网络问题加载不出 PyVis 库导致白屏
    net = Network(height="750px", width="100%", bgcolor="#222222", font_color="white", cdn_resources='in_line')

    # === 3. 构建网络 ===
    print("正在构建网络结构...")

    # 添加 GO 节点 (功能中心)
    for go in top_go_ids:
        net.add_node(go, label=go, title=f"Function: {go}", color='#FF9F40', size=25, shape='square')

    # 添加 Gene 节点 (去重)
    # 注意：subset 中同一个 gene 可能出现多次（对应不同 GO），但节点只能加一次
    unique_genes = subset[['Gene_Symbol', 'UniProt_ID', 'Protein_Name']].drop_duplicates(subset=['Gene_Symbol'])

    for _, row in unique_genes.iterrows():
        gene = row['Gene_Symbol']
        # 构造悬停信息
        hover_info = (
            f"<b>Gene:</b> {gene}<br>"
            f"<b>UniProt:</b> {row['UniProt_ID']}<br>"
            f"<b>Name:</b> {row['Protein_Name']}"
        )
        net.add_node(gene, label=gene, title=hover_info, color='#4BC0C0', size=10)

    # 添加连线
    for _, row in subset.iterrows():
        net.add_edge(row['Gene_Symbol'], row['GO_ID'], color='#555555')

    # === 4. 物理引擎稳定 (解决“炸飞”问题) ===
    # 强制使用 BarnesHut 算法，并将引力常数设为较小的值，防止初始速度过快
    net.barnes_hut(gravity=-8000, central_gravity=0.3, spring_length=95, spring_strength=0.04, damping=0.8)
    output_file = "interactive_bipartite.html"
    output_path = os.path.join(species_config.output_figure_dir, output_file)
    print(f"正在保存文件到: {output_path}")

    # 1. 先生成 HTML 字符串
    # 注意：generate_html() 需要传入文件名作为参数，以便它处理路径引用，但这不会真的写文件
    html_content = net.generate_html(output_path)

    # 2. 手动用 utf-8 格式打开并写入
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)



from matplotlib_venn import venn2
import re
#计算protein覆盖情况
def calculate_coverage_metrics(df_NCBI, df_match, df_string_aliases):
    def normalize(tag):
        return tag.replace('_', '') if tag else None
    # 计算基因集

    genome_genes = set(df_NCBI[
        (df_NCBI['# feature'] == 'gene') & 
        (df_NCBI['class'] == 'protein_coding')
    ]['locus_tag'])
    print(f"  -> 全基因组共发现 {len(genome_genes)} 个蛋白编码基因")
    
    # 取 GO 基因集
    # 取交集: 只有在 NCBI 全集里的才算有效覆盖
    #go_genes = set(df_match['String_ID'].dropna().apply(lambda x: str(x).split('.')[1] if '.' in str(x) else None))
    #go_genes = go_genes.intersection(genome_genes)
    genome_genes_norm = set(normalize(t) for t in genome_genes if t)

    go_genes_raw = set(df_match['String_ID'].dropna()
    .apply(lambda x: str(x).split('.')[1] if '.' in str(x) else None))
    
    go_genes = set(t for t in go_genes_raw if normalize(t) in genome_genes_norm)
    print(f"  -> GO 数据库有效覆盖 {len(go_genes)} 个基因")



    # 取 STRING 基因集
    #string_genes = set(df_string_aliases['#string_protein_id'].apply(lambda x: str(x).split('.')[1] if '.' in str(x) else None))
    # 同样取交集
    #string_genes = string_genes.intersection(genome_genes)
    string_genes_raw = set(df_string_aliases['#string_protein_id']
        .apply(lambda x: str(x).split('.')[1] if '.' in str(x) else None))
    string_genes = set(t for t in string_genes_raw if normalize(t) in genome_genes_norm)
    print(f"  -> STRING 数据库有效覆盖 {len(string_genes)} 个基因")

    return genome_genes, go_genes, string_genes


#计算pairs覆盖情况
def calculate_pair_metrics(df_NCBI, df_pairs):
    """
    计算三个集合的大小：
    1. Total Pairs (理论所有C2^2的最大值)
    2. GO Supported Pairs (pairs表中被GO支持的)
    3. STRING Supported Pairs (pairs表中被STRING支持的)
    """
    metrics = {}
    # 剔除df_pairs里strand1和strand2不同的数据
    df_pairs = df_pairs[df_pairs['Strand1'] == df_pairs['Strand2']]
    
    # --- A. 计算 NCBI 全组合 (Universe) ---
    print(f"正在读取 NCBI 全基因组: {df_NCBI} ...")
    # 筛选蛋白编码基因
    valid_genes = df_NCBI[
        (df_NCBI['# feature'] == 'gene') & 
        (df_NCBI['class'] == 'protein_coding')
    ]
    N = len(valid_genes)
    total_genome_pairs = (N * (N - 1)) // 2
    print(f"  -> 基因数 N={N}, 理论全组合数 = {total_genome_pairs:,}")
 
    # --- B. 读取 Pairs 表 ---
    # 构造唯一 Pair ID (排序去重)
   
    df_pairs['Pair_ID'] = df_pairs.apply(lambda row: '_'.join(sorted([str(row['Protein1']), str(row['Protein2'])])), axis=1)
    # --- C. 提取 GO 集合 ---
    # 逻辑: GO_Strength 不为空且不为0
    # 先处理优先级 (High > Medium > Low)，确保每对只算一次
    strength_map = {'High': 3, 'Medium': 2, 'Low': 1}
    df_pairs['go_score'] = df_pairs['GO_Strength'].map(strength_map).fillna(0)

    # 提取支持 GO 的唯一 Pair ID
    go_pairs = set(df_pairs[df_pairs['go_score'] > 0]['Pair_ID'])
    print(df_pairs)
    print(f"  -> GO 支持的对数: {len(go_pairs):,}")

    # --- D. 提取 STRING 集合 ---
    if 'combined_score' in df_pairs.columns:
        string_pairs = set(df_pairs[df_pairs['combined_score'] > 0]['Pair_ID'])
    else:
        # 兜底: 如果没有 Score 列，看 ID 列
        string_pairs = set(df_pairs[df_pairs['Protein1'].notna()]['Pair_ID']) 
        
    print(f"  -> STRING 支持的对数: {len(string_pairs):,}")
    
    return total_genome_pairs, go_pairs, string_pairs



def calculate_high_confidence_pairs(df_NCBI, df_pairs):
    """
    计算三个集合的大小：
    1. Total Pairs (理论所有C2^2的最大值)
    2. GO Supported Pairs (pairs表中被GO支持的)
    3. STRING Supported Pairs (pairs表中被STRING支持的)
    """
    # 剔除df_pairs里strand1和strand2不同的数据
    df_pairs = df_pairs[df_pairs['Strand1'] == df_pairs['Strand2']]

    metrics = {}
    
    # --- A. 计算 NCBI 全组合 (Universe) ---
    # 筛选蛋白编码基因
    valid_genes = df_NCBI[
        (df_NCBI['# feature'] == 'gene') & 
        (df_NCBI['class'] == 'protein_coding')
    ]
    N = len(valid_genes)
    total_genome_pairs = (N * (N - 1)) // 2
    print(f"  -> 基因数 N={N}, 理论全组合数 = {total_genome_pairs:,}")
 
    # --- B. 读取 Pairs 表 ---
    # 构造唯一 Pair ID (排序去重)
   
    df_pairs['Pair_ID'] = df_pairs.apply(lambda row: '_'.join(sorted([str(row['Protein1']), str(row['Protein2'])])), axis=1)
    # --- C. 提取 GO 集合 ---
    # 逻辑: GO_Strength 不为空且不为0
    # 先处理优先级 (High > Medium > Low)，确保每对只算一次
    strength_map = {'High': 3, 'Medium': 2, 'Low': 1}
    df_pairs['go_score'] = df_pairs['GO_Strength'].map(strength_map).fillna(0)

    #提取高置信度的蛋白支持对
    df_pairs = df_pairs[(df_pairs['go_score'] == 3) | (df_pairs['STRING_Strength'] == 'High')]

    # 提取支持 GO 的唯一 Pair ID
    go_pairs = set(df_pairs[df_pairs['go_score'] > 0]['Pair_ID'])

    # --- D. 提取 STRING 集合 ---
    if 'combined_score' in df_pairs.columns:
        string_pairs = set(df_pairs[df_pairs['combined_score'] > 0]['Pair_ID'])
    else:
        # 兜底: 如果没有 Score 列，看 ID 列
        string_pairs = set(df_pairs[df_pairs['Protein1'].notna()]['Pair_ID']) 
        
    
    return total_genome_pairs, go_pairs, string_pairs

#画韦恩图
def plot_venn(total_set, go_set, string_set, species_config):
    # 准备画布
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # --- 图 1: 覆盖率柱状图 ---
    total_count = len(total_set)
    ratios = [len(go_set)/total_count, len(string_set)/total_count]
    labels = ['GO', 'STRING']
    colors = ['#4e79a7', '#f28e2b']
    
    bars = axes[0].bar(labels, ratios, color=colors, alpha=0.8, width=0.5)
    
    # 装饰 - 增加 Y 轴上限留出空间
    axes[0].axhline(1.0, color='gray', linestyle='--', label='NCBI Total Genome (100%)')
    axes[0].set_ylim(0, 1.25)  # 增加上限
    axes[0].set_ylabel('Coverage Ratio', fontsize=12)
    axes[0].set_title(f'Annotation Coverage (Total Genes: {total_count})', fontsize=14, pad=15)
    axes[0].legend(loc='lower right')  # 图例移到右下角避免遮挡
    
    # 标数值 - 增加间距
    for bar in bars:
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2., height + 0.05,
                     f'{height:.1%}\n({int(height*total_count)})',
                     ha='center', va='bottom', fontsize=12, fontweight='bold')

    # --- 图 2: 韦恩图 (GO vs STRING) ---

    venn = venn2([go_set, string_set], 
                 set_labels=('GO Annotated', 'STRING Annotated'), 
                 set_colors=('#4e79a7', '#f28e2b'),
                 alpha=0.6,
                 ax=axes[1])
    
    axes[1].set_title('Overlap of Annotated Proteins', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], "proteins overlap.png"), dpi=300)
    plt.close()

def plot_pairs_venn(total, go_set, string_set, species_config):

    # 准备画布: 左边是覆盖率(Bar)，右边是重叠(Venn)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # --- 图 1: 基因组空间覆盖率 (Bar Chart) ---
    # 用 Log Scale 或者断轴，或者直接标数字
    
    counts = [len(go_set), len(string_set)]
    ratios = [c / total for c in counts]
    labels = ['GO Supported', 'STRING Supported']
    colors = ['#4e79a7', '#f28e2b']
    
    bars = axes[0].bar(labels, ratios, color=colors, alpha=0.8, width=0.5)
    
    # 装饰 - 增加 Y 轴上限留出空间
    max_ratio = max(ratios)
    axes[0].set_ylim(0, max_ratio * 1.35)  # 增加 35% 空间
    axes[0].set_ylabel('Coverage Ratio (vs Total Genome Pairs)', fontsize=12)
    axes[0].set_title(f'Sparsity of Functional Network\n(Total Possible Pairs: {total:,})', fontsize=14, pad=15)
    
    # 标注极小的数值 - 增加间距
    for bar in bars:
        height = bar.get_height()
        count = int(height * total)
        # 在柱子上方写百分比和绝对数量
        axes[0].text(bar.get_x() + bar.get_width()/2., height + max_ratio * 0.02,
                     f'{height:.4%}\n(n={count:,})',
                     ha='center', va='bottom', fontsize=12, fontweight='bold')

    # --- 图 2: 验证集韦恩图 (Venn Diagram) ---
    # 这里不画 Total，只画 GO 和 STRING 的交集
    # 因为 Total 是背景板，画进韦恩图会让另外两个圈缩成像素点
    
    venn2([go_set, string_set],
          set_labels=('GO Supported', 'STRING Supported'),
          set_colors=('#4e79a7', '#f28e2b'),
          alpha=0.6,
          ax=axes[1])
    
    axes[1].set_title('Overlap of Validated Pairs', fontsize=14, pad=15)
    
    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(species_config['output_figure_dir'], "pairs overlap.png"), dpi=300, bbox_inches='tight')
    plt.close()

def plot_high_confidence_pairs_venn(total, go_set, string_set, high_go_set, high_string_set, species_config):

    # 准备画布: 左边是覆盖率(Bar)，右边是重叠(Venn)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # --- 图 1: 基因组空间覆盖率 (Bar Chart) ---
    # 注意: 因为 Total 巨大 (800万)，而 GO/STRING 可能只有几千
    # 普通柱状图根本看不见，我们需要用 Log Scale 或者断轴，或者直接标数字
    
    # High Confidence Counts
    high_counts = [len(high_go_set), len(high_string_set)]
    high_ratios = [c / total for c in high_counts]
    
    # All Pairs Counts
    all_counts = [len(go_set), len(string_set)]
    all_ratios = [c / total for c in all_counts]
    
    labels = ['GO Supported', 'STRING Supported']
    
    # Plot All Pairs (Background, Lighter)
    bars_all = axes[0].bar(labels, all_ratios, color=['#4e79a7', '#f28e2b'], alpha=0.6, width=0.5, label='All Pairs')
    
    # Plot High Confidence Pairs (Foreground, Darker)
    bars_high = axes[0].bar(labels, high_ratios, color=["#24059e", "#e10404"], alpha=1.0, width=0.5, label='High Confidence')
    
    # 装饰 - 增加 Y 轴上限留出空间
    max_all_ratio = max(all_ratios)
    axes[0].set_ylim(0, max_all_ratio * 1.35)  # 增加 35% 空间
    axes[0].set_ylabel('Coverage Ratio (vs Total Genome Pairs)', fontsize=12)
    axes[0].set_title(f'Sparsity of Functional Network\n(Total Possible Pairs: {total:,})', fontsize=14, pad=15)
    axes[0].legend(loc='upper left')
    
    # 标注极小的数值 - 增加间距
    for i, (bar_all, bar_high) in enumerate(zip(bars_all, bars_high)):
        height_all = bar_all.get_height()
        
        count_all = all_counts[i]
        count_high = high_counts[i]
        ratio_in_group = count_high / count_all if count_all > 0 else 0
        
        # 在柱子上方写百分比和绝对数量
        axes[0].text(bar_all.get_x() + bar_all.get_width()/2., height_all + max_all_ratio * 0.02,
                     f'All: {count_all:,}\nHigh: {count_high:,}\n({ratio_in_group:.1%})',
                     ha='center', va='bottom', fontsize=12, fontweight='bold')

    # --- 图 2: 验证集韦恩图 (Venn Diagram) ---
    # 这里不画 Total，只画 GO 和 STRING 的交集
    # 因为 Total 是背景板，画进韦恩图会让另外两个圈缩成像素点
    
    venn2([high_go_set, high_string_set],
          set_labels=('GO Supported', 'STRING Supported'),
          set_colors=("#24059e", "#e10404"),
          alpha=0.6,
          ax=axes[1])
    
    axes[1].set_title('Overlap of High Confidence Validated Pairs', fontsize=14, pad=15)

    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(species_config['output_figure_dir'], "high confidence pairs overlap.png"), dpi=300, bbox_inches='tight')
    plt.close()