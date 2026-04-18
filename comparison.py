import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import config
from matplotlib_venn import venn2
from scipy.stats import spearmanr
from scipy.stats import kruskal
import scikit_posthocs as sp
import os


#把 STRING 的分数和保守性做比对
reference_count = 1308  # 代表 1308 个物种的保守矩阵总分是 1308

def plot_spacedust_pairs_venn(total, go_set, string_set, species_config):

    # 准备画布: 左边是覆盖率(Bar)，右边是重叠(Venn)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # --- 图 1: 基因组空间覆盖率 (Bar Chart) ---
    # 注意: 因为 Total 巨大 (800万)，而 GO/STRING 可能只有几千
    # 普通柱状图根本看不见，我们需要用 Log Scale 或者断轴，或者直接标数字
    
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

    # --- 图 2: 韦恩图 (带全集圆) ---
    from matplotlib.patches import Circle
    
    ax = axes[1]
    ax.set_aspect('equal')
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2, 2)
    ax.axis('off')
    
    # 计算集合统计
    go_keys = set(go_set['pair_key'])
    string_keys = set(string_set['pair_key'])
    go_total = len(go_keys)
    string_total = len(string_keys)
    go_only = len(go_keys - string_keys)
    string_only = len(string_keys - go_keys)
    intersection = len(go_keys & string_keys)
    universe_only = total - len(go_keys | string_keys)
    
    # 画全集圆 (灰色背景)
    universe_circle = Circle((0, 0), 1.5, fill=True, color='lightgray', alpha=0.5)
    ax.add_patch(universe_circle)
    ax.text(0, 1.9, f'Universe (n={total:,})', ha='center', va='bottom', fontsize=11,
            fontweight='bold', color='gray')
    
    # 根据数据量比例设置圆的大小 (面积与数量成正比)
    max_radius = 0.9
    go_radius = max_radius * (go_total / string_total) ** 0.5  # GO 圆半径按比例缩小
    string_radius = max_radius
    
    # 调整圆心位置，确保有适当重叠
    go_center_x = -0.5
    string_center_x = 0.5
    
    # 画 GO 圆 (蓝色) - 较小的圆
    go_circle = Circle((go_center_x, 0), go_radius, fill=True, color='#4e79a7', alpha=0.6)
    ax.add_patch(go_circle)
    
    # 画 STRING 圆 (橙色) - 较大的圆
    string_circle = Circle((string_center_x, 0), string_radius, fill=True, color='#f28e2b', alpha=0.6)
    ax.add_patch(string_circle)
    
    # 标注各区域数量 - 根据圆大小调整位置
    # GO only: 在 GO 圆左侧非重叠区域
    ax.text(go_center_x - go_radius * 0.5, 0, f'{go_only:,}', ha='center', va='center', fontsize=12,
            fontweight='bold', color='white')
    # STRING only: 在 STRING 圆右侧非重叠区域
    ax.text(string_center_x + string_radius * 0.4, 0, f'{string_only:,}', ha='center', va='center', fontsize=14,
            fontweight='bold', color='white')
    # Intersection: 在重叠区域中心偏上，避免与圆边界重合
    overlap_center_x = go_center_x + go_radius - (go_radius + string_radius - (string_center_x - go_center_x)) / 2
    ax.text(overlap_center_x, 0.15, f'{intersection:,}', ha='center', va='center', fontsize=12,
            fontweight='bold', color='white')
    
    # 标注全集独有部分
    ax.text(0, -1.3, f'Not in GO/STRING: {universe_only:,}', ha='center', va='center',
            fontsize=10, style='italic', color='gray')
    
    # 图例 - 显示实际数量
    ax.text(-1.8, -1.7, f'GO (n={go_total:,})', color='#4e79a7', fontsize=10, fontweight='bold')
    ax.text(0.8, -1.7, f'STRING (n={string_total:,})', color='#f28e2b', fontsize=10, fontweight='bold')
    
    axes[1].set_title('Overlap of Validated Pairs', fontsize=14, pad=15)
    
    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(species_config['output_figure_dir'], "spacedust pairs overlap.png"), dpi=300, bbox_inches='tight')
    plt.close()


def compare_pairs(df_matrix, df_lookup, df_pairs, species_config):

    #CORE Algorithm：把 1308 个物种的分数加起来，得到每一对基因的 Total Score
    # matrix_df.sum(axis=0) 会对每一列求和，返回一个 Series：索引是列名 (0,1,2...)，值是总分
    spacedust_scores = df_matrix.sum(axis=0).reset_index()
    spacedust_scores.columns = ['matrix_idx', 'spacedust_score']
    spacedust_scores['spacedust_score_proportion'] = spacedust_scores['spacedust_score']
    spacedust_scores['matrix_idx'] = spacedust_scores['matrix_idx'].astype(int)

    print("正在翻译矩阵索引...")

    # 构建字典：{矩阵索引 (idx): tag(511145.b0001)}
    print(df_lookup.head())
    idx_to_id = dict(zip(df_lookup['idx'], df_lookup['STRING_ID']))

    # 遍历每一个矩阵列（基因对），提取出真正的蛋白对名字
    protein1_list = []
    protein2_list = []

    for idx in spacedust_scores['matrix_idx']:
        # 根据我们的物理映射法则：矩阵列 X，对应 lookup 里的索引 X 和 X+1
        p1 = idx_to_id.get(idx, np.nan)
        p2 = idx_to_id.get(idx + 1, np.nan)
        protein1_list.append(p1)
        protein2_list.append(p2)

    spacedust_scores['protein1'] = protein1_list
    spacedust_scores['protein2'] = protein2_list
    spacedust_scores = spacedust_scores.dropna(subset=['protein1', 'protein2'])

    spacedust_scores['pair_key'] = spacedust_scores.apply(
        lambda row: "_".join(sorted([str(row['protein1']), str(row['protein2'])])), axis=1
    )


    print("正在加载并合并 STRING 互作数据...")
    STRING_P1_COL = "Protein1"
    STRING_P2_COL = "Protein2"
    string_pairs = df_pairs
    string_pairs['pair_key'] = string_pairs.apply(
        lambda row: "_".join(sorted([str(row[STRING_P1_COL]), str(row[STRING_P2_COL])])), axis=1
    )

    string_pairs = string_pairs.drop_duplicates(subset=['pair_key'], keep='first')

    # 定义要分析的 STRING 分数列
    string_score_columns = ['original_combined_score', 'combined_score', 'fusion_score', 'experimental_score',
                             'database_score', 'textmining_score', 'neighborhood_score', 'cooccurence_score']
    
    # 检查哪些列存在于数据中
    available_score_cols = [col for col in string_score_columns if col in string_pairs.columns]
    print(f"可用的 STRING 分数列：{available_score_cols}")
    
    # 利用 pair_key 将进化矩阵分数与 STRING 分数拼在一起！
    # 使用 outer join 保留所有数据，包括 Spacedust 独有和 STRING 独有的蛋白对
    merged_df = pd.merge(
        spacedust_scores,
        string_pairs[['pair_key'] + available_score_cols],
        on='pair_key',
        how='outer'
    )

    # 在 fillna 之前先标记数据来源，用于区分交集、仅 Spacedust、仅 STRING
    merged_df['has_spacedust'] = merged_df['spacedust_score'].notna()
    merged_df['has_string'] = merged_df['combined_score'].notna() if 'combined_score' in merged_df.columns else False
    
    # 将所有 STRING 分数从 0-1000 缩放到 0-1 (概率)，缺失值设为 0
    for col in available_score_cols:
        merged_df[col] = (merged_df[col] / 1000.0).fillna(0)
    # spacedust_score 缺失值也设为 0
    merged_df['spacedust_score'] = merged_df['spacedust_score'].fillna(0)

    # 数据统计
    total_spacedust_pairs = len(spacedust_scores)
    total_string_pairs = len(string_pairs)
    # 交集：同时有 spacedust 和 STRING 数据
    intersection_pairs = merged_df[(merged_df['has_spacedust']) & (merged_df['has_string'])]
    # 仅 Spacedust：有 spacedust 但没有 STRING
    only_spacedust = merged_df[(merged_df['has_spacedust']) & (~merged_df['has_string'])]
    # 仅 STRING：有 STRING 但没有 spacedust
    only_string = merged_df[(~merged_df['has_spacedust']) & (merged_df['has_string'])]
    # 交集中 STRING 无数据的样本量（即 combined_score 为 0 的数量）
    intersection_no_string_data = intersection_pairs[intersection_pairs['combined_score'] == 0] if 'combined_score' in intersection_pairs.columns else pd.DataFrame()

    # 画图时使用交集数据（两者都有的）
    plot_df = intersection_pairs[intersection_pairs['original_combined_score'] > 0].copy()
    print(f"散点图将使用 {len(plot_df):,} 对交集数据进行绘制。")

    # ==========================================
    # 第四步：绘制生物学发现（散点图与统计学检验）- 多个 STRING 分数版本
    # ==========================================
    # 分数列名称映射（用于图表标题）
    score_labels = {
        'original_combined_score': 'Original Combined Score',
        'combined_score': 'Revised Combined Score',
        'fusion_score': 'Fusion Score',
        'experimental_score': 'Experimental Score',
        'database_score': 'Database Score',
        'textmining_score': 'Textmining Score',
        'coexpression_score': 'Coexpression Score',
        'cooccurence_score' : 'Cooccurence Score',
        'neighborhood_score' : 'Neighborhood Score'
    }
    
    for score_col in available_score_cols:
        # 计算斯皮尔曼相关系数
        valid_data = plot_df[(plot_df['spacedust_score_proportion'].notna()) & (plot_df[score_col].notna())]
        if len(valid_data) < 2:
            print(f"跳过 {score_col}：有效数据不足")
            continue
            
        corr, pval = spearmanr(valid_data['spacedust_score_proportion'], valid_data[score_col])
        print(f"[{score_col}] Spearman 相关系数：{corr:.3f}, p-value: {pval:.2e}")

        # 画散点图
        plt.figure(figsize=(10, 7))
        sns.regplot(x='spacedust_score_proportion', y=score_col, data=valid_data,
                    scatter_kws={'alpha':0.4, 's':20, 'color':'#2c7bb6'},
                    line_kws={'color':'#d7191c', 'linewidth':2})

        # 装饰图表
        score_label = score_labels.get(score_col, score_col.replace('_', ' ').title())
        plt.title(f'Evolutionary Conservation vs. STRING {score_label}\n(Spearman r={corr:.3f}, p={pval:.2e})', fontsize=14)
        # 添加样本量标注
        plt.text(0.98, 0.02, f'n = {len(valid_data):,}', transform=plt.gca().transAxes,
                 fontsize=12, ha='right', va='bottom', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.xlabel('Spacedust Conservation Score (Evolutionary Evidence)', fontsize=12)
        plt.ylabel(f'STRING {score_label} (0-1)', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.5)
        sns.despine()

        # 保存文件
        filename = f'spacedust_vs_string_{score_col.replace("_score", "")}_scatter.png'
        plt.savefig(os.path.join(species_config['output_figure_dir'], filename), dpi=300, bbox_inches='tight')
        plt.close()


    # ==========================================
    # 第五步：可视化 Spacedust 独有数据的分布
    # ==========================================
    # 把交集以外的数据分离出来，分析 Spacedust 独有的数据分布
    print(f"\n正在分析 {len(only_spacedust):,} 对 Spacedust 独有（STRING 无）的蛋白对...")
    
    
    # 画 Spacedust 独有数据的分数分布直方图
    plt.figure(figsize=(10, 6))
    plt.hist(only_spacedust['spacedust_score_proportion'], bins=20, color='#d7191c', alpha=0.7, edgecolor='black')
    plt.title(f'Spacedust-Only Pairs Conservation Score Distribution\n(n = {len(only_spacedust):,})', fontsize=14)
    plt.xlabel('Spacedust Conservation Score', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'spacedust_only_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # 画 原始的保守分数分布直方图
    plt.figure(figsize=(10, 6))
    plt.hist(spacedust_scores['spacedust_score_proportion'], bins=20, color='#d7191c', alpha=0.7, edgecolor='black')
    plt.title(f'Original Pairs Conservation Score Distribution\n(n = {len(spacedust_scores):,})', fontsize=14)
    plt.xlabel('Spacedust Conservation Score', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'conserved_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()

    
    # 画箱线图对比交集和 Spacedust 独有的分数分布
    plt.figure(figsize=(8, 6))
    compare_data = [
        plot_df['spacedust_score_proportion'].dropna(),
        only_spacedust['spacedust_score_proportion'].dropna()
    ]
    bp = plt.boxplot(compare_data, labels=['Intersection\n(with STRING)', 'Spacedust Only\n(no STRING)'],
                     patch_artist=True, widths=0.5)
    colors = ['#2c7bb6', '#d7191c']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    plt.title('Spacedust Score Distribution: Intersection vs. Spacedust-Only', fontsize=14)
    plt.ylabel('Spacedust Conservation Score', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5, axis='y')
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'intersection_vs_spacedust_only_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    

    go_pairs = df_pairs



    GO_P1_COL = "Protein1"
    GO_P2_COL = "Protein2"
    GO_STRENGTH_COL = "GO_Strength"
    # 先用 'None' 填充 NaN 值，再转换为字符串，确保所有空值都正确转换为 'None'
    go_pairs[GO_STRENGTH_COL] = go_pairs[GO_STRENGTH_COL].fillna('None').astype(str)
    level_to_num = {'None': 0, 'Low': 1, 'Medium': 2, 'High': 3}
    num_to_level = {0: 'None', 1: 'Low', 2: 'Medium', 3: 'High'}
    #pair_key
    go_pairs['pair_key'] = go_pairs.apply(
        lambda row: "_".join(sorted([str(row[GO_P1_COL]), str(row[GO_P2_COL])])), axis=1
        )
    go_pairs[GO_STRENGTH_COL] = go_pairs[GO_STRENGTH_COL].replace({'nan': 'None', '': 'None'})
    go_pairs['Strength_Num'] = go_pairs[GO_STRENGTH_COL].map(level_to_num)


    go_unique_df = go_pairs.groupby('pair_key')['Strength_Num'].max().reset_index()

    go_unique_df[GO_STRENGTH_COL] = go_unique_df['Strength_Num'].map(num_to_level)

    # 使用 left join 保留所有 Spacedust 数据，GO 无注释的会自动为 NaN，后续填充为'None'
    merged_go_df = pd.merge(spacedust_scores, go_unique_df[['pair_key', 'GO_Strength', 'Strength_Num']], on='pair_key', how='left')
    
    # 将 NaN 填充为'None'，表示无 GO 注释
    merged_go_df[GO_STRENGTH_COL] = merged_go_df[GO_STRENGTH_COL].fillna('None')
    merged_go_df['Strength_Num'] = merged_go_df['Strength_Num'].fillna(0).astype(int)

    print(f"共找到 {len(merged_go_df)} 对同时具有进化分数和 GO 注释的基因对。")

    # 交集中 GO 无数据的样本量（即 GO_Strength 为 None 的数量）
    intersection_no_go_data = merged_go_df[merged_go_df[GO_STRENGTH_COL] == 'None']
    intersection_go_data = merged_go_df[merged_go_df[GO_STRENGTH_COL] != 'None']


    # ==========================================
    # Kruskal-Wallis 检验 + Dunn's 事后检验
    # ==========================================
    groups_dict = dict(list(merged_go_df.groupby(GO_STRENGTH_COL)))
    groups = [(name, group['spacedust_score'].values) for name, group in groups_dict.items() if len(group) > 0]
    
    kw_pval = None
    dunn_results = None
    pairwise_comparisons = []
    
    # Kruskal-Wallis 检验
    if len(groups) >= 2:
        group_values = [g[1] for g in groups]
        kw_stat, kw_pval = kruskal(*group_values)
        print(f"\n多组差异非参数检验 (Kruskal-Wallis): χ² = {kw_stat:.4f}, p = {kw_pval:.2e}")
    else:
        print(f"\n警告：只有 {len(groups)} 个 GO_Strength 分组，无法进行 Kruskal-Wallis 检验")
        kw_stat, kw_pval = None, None
    
    # Dunn's 事后检验（仅对 High, Medium, Low 三组）
    target_groups = ['High', 'Medium', 'Low']
    valid_groups = [g for g in target_groups if g in groups_dict and len(groups_dict[g]) >= 2]
    
    if len(valid_groups) >= 2:
        print(f"\nDunn's 事后检验 (目标组：{valid_groups})")
        try:
            dunn_df = merged_go_df[merged_go_df[GO_STRENGTH_COL].isin(valid_groups)]
            # scikit-posthocs 参数：a=DataFrame, val_col=因变量列，group_col=分组列
            dunn_results = sp.posthoc_dunn(dunn_df, val_col='spacedust_score', group_col=GO_STRENGTH_COL, p_adjust='bonferroni')
            print("Dunn's test 结果 (Bonferroni 校正后 P 值):")
            print(dunn_results)
            
            # 提取两两比较结果
            for i, g1 in enumerate(valid_groups):
                for g2 in valid_groups[i+1:]:
                    if g1 in dunn_results.index and g2 in dunn_results.columns:
                        p_val = dunn_results.loc[g1, g2]
                        pairwise_comparisons.append((g1, g2, p_val))
                        sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                        print(f"  {g1} vs {g2}: p = {p_val:.4f} [{sig}]")
        except Exception as e:
            print(f"Dunn's test 执行失败：{e}")
            dunn_results = None
    else:
        print(f"\n警告：有效组数 ({len(valid_groups)}) < 2，无法进行 Dunn's 事后检验")
    
    # ==========================================
    # 绘制箱线图 + 显著性标记
    # ==========================================
    fig, ax = plt.subplots(figsize=(10, 7))

    sns.boxplot(x=GO_STRENGTH_COL, y='spacedust_score', data=merged_go_df, hue=GO_STRENGTH_COL,
                order=['High', 'Medium', 'Low', 'None'], palette='YlGnBu', width=0.5, showfliers=False, ax=ax)

    sns.stripplot(x=GO_STRENGTH_COL, y='spacedust_score', data=merged_go_df, order=['High', 'Medium', 'Low', 'None'],
                color='black', alpha=0.15, jitter=True, size=3, ax=ax)

    # 添加每组样本量标注
    group_counts = merged_go_df.groupby(GO_STRENGTH_COL).size()
    ylim = ax.get_ylim()
    y_offset = ylim[1] * 0.15
    for i, level in enumerate(['High', 'Medium', 'Low', 'None']):
        n = group_counts.get(level, 0)
        ax.text(i, ylim[1] + y_offset, f'n={n}', ha='center', fontsize=12, fontweight='bold')
    
    # 扩展 Y 轴上限以容纳显著性标记文字
    ax.set_ylim(ylim[0], ylim[1] + y_offset * 3.5)
    
    # ==========================================
    # 添加显著性文字标注（直接显示 P 值）
    # ==========================================
    if pairwise_comparisons:
        y_pos_base = ylim[1] + y_offset * 1.3  # 文字标注的基准 Y 位置
        y_step = y_offset * 0.35  # 每行的高度间隔
        
        for idx, (g1, g2, p_val) in enumerate(pairwise_comparisons):
            y_pos = y_pos_base + idx * y_step
            
            # 格式化 P 值
            if p_val < 0.001:
                p_str = 'p < 0.001'
            elif p_val < 0.01:
                p_str = 'p = {:.3f}'.format(p_val)
            elif p_val < 0.05:
                p_str = 'p = {:.3f}'.format(p_val)
            else:
                p_str = 'p = {:.3f} (ns)'.format(p_val)
            
            # 显示比较结果文字
            ax.text(-0.35, y_pos,
                   f'{g1} vs {g2}:  {p_str}',
                   ha='left', va='bottom', fontsize=10,
                   fontfamily='monospace',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.7))
    
    # 标题
    if kw_pval is not None:
        ax.set_title(f'Evolutionary Conservation vs. GO Evidence Strength\n(Kruskal-Wallis p={kw_pval:.2e})',
                    fontsize=14, pad=20)
    else:
        ax.set_title(f'Evolutionary Conservation vs. GO Evidence Strength\n(Kruskal-Wallis: N/A)',
                    fontsize=14, pad=20)
    
    ax.set_xlabel('Highest GO Evidence Strength (Max Pooled)', fontsize=12)
    ax.set_ylabel('Spacedust Conservation Score', fontsize=12)
    sns.despine()

    plt.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'spacedust_vs_go_strength.png'), dpi=300, bbox_inches='tight')
    plt.close()




    #merged_go_df['Spacedust_Tier'] = pd.cut(
    #    merged_go_df['spacedust_score'],
    #    bins=[float('-inf'), 3, 20, 40, 60, 80, 100, float('inf')],
    #    labels=['0-3', '3-20', '20-40', '40-60', '60-80', '80-100', '100+'],
    #    right=False
    #)

    merged_go_df['Spacedust_Tier'] = pd.cut(
        merged_go_df['spacedust_score'],
        bins=[0, 1e-5, 3, 10, 20, 40, 60, 80, 110, 140, 170, float('inf')],
        labels=['0', '1~2', '3~9', '10~19', '20~39', '40~59', '60~79', '80~109', '110~139', '140~169', '170+'],
        right=False
    )

    #统计每个梯队里，None, Low, Medium, High 的构成比例
    tier_go_counts = merged_go_df.groupby(['Spacedust_Tier', GO_STRENGTH_COL]).size().unstack(fill_value=0)

    # 把顺序排好，确保所有预期列都存在（缺失的列用 0 填充）
    expected_columns = ['None', 'Low', 'Medium', 'High']
    tier_go_counts = tier_go_counts.reindex(columns=expected_columns, fill_value=0)

    # 转化为百分比 (每一行加起来是 100%)
    tier_go_pct = tier_go_counts.div(tier_go_counts.sum(axis=1), axis=0) * 100

    # 3. 画堆叠柱状图
    fig, ax = plt.subplots(figsize=(9, 7))
    tier_go_pct.plot(kind='bar', stacked=True, ax=ax, colormap='YlGnBu', edgecolor='black')

    plt.title('GO Evidence Distribution across Spacedust Tiers', fontsize=14, pad=20)
    plt.xlabel('Spacedust Conservation Tier', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.xticks(rotation=0)

    # 添加每层样本量标注
    tier_totals = tier_go_counts.sum(axis=1)
    for i, (tier, total) in enumerate(tier_totals.items()):
        ax.text(i, 108, f'n={total}', ha='center', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 110)  # 扩展 Y 轴上限以容纳标注

    # 把图例移到外面
    plt.legend(title='GO Strength', bbox_to_anchor=(1.05, 1), loc='upper left')
    sns.despine()

    plt.savefig(os.path.join(species_config['output_figure_dir'], 'spacedust_tier_go_enrichment.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print("各梯队内 High 置信度的真实占比：")
    print(tier_go_pct['High'])







    # ==========================================
    # 新增：高保守性蛋白对 的 GO 和 STRING 置信度分布
    # ==========================================
    go_string_pairs = df_pairs
    go_string_pairs['pair_key'] = go_string_pairs.apply(
        lambda row: "_".join(sorted([str(row[STRING_P1_COL]), str(row[STRING_P2_COL])])), axis=1
    )
    STRING_SCORE_COL = "original_combined_score"
    go_string_pairs = go_string_pairs.drop_duplicates(subset=['pair_key'], keep='first')
    merged_go_string_df = pd.merge(
        spacedust_scores,
        go_string_pairs[['pair_key', STRING_SCORE_COL, GO_STRENGTH_COL, 'combined_score']],
        on='pair_key',
        how='outer'
    )
    merged_go_string_df[STRING_SCORE_COL] = (merged_go_string_df[STRING_SCORE_COL] / 1000.0).fillna(0)
    merged_go_string_df[GO_STRENGTH_COL] = merged_go_string_df[GO_STRENGTH_COL].fillna('None').astype(str)
    # combined_score 列也需要 fillna(0)，避免在 groupby 时出现问题
    merged_go_string_df['combined_score'] = (merged_go_string_df['combined_score'] / 1000.0).fillna(0)

    high_conservation_df = merged_go_string_df[merged_go_string_df['spacedust_score'] > 150].copy()
    low_conservation_df = merged_go_string_df[merged_go_string_df['spacedust_score'] <= 150].copy()
    
    print(f"\n高保守性蛋白对 (spacedust_score > 150): {len(high_conservation_df):,} 对")
    print(f"低保守性蛋白对 (spacedust_score ≤ 150): {len(low_conservation_df):,} 对")
    
    

    # STRING 分数分布对比（箱线图）
    plt.figure(figsize=(10, 6))
    
    high_string_scores = high_conservation_df[STRING_SCORE_COL]
    low_string_scores = low_conservation_df[STRING_SCORE_COL]
    
    bp = plt.boxplot([high_string_scores.dropna(), low_string_scores.dropna()],
                     labels=[f'High Conservation\n(n={len(high_string_scores):,})',
                             f'Low Conservation\n(n={len(low_string_scores):,})'],
                     patch_artist=True, widths=0.5)
    
    colors = ['#005A32', '#D9D9D9']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    plt.ylabel('STRING Combined Score (0-1)', fontsize=12)
    plt.title('STRING Functional Evidence: High vs. Low Conservation Pairs(Threshold: 150)', fontsize=14, pad=15)
    plt.grid(True, linestyle='--', alpha=0.5, axis='y')
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'high_vs_low_conservation_string_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # ==========================================
    # GO x STRING 置信度热力图
    # ==========================================
    print("\n正在生成 GO x STRING 置信度热力图...")
    
    # 定义 STRING 置信度区间
    def string_confidence_level(score):
        if score <= 0.5:
            return '0-0.5'
        elif score <= 0.6:
            return '0.5-0.6'
        elif score <= 0.7:
            return '0.6-0.7'
        elif score <= 0.8:
            return '0.7-0.8'
        elif score <= 0.9:
            return '0.8-0.9'
        else:
            return '0.9-1.0'
    
    # 创建 STRING 置信度列
    merged_go_string_df['string_confidence'] = merged_go_string_df[STRING_SCORE_COL].apply(string_confidence_level)
    
    # 定义顺序
    go_order = ['None', 'Low', 'Medium', 'High']
    string_order = ['0-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0']
    
    # 使用 merged_go_string_df 并过滤出在 intersection_pairs 中的 pair_key
    # 这样确保有 Spacedust、STRING 和 GO 三个数据源
    intersection_keys = set(intersection_pairs['pair_key'])
    heatmap_data = merged_go_string_df[
        (merged_go_string_df['pair_key'].isin(intersection_keys)) &
        (merged_go_string_df[GO_STRENGTH_COL].notna())
    ].copy()
    
    print(f"热力图使用数据：{len(heatmap_data):,} 对（同时有 Spacedust、STRING 和 GO 数据）")
    
    # 按 GO 和 STRING 置信度分组，计算样本量和 spacedust 平均分
    heatmap_summary = heatmap_data.groupby(
        [GO_STRENGTH_COL, 'string_confidence']
    ).agg(
        count=('pair_key', 'count'),
        mean_spacedust=('spacedust_score', 'mean')
    ).reset_index()
    
    # 创建透视表
    count_pivot = heatmap_summary.pivot(
        index=GO_STRENGTH_COL,
        columns='string_confidence',
        values='count'
    ).reindex(index=go_order, columns=string_order).fillna(0).astype(int)
    
    mean_pivot = heatmap_summary.pivot(
        index=GO_STRENGTH_COL,
        columns='string_confidence',
        values='mean_spacedust'
    ).reindex(index=go_order, columns=string_order).fillna(0)
    
    # ==========================================
    # 卡方检验：GO_Strength 与 STRING_Confidence 独立性检验
    # ==========================================
    from scipy.stats import chi2_contingency
    
    # 构建 4x5 列联表（合并 STRING 最后两个区间为 0.8-1.0）
    string_conf_merged = heatmap_data['string_confidence'].apply(
        lambda x: '0.8-1.0' if x in ['0.8-0.9', '0.9-1.0'] else x
    )
    contingency_table = pd.crosstab(
        heatmap_data[GO_STRENGTH_COL],
        string_conf_merged
    ).reindex(index=go_order, columns=['0-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-1.0']).fillna(0).astype(int)
    
    print("\n" + "="*60)
    print("卡方检验数据准备 (GO x STRING Independence Test)")
    print("="*60)
    print(f"列联表维度：4×5")
    print(contingency_table)
    
    # 检查并删除全零的行和列
    row_sums = contingency_table.sum(axis=1)
    col_sums = contingency_table.sum(axis=0)
    
    valid_rows = row_sums > 0
    valid_cols = col_sums > 0
    
    print(f"\n有效行数：{valid_rows.sum()}/{len(valid_rows)}")
    print(f"有效列数：{valid_cols.sum()}/{len(valid_cols)}")
    
    # 检查是否可以进行卡方检验
    chi2_result = None
    p_value = None
    dof = None
    cramers_v = None
    effect_interpretation = None
    chi2_failed = False
    failure_reason = None
    
    if valid_rows.sum() >= 2 and valid_cols.sum() >= 2:
        reduced_table = contingency_table.loc[valid_rows, valid_cols]
        print(f"\n简化后的列联表维度：{reduced_table.shape[0]}×{reduced_table.shape[1]}")
        
        try:
            chi2, p_value, dof, expected = chi2_contingency(reduced_table)
            
            n = reduced_table.sum().sum()
            if n > 0:
                phi2 = chi2 / n
                r, k = reduced_table.shape
                cramers_v = np.sqrt(phi2 / min(k-1, r-1))
                
                if cramers_v < 0.1:
                    effect_interpretation = "Negligible"
                elif cramers_v < 0.3:
                    effect_interpretation = "Small"
                elif cramers_v < 0.5:
                    effect_interpretation = "Medium"
                else:
                    effect_interpretation = "Large"
            
            chi2_result = chi2
            print(f"\n✓ 卡方检验成功")
            print(f"χ² = {chi2:.4f}, df = {dof}, p = {p_value:.4e}")
            print(f"Cramér's V = {cramers_v:.4f} ({effect_interpretation})")
            
        except ValueError as e:
            chi2_failed = True
            failure_reason = "Data sparsity (zero expected frequency)"
            print(f"\n✗ 卡方检验失败：期望频数为零")
    else:
        chi2_failed = True
        failure_reason = "Insufficient data (valid rows<2 or cols<2)"
        print(f"\n✗ 无法进行卡方检验：{failure_reason}")
    
    print("="*60)
    
    # 绘制热力图
    fig, ax = plt.subplots(figsize=(11, 9))
    
    sns.heatmap(mean_pivot, annot=True, fmt='.1f', cmap='YlOrRd',
                cbar_kws={'label': 'Mean Spacedust Score'}, ax=ax,
                annot_kws={'size': 14, 'weight': 'bold'})
    
    for i, go_level in enumerate(go_order):
        for j, string_level in enumerate(string_order):
            count_val = count_pivot.loc[go_level, string_level]
            if count_val > 0:
                ax.text(j + 0.5, i + 0.85, f'n={count_val}',
                       ha='center', va='center', fontsize=9, color='black', alpha=0.7)
    
    ax.set_xlabel('STRING Confidence Level', fontsize=12, fontweight='bold')
    ax.set_ylabel('GO Evidence Strength', fontsize=12, fontweight='bold')
    ax.set_title('GO x STRING Confidence Heatmap\n(Mean Spacedust Score & Sample Size)',
                fontsize=14, fontweight='bold', pad=15)
    ax.set_xticklabels(string_order, rotation=0, ha='center')
    ax.set_yticklabels(go_order, rotation=0, va='center')
    
    # 添加文本框：成功显示统计量，失败显示原因
    if not chi2_failed and chi2_result is not None:
        stats_text = (
            f"Chi-Square Test: χ² = {chi2_result:.2f}, df = {dof}, p = {p_value:.2e}  |  "
            f"Cramér's V = {cramers_v:.3f} ({effect_interpretation})"
        )
    else:
        stats_text = f"Chi-Square Test: Not applicable — {failure_reason}"
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(0.5, 0.02, stats_text, fontsize=11, ha='center', va='bottom',
             bbox=props, fontfamily='monospace')
    
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'go_string_confidence_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"热力图已保存：go_string_confidence_heatmap.png")
    print(f"总样本量（intersection_pairs）：{len(heatmap_data):,}")
    

 
    

    # ==========================================
    # Spacedust Score 区间 vs STRING Score 箱线图
    # ==========================================
    print("\n正在生成 Spacedust Score 区间 vs STRING Score 箱线图...")
    
    # 使用与上面相同的梯队划分方法
    merged_go_string_df['Spacedust_Tier'] = pd.cut(
        merged_go_string_df['spacedust_score'],
        bins=[0, 1e-5, 3, 10, 20, 40, 60, 80, 110, 140, 170, float('inf')],
        labels=['0', '1~2', '3~9', '10~19', '20~39', '40~59', '60~79', '80~109', '110~139', '140~169', '170+'],
        right=False
    )
    
    # 将 STRING 分数从 0-1 转换回 0-1000 以便阅读
    merged_go_string_df['string_score_scaled'] = merged_go_string_df[STRING_SCORE_COL] * 1000
    
    # 绘制箱线图
    plt.figure(figsize=(12, 7))
    
    # 使用箱线图展示 STRING 分数在每个 Spacedust 区间的分布
    order = ['0', '1~2', '3~9', '10~19', '20~39', '40~59', '60~79', '80~109', '110~139', '140~169', '170+']
    available_tiers = [t for t in order if t in merged_go_string_df['Spacedust_Tier'].values]
    
    sns.boxplot(x='Spacedust_Tier', y='original_combined_score', data=merged_go_string_df,
                order=available_tiers, palette='Blues', width=0.6, showfliers=False)
    
    # 添加样本量标注
    tier_counts = merged_go_string_df.groupby('Spacedust_Tier').size()
    for i, tier in enumerate(available_tiers):
        n = tier_counts.get(tier, 0)
        plt.text(i, plt.gca().get_ylim()[1] * 1.05, f'n={n:,}', ha='center', fontsize=10, fontweight='bold')
    
    plt.title('STRING Original Combined Score Distribution across Spacedust Conservation Tiers', fontsize=14, pad=20)
    plt.xlabel('Spacedust Conservation Score Tier', fontsize=12)
    plt.ylabel('STRING Combined Score (0-1)', fontsize=12)
    plt.xticks(rotation=45)
    plt.ylim(0, 1.1)  # STRING 分数范围是 0-1
    
    sns.despine()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'spacedust_tier_string_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()




    #再画一个剔除后的score箱线图
    # 绘制箱线图
    plt.figure(figsize=(12, 7))
    
    sns.boxplot(x='Spacedust_Tier', y='combined_score', data=merged_go_string_df,
                order=available_tiers, palette='Blues', width=0.6, showfliers=False)
    
    # 添加样本量标注
    tier_counts = merged_go_string_df.groupby('Spacedust_Tier').size()
    for i, tier in enumerate(available_tiers):
        n = tier_counts.get(tier, 0)
        plt.text(i, plt.gca().get_ylim()[1] * 1.05, f'n={n:,}', ha='center', fontsize=10, fontweight='bold')
    
    plt.title('STRING Revised Combined Score Distribution across Spacedust Conservation Tiers', fontsize=14, pad=20)
    plt.xlabel('Spacedust Conservation Score Tier', fontsize=12)
    plt.ylabel('STRING Combined Score (0-1)', fontsize=12)
    plt.xticks(rotation=45)
    plt.ylim(0, 1.1)  # STRING 分数范围是 0-1
    
    sns.despine()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'spacedust_tier_revised_string_boxplot.png'), dpi=300, bbox_inches='tight')
    plt.close()




       
    # ==========================================
    # 双 Y 轴图：STRING 均分 vs GO 注释比例
    # ==========================================
    print("\n正在生成双 Y 轴图：STRING 均分 vs GO 注释比例...")
    
    # 按 Spacedust 梯队分组统计
    
    dual_axis_summary = merged_go_string_df.groupby('Spacedust_Tier').agg(
        string_mean=(STRING_SCORE_COL, 'mean'),  # STRING 平均分 (0-1)
        go_annotated=('GO_Strength', lambda x: (x != 'None').sum() / len(x)),  # GO 有注释比例
        go_high=('GO_Strength', lambda x: (x == 'High').sum() / len(x)),  # GO 高置信度比例
        count=('pair_key', 'count')  # 样本量
    ).reset_index()
    
    # 创建图形和双 Y 轴
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    # 左 Y 轴 - STRING 平均分（柱状图）
    bar_width = 0.35
    x = range(len(dual_axis_summary))
    bars = ax1.bar(x, dual_axis_summary['string_mean'], bar_width,
                   color='#2c7bb6', alpha=0.8, label='STRING Mean Score')
    
    ax1.set_xlabel('Spacedust Conservation Tier', fontsize=12, fontweight='bold')
    ax1.set_ylabel('STRING Mean Score (0-1)', fontsize=12, fontweight='bold', color='#2c7bb6')
    ax1.tick_params(axis='y', labelcolor='#2c7bb6')
    ax1.set_ylim(0, 1.1)
    
    # 在柱子上方标注数值
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9, color='#2c7bb6', fontweight='bold')
    
    # 右 Y 轴 - GO 有注释比例（折线图）
    ax2 = ax1.twinx()
    
    # 第一条折线：GO 有注释比例
    line1 = ax2.plot(x, dual_axis_summary['go_annotated'], color='#d7191c',
                    marker='o', markersize=10, linewidth=3, label='GO Annotated %')
    
    ax2.set_ylabel('GO Proportion', fontsize=12, fontweight='bold', color='#d7191c')
    ax2.tick_params(axis='y', labelcolor='#d7191c')
    ax2.set_ylim(0, 1.1)
    
    # 在折线点上方标注数值
    for i, val in enumerate(dual_axis_summary['go_annotated']):
        ax2.text(i, val + 0.03, f'{val:.1%}', ha='center', va='bottom',
                fontsize=9, color='#d7191c', fontweight='bold')
    
    # 第二条折线：GO 高置信度比例
    line2 = ax2.plot(x, dual_axis_summary['go_high'], color="#d9f009",
                    marker='s', markersize=8, linewidth=2, linestyle='--', label='GO High Confidence %')
    
    # 在第二条折线点上方标注数值
    for i, val in enumerate(dual_axis_summary['go_high']):
        ax2.text(i, val + 0.03, f'{val:.1%}', ha='center', va='bottom',
                fontsize=8, color="#d9f009", fontweight='bold')
    
    # X 轴标签
    tier_labels = dual_axis_summary['Spacedust_Tier'].astype(str)
    ax1.set_xticks(x)
    ax1.set_xticklabels(tier_labels, rotation=45, ha='right')
    
    # 标题和图例
    plt.title('STRING Original Score & GO Annotation Rate across Spacedust Tiers',
              fontsize=14, fontweight='bold', pad=15)
    
    # 合并图例
    bars_legend = ax1.bar([], [], color='#2c7bb6', alpha=0.8, label='STRING Mean Score')
    line1_legend = ax2.plot([], [], color='#d7191c', marker='o', markersize=10, linewidth=3,
                          label='GO Annotated %')[0]
    line2_legend = ax2.plot([], [], color='#1f78b4', marker='s', markersize=8, linewidth=2, linestyle='--',
                          label='GO High Confidence %')[0]
    ax1.legend(handles=[bars_legend, line1_legend, line2_legend], loc='upper left')
    
    # 样本量标注
    for i, n in enumerate(dual_axis_summary['count']):
        ax1.text(i, -0.08, f'n={n:,}', ha='center', va='top', fontsize=8, color='gray')
    
    fig.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'string_go_dual_axis_trend.png'), dpi=300, bbox_inches='tight')
    plt.close()



    #再画一个剔除后的score双Y轴图
    dual_axis_summary = merged_go_string_df.groupby('Spacedust_Tier').agg(
        string_mean=("combined_score", 'mean'),  # STRING 平均分 (0-1)
        go_annotated=('GO_Strength', lambda x: (x != 'None').sum() / len(x)),  # GO 有注释比例
        go_high=('GO_Strength', lambda x: (x == 'High').sum() / len(x)),  # GO 高置信度比例
        count=('pair_key', 'count')  # 样本量
    ).reset_index()
    # 创建图形和双 Y 轴
    fig, ax1 = plt.subplots(figsize=(12, 7))
    
    # 左 Y 轴 - STRING 平均分（柱状图）
    bar_width = 0.35
    x = range(len(dual_axis_summary))
    bars = ax1.bar(x, dual_axis_summary['string_mean'], bar_width,
                   color='#2c7bb6', alpha=0.8, label='STRING Mean Score')
    
    ax1.set_xlabel('Spacedust Conservation Tier', fontsize=12, fontweight='bold')
    ax1.set_ylabel('STRING Mean Score (0-1)', fontsize=12, fontweight='bold', color='#2c7bb6')
    ax1.tick_params(axis='y', labelcolor='#2c7bb6')
    ax1.set_ylim(0, 1.1)
    
    # 在柱子上方标注数值
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9, color='#2c7bb6', fontweight='bold')
    
    # 右 Y 轴 - GO 有注释比例（折线图）
    ax2 = ax1.twinx()
    
    # 第一条折线：GO 有注释比例
    line1 = ax2.plot(x, dual_axis_summary['go_annotated'], color='#d7191c',
                    marker='o', markersize=10, linewidth=3, label='GO Annotated %')
    
    ax2.set_ylabel('GO Proportion', fontsize=12, fontweight='bold', color='#d7191c')
    ax2.tick_params(axis='y', labelcolor='#d7191c')
    ax2.set_ylim(0, 1.1)
    
    # 在折线点上方标注数值
    for i, val in enumerate(dual_axis_summary['go_annotated']):
        ax2.text(i, val + 0.03, f'{val:.1%}', ha='center', va='bottom',
                fontsize=9, color='#d7191c', fontweight='bold')
    
    # 第二条折线：GO 高置信度比例
    line2 = ax2.plot(x, dual_axis_summary['go_high'], color="#d9f009",
                    marker='s', markersize=8, linewidth=2, linestyle='--', label='GO High Confidence %')
    
    # 在第二条折线点上方标注数值
    for i, val in enumerate(dual_axis_summary['go_high']):
        ax2.text(i, val + 0.03, f'{val:.1%}', ha='center', va='bottom',
                fontsize=8, color="#d9f009", fontweight='bold')
    
    # X 轴标签
    tier_labels = dual_axis_summary['Spacedust_Tier'].astype(str)
    ax1.set_xticks(x)
    ax1.set_xticklabels(tier_labels, rotation=45, ha='right')
    
    # 标题和图例
    plt.title('STRING Revised Score & GO Annotation Rate across Spacedust Tiers',
              fontsize=14, fontweight='bold', pad=15)
    
    # 合并图例
    bars_legend = ax1.bar([], [], color='#2c7bb6', alpha=0.8, label='STRING Mean Score')
    line1_legend = ax2.plot([], [], color='#d7191c', marker='o', markersize=10, linewidth=3,
                          label='GO Annotated %')[0]
    line2_legend = ax2.plot([], [], color='#1f78b4', marker='s', markersize=8, linewidth=2, linestyle='--',
                          label='GO High Confidence %')[0]
    ax1.legend(handles=[bars_legend, line1_legend, line2_legend], loc='upper left')
    
    # 样本量标注
    for i, n in enumerate(dual_axis_summary['count']):
        ax1.text(i, -0.08, f'n={n:,}', ha='center', va='top', fontsize=8, color='gray')
    
    fig.tight_layout()
    plt.savefig(os.path.join(species_config['output_figure_dir'], 'revised_string_go_dual_axis_trend.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    
    print(f"双 Y 轴图已保存：revised_string_go_dual_axis_trend.png")
    merged_go_string_df.to_csv(os.path.join(species_config['output_figure_dir'], 'merged_go_string_data.csv'), index=False)
    






    # 计算并打印每个区间的 STRING 分数统计信息
    print("\n各 Spacedust 梯队内 STRING 分数的统计信息：")
    tier_string_stats = merged_go_string_df.groupby('Spacedust_Tier')['string_score_scaled'].agg(['count', 'mean', 'median', 'std'])
    print(tier_string_stats)




    # 从 intersection_pairs 中减去 STRING 无数据的样本，得到真正有 STRING 数据的蛋白对
    intersection_with_string_data = intersection_pairs[intersection_pairs['combined_score'] > 0]
    plot_spacedust_pairs_venn(len(spacedust_scores), intersection_go_data, intersection_with_string_data, species_config)
    
    





    # 输出统计摘要
    print("\n" + "="*60)
    print("高保守性蛋白对 (spacedust_score > 150) 详细统计")
    print("="*60)
    print(f"总数：{len(high_conservation_df):,} 对")
    print(f"\nGO 证据强度分布:")
    for level in ['High', 'Medium', 'Low', 'None']:
        count = len(high_conservation_df[high_conservation_df[GO_STRENGTH_COL] == level])
        pct = count / len(high_conservation_df) * 100
        print(f"  {level}: {count:,} ({pct:.1f}%)")
    
    print(f"\nSTRING 功能证据分布:")
    string_annotated = high_conservation_df[high_conservation_df['combined_score'] > 0]
    print(f"  有 STRING 数据：{len(string_annotated):,} ({len(string_annotated)/len(high_conservation_df)*100:.1f}%)")
    if len(string_annotated) > 0:
        print(f"  STRING 平均分：{string_annotated[STRING_SCORE_COL].mean():.3f}")
        print(f"  STRING 中位数：{string_annotated[STRING_SCORE_COL].median():.3f}")
        high_string = len(string_annotated[string_annotated[STRING_SCORE_COL] > 0.7])
        print(f"  高置信度 STRING(>0.7): {high_string:,} ({high_string/len(string_annotated)*100:.1f}%)")
    
    # 计算高保守性组中 GO 和 STRING 同时高置信度的比例
    go_high = high_conservation_df[high_conservation_df[GO_STRENGTH_COL] == 'High']
    go_high_pct = len(go_high) / len(high_conservation_df) * 100
    string_high = string_annotated[string_annotated[STRING_SCORE_COL] > 0.7] if len(string_annotated) > 0 else pd.DataFrame()
    string_high_pct = len(string_high) / len(string_annotated) * 100 if len(string_annotated) > 0 else 0
    
    print(f"\n高置信度对比:")
    print(f"  GO High 置信度：{go_high_pct:.1f}%")
    print(f"  STRING 高置信度 (>0.7): {string_high_pct:.1f}%")
    print("="*60 + "\n")



    print("\n" + "="*60)
    print("数据统计报告")
    print("="*60)
    print(f"Spacedust 预测的蛋白对总数：{total_spacedust_pairs:,}")
    print(f"pairs数据库的蛋白对总数：{len(string_pairs):,}")
    print(f"两者交集（同时存在）：{len(intersection_pairs):,} ({len(intersection_pairs)/total_spacedust_pairs*100:.1f}% of Spacedust)")
    print(f"仅 Spacedust 有（Pairs 无）：{len(only_spacedust):,} ({len(only_spacedust)/total_spacedust_pairs*100:.1f}% of Spacedust)")
    print(f"仅 Pairs 有（Spacedust 无）：{len(only_string):,}")
    print(f"交集中STRING无数据的样本量：{len(intersection_no_string_data):,}")
    print(f"交集中GO无数据的样本量：{len(intersection_no_go_data):,}")
    print("="*60 + "\n")
    


