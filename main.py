"""
多物种批量分析主程序
整合 main-1.py 和 main-2.py 的完整分析流程
遍历 8 个物种，自动加载数据、执行分析、保存结果
"""

import os
import pandas as pd
import data_loader
import data_process
import querymapper
import comparison
import analyzer
import visualizer_batch
import species_config
from datetime import datetime


def analysis(species_key, config):
    print(f"\n{'='*60}")
    print(f"开始处理：{config['name']} ({species_key})")
    print(f"{'='*60}")

    df_goa = data_loader.load_goa(config['goa_path'])
    df_alias = data_loader.load_alias(config['aliases_path'])
    feature_table = pd.read_csv(config['feature_table_path'], sep='\t')
    df_links = data_loader.load_string_links(config['links_path'])

    match_table = data_process.match_bp_table(df_goa, df_alias, feature_table, config['taxon_id'], df_links)
    print(f"匹配后的 GOA 数据：{len(match_table)} 条记录")
    match_table.to_csv(config['output_data_dir'] + '/match_table.csv', index=False)

    #main2 做核心分析,生成GO-STRING蛋白对数据表 
    keep_channels = ['fusion', 'coexpression', 'experimental', 'database', 'textmining']

    df_links['original_combined_score'] = df_links['combined_score'] #保留原始综合分数，方便后续比较
    df_links['combined_score'] = df_links.apply(lambda row: data_process.recalculate_string_score(row, keep_channels), axis=1) #剔除部分分数后重新计算综合分数
    df_links['fusion_score'] = df_links['fusion']
    df_links['experimental_score'] = df_links['experimental']
    df_links['database_score'] = df_links['database']
    df_links['textmining_score'] = df_links['textmining']
    df_links['coexpression_score'] = df_links['coexpression']
    df_links['neighborhood_score'] = df_links['neighborhood']
    df_links['cooccurence_score'] = df_links['cooccurence']
    df_links = df_links[['protein1', 'protein2', 'original_combined_score', 'combined_score', 'fusion_score', 'experimental_score', 'database_score', 'textmining_score', 'neighborhood_score', 'cooccurence_score', 'coexpression_score']]
    df_links.to_csv(os.path.join(config['output_data_dir'], 'string_links_recalculated.csv'), index=False) #保存重新计算分数后的链接数据，方便后续分析和检查



    df_merged = analyzer.analyze_links(match_table, df_links)
    df_pairs = df_merged
    df_pairs = analyzer.add_gene_names(df_merged, feature_table, config['taxon_id'])
    df_pairs.to_csv(os.path.join(config['output_data_dir'], 'go-string_match_pairs.csv'), index=False)


    #main3 画OVERLAP评估基因覆盖情况
    df_aliases = data_loader.load_alias(config['aliases_path'])

    total, go, string = visualizer_batch.calculate_coverage_metrics(feature_table, match_table, df_aliases)
    visualizer_batch.plot_venn(total, go, string, config)
    total_pairs, go_pairs, string_pairs = visualizer_batch.calculate_pair_metrics(feature_table, df_pairs)
    visualizer_batch.plot_pairs_venn(total_pairs, go_pairs, string_pairs, config)
    total_pairs, high_go_pairs, high_string_pairs = visualizer_batch.calculate_high_confidence_pairs(feature_table, df_pairs)
    visualizer_batch.plot_high_confidence_pairs_venn(total_pairs, go_pairs, string_pairs, high_go_pairs, high_string_pairs, config)


    #main4和spacedust比对
    df_query = data_loader.load_query_lookup(config['querydb_lookup_path'])
    df_matrix = pd.read_csv(os.path.join(config['matrix_path']), index_col=0)

    df_lookup = querymapper.query_mapping(df_query, feature_table, config['taxon_id'], df_links)
    df_lookup.to_csv(os.path.join(config['output_data_dir'], 'query_mapping.csv'), index=False) #保存查询映射结果，方便后续分析和检查
    comparison.compare_pairs(df_matrix, df_lookup, df_pairs, config)



def main():
    """主函数：批量处理所有物种"""
    print("="*60)
    print("多物种批量分析程序")
    print(f"开始时间：{datetime.now().isoformat()}")
    print("="*60)
    
    # 获取所有物种配置
    all_species = list(species_config.SPECIES_CONFIG.keys())
    print(f"\n共发现 {len(all_species)} 个物种：{', '.join(all_species)}")
    
    # 遍历处理每个物种
    for i, species_key in enumerate(all_species, 1):
        config = species_config.SPECIES_CONFIG[species_key]
        print(f"\n>>> 进度：{i}/{len(all_species)}")
        
        try:
            analysis(species_key, config)
        except Exception as e:
            print(f"\n!!! 错误：处理 {species_key} 时发生异常：{str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    print(f"\n{'='*60}")
    print(f"批量处理完成！")
    print(f"结束时间：{datetime.now().isoformat()}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
