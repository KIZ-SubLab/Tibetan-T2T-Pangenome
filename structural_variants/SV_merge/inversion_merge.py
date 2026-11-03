import os
import glob
import pandas as pd
import numpy as np

def calculate_reciprocal_overlap(s1, e1, s2, e2):
    """计算两个区间的相互重叠比例 (Reciprocal Overlap >= 50%)"""
    overlap = min(e1, e2) - max(s1, s2)
    if overlap <= 0:
        return 0.0
    len1 = e1 - s1
    len2 = e2 - s2
    return min(overlap / len1, overlap / len2)

def find_connected_components(nodes, edges):
    """基于广度优先搜索(BFS)寻找图中的所有连通组件"""
    visited = set()
    components = []
    
    # 建立邻接表
    adj = {i: [] for i in range(len(nodes))}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        
    for node_idx in range(len(nodes)):
        if node_idx not in visited:
            component = []
            queue = [node_idx]
            visited.add(node_idx)
            
            while queue:
                curr = queue.pop(0)
                component.append(nodes[curr])
                for neighbor in adj[curr]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            components.append(component)
    return components

def merge_to_non_redundant_sites(all_sample_variants, sample_order):
    """
    通过图论连通组件（Connected Components）进行群体变异彻底去重。
    确保任何相互重叠 >= 50% 的区间都能被安全聚类到一起。
    """
    if not all_sample_variants:
        return []
    
    # 汇总所有记录
    df_all = pd.concat(all_sample_variants).reset_index(drop=True)
    
    # 建立样本优先级字典 (用于后续选择基准坐标)
    sample_priority = {s: i for i, s in enumerate(sample_order)}
    df_all['priority'] = df_all['sample'].map(sample_priority)
    
    unique_sites = []
    
    # 按染色体分开处理，防止跨染色体比对
    for chrom, df_chrom in df_all.groupby('chr'):
        nodes = df_chrom.to_dict('records')
        edges = []
        n_nodes = len(nodes)
        
        # 1. 两两比对，如果相互重叠 >= 50%，则建立一条边
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                s1, e1 = nodes[i]['start'], nodes[i]['start'] + nodes[i]['length']
                s2, e2 = nodes[j]['start'], nodes[j]['start'] + nodes[j]['length']
                
                if calculate_reciprocal_overlap(s1, e1, s2, e2) >= 0.5:
                    edges.append((i, j))
                    
        # 2. 寻找所有的连通组件
        components = find_connected_components(nodes, edges)
        
        # 3. 对每个组件进行塌陷，提取主坐标与属性
        for component in components:
            unique_sites.append(collapse_component_to_site(component, sample_order))
            
    return unique_sites

def collapse_component_to_site(component, sample_order):
    """
    在一个连通组件内：
    1. 严格筛选出符合 sample_order 优先级最高（最靠前）的样本坐标作为基准坐标。
    2. 如果最高优先级的样本有多个记录/Caller，取其中 start 最小的第一个。
    """
    sorted_component = sorted(component, key=lambda x: (x['priority'], x['start']))
    base = sorted_component[0]
    
    # 收集当前组件包含的所有非冗余 TOOL/Caller 名
    all_callers = set()
    for r in component:
        for c in str(r['SVcaller']).split(','):
            if c.strip():
                all_callers.add(c.strip())
                
    return {
        'chr': base['chr'],
        'start': base['start'],
        'end': base['start'] + base['length'],
        'length': base['length'],
        'TOOL': ",".join(sorted(list(all_callers))),
        'cluster_members': component
    }

def load_and_split_sample(sample_name, caller_dict):
    """从5个文件夹中读取单个样本并拆分 Hap，兼容多列输入"""
    all_caller_rows = []
    for caller_name, folder_path in caller_dict.items():
        bed_matches = glob.glob(os.path.join(folder_path, f"{sample_name}.bed*"))
        if not bed_matches: continue
        if caller_name.lower() in ['lgvar', 'pav']:
            target_cols = [0, 1, 2, 4]  # 读取 0, 1, 2, 4 列
        else:
            target_cols = [0, 1, 2, 3]
        df_c = pd.read_csv(bed_matches[0], sep=r'\s+', header=None, usecols=target_cols,
                           names=['chr', 'start', 'length', 'genotype'])
        df_c['SVcaller'] = caller_name
        all_caller_rows.append(df_c)
        
    if not all_caller_rows: return pd.DataFrame()
    
    df_all_callers = pd.concat(all_caller_rows, ignore_index=True)
    
    hap1_rows, hap2_rows = [], []
    for _, row in df_all_callers.iterrows():
        gt = str(row['genotype'])
        if gt in ['1/1', '1|1']:
            hap1_rows.append(row.copy())
            hap2_rows.append(row.copy())
        elif gt in ['1/0', '1|0', './1', '.|1']:
            hap1_rows.append(row.copy())
        elif gt in ['0/1', '0|1', '1/.', '1|.', '0/2', '1/2']: 
            hap2_rows.append(row.copy())
            
    df_h1 = pd.DataFrame(hap1_rows)
    df_h2 = pd.DataFrame(hap2_rows)
    if not df_h1.empty: df_h1['hap'] = 'Hap1'
    if not df_h2.empty: df_h2['hap'] = 'Hap2'
    
    df_sample = pd.concat([df_h1, df_h2], ignore_index=True)
    if not df_sample.empty:
        df_sample['sample'] = sample_name
    return df_sample

# ==================== 主程序流程 ====================
if __name__ == "__main__":
    
    # 1. 从外部文件读取样本名 List
    with open("sample_list.txt", "r") as f:
        sample_order = [line.strip() for line in f if line.strip()]
    
    # 需求 2: 在这里提供 5 个 Caller 的目录路径
    caller_folders = {
        "lgvar": "/public/home/zhoufeifei/project/APG_Xizang/SV_calling/inversion/01-lgvar",
        "pav": "/public/home/zhoufeifei/project/APG_Xizang/SV_calling/inversion/02-pav",
        "svimasm": "/public/home/zhoufeifei/project/APG_Xizang/SV_calling/inversion/03-svimasm",
        "pbsv": "/public/home/zhoufeifei/project/APG_Xizang/SV_calling/inversion/04-pbsv",
        "sniffles2": "/public/home/zhoufeifei/project/APG_Xizang/SV_calling/inversion/05-sniffles"  # 请替换为第5个caller的实际名称和路径
    }
    
    all_sample_variants = []
    print(">>> 步骤 1: 正在从 5 个文件夹中读取样本并切分 Haplotype (兼容多列数据)...")
    for s_name in sample_order:
        df_sample = load_and_split_sample(s_name, caller_folders)
        if not df_sample.empty:
            all_sample_variants.append(df_sample)
            
    print(">>> 步骤 2: 正在利用图论连通组件算法进行全群体合并去重...")
    print(sample_order)
    unique_sites = merge_to_non_redundant_sites(all_sample_variants, sample_order)
    
    print(">>> 步骤 3: 正在构建群体基因型矩阵...")
    matrix_rows = []
    hap_columns = []
    for s_name in sample_order:
        hap_columns.extend([f"{s_name}.1", f"{s_name}.2"])
        
    for site in unique_sites:
        row_data = {
            'TOOL': site['TOOL'],
            'CHR': site['chr'],
            'START': site['start'],
            'END': site['end'],
            'LENGTH': site['length']
        }
        
        for col in hap_columns:
            row_data[col] = 0
            
        for member in site['cluster_members']:
            s = member['sample']
            h = "1" if member['hap'] == "Hap1" else "2"
            col_name = f"{s}.{h}"
            row_data[col_name] = 1
            
        total_haps = len(hap_columns)
        count_1 = sum([row_data[col] for col in hap_columns])
        row_data['AF'] = count_1 / total_haps if total_haps > 0 else 0.0
        
        matrix_rows.append(row_data)
        
    df_matrix = pd.DataFrame(matrix_rows)
    
    if not df_matrix.empty:
        final_cols = ['TOOL', 'CHR', 'START', 'END', 'LENGTH', 'AF'] + hap_columns
        df_matrix = df_matrix[final_cols]
        df_matrix = df_matrix.sort_values(['CHR', 'START']).reset_index(drop=True)
        
        # 需求 3: 对最终的矩阵数据进行分流过滤
        
        # 3.1 过滤保留至少两个 Caller 支持的记录 (TOOL 列包含逗号说明有多个工具)
        df_multi_support = df_matrix[df_matrix['TOOL'].str.contains(',', na=False)].copy()
        
        # 3.2 过滤出仅由 lgvar 独有支持的记录 (TOOL 严格等于 'lgvar')
        df_lgvar_unique = df_matrix[df_matrix['TOOL'] == 'lgvar'].copy()
        
        # 保存主文件
        output_main = "population_inversion_multi_support_matrix.tsv"
        df_multi_support.to_csv(output_main, sep="\t", index=False)
        print(f"\n>>> [主文件] 至少2个Caller支持的矩阵已写入：{output_main} (共 {len(df_multi_support)} 条记录)")
        
        # 保存 lgvar 独有文件
        output_lgvar = "population_inversion_lgvar_unique_matrix.tsv"
        df_lgvar_unique.to_csv(output_lgvar, sep="\t", index=False)
        print(f">>> [独有文件] lgvar独有支持的矩阵已写入：{output_lgvar} (共 {len(df_lgvar_unique)} 条记录)")
        
    else:
        print("\n 错误: 未能在输入的 BED 文件中聚类出任何有效变异位点。")
