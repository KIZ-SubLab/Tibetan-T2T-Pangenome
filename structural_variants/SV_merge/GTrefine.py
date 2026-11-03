import os
import sys
import pysam
from intervaltree import IntervalTree

def load_bed_to_tree(bed_path):
    """将单个 BED 文件加载为 Chromosome -> IntervalTree 的字典"""
    tree_dict = {}
    if not os.path.exists(bed_path):
        return None
        
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            lines = line.strip().split('\t')
            chrom = lines[0]
            start = int(lines[1])
            end = int(lines[2])
            
            if chrom not in tree_dict:
                tree_dict[chrom] = IntervalTree()
            tree_dict[chrom].addi(start, end)
    return tree_dict

def check_overlap(chrom, start, end, tree_dict, min_overlap_ratio=0.8):
    """检查变异区间是否与 BED 区间满足设定的最小重叠度比例"""
    if not tree_dict or chrom not in tree_dict:
        return False
        
    sv_len = end - start
    if sv_len <= 0:
        return True
        
    overlaps = tree_dict[chrom].overlap(start, end)
    if not overlaps:
        return False
        
    total_intersection = 0
    for interval in overlaps:
        o_start = max(start, interval.begin)
        o_end = min(end, interval.end)
        if o_end > o_start:
            total_intersection += (o_end - o_start)
            
    return (total_intersection / sv_len) >= min_overlap_ratio

def batch_refine_vcf(in_vcf_path, out_vcf_path, bed_dir, min_overlap=0.8):
    """
    全自动批量修正 VCF 中所有样本的基因型，并规范化分隔符
    """
    # 1. 打开 VCF，获取所有样本名
    in_vcf = pysam.VariantFile(in_vcf_path)
    all_samples = list(in_vcf.header.samples)
    print(f"Found {len(all_samples)} samples in VCF header.")
    
    # 2. 批量预加载所有样本的 IntervalTree 到内存
    print("Pre-loading all sample BED files into memory...")
    sample_trees = {}
    active_samples = []
    
    for sample in all_samples:
        hap1_path = os.path.join(bed_dir, f"{sample}_hap1.bed")
        hap2_path = os.path.join(bed_dir, f"{sample}_hap2.bed")
        
        t1 = load_bed_to_tree(hap1_path)
        t2 = load_bed_to_tree(hap2_path)
        
        if t1 and t2:
            sample_trees[sample] = {'hap1': t1, 'hap2': t2}
            active_samples.append(sample)
        else:
            print(f" -> Missing BED for {sample} ({sample}_hap1.bed or _hap2.bed not found). Skipping refinement.")
            
    print(f"Initialization done. Total {len(active_samples)} / {len(all_samples)} samples will be processed.")
    
    # 3. 设置输出 VCF 头文件
    header = in_vcf.header
    out_vcf = pysam.VariantFile(out_vcf_path, 'w', header=header)
    
    # 统计每个样本纠正了多少个变异点
    refine_counts = {sample: 0 for sample in active_samples}
    
    print("Processing VCF lines...")
    line_count = 0
    
    for record in in_vcf:
        line_count += 1
        if line_count % 10000 == 0:
            print(f"Processed {line_count} variants...")
            
        chrom = record.chrom
        pos = record.pos - 1  # pysam 读入自动转为 0-based
        
        # 提取变异的结束位置 (支持精确解析 SVLEN)
        sv_len = 1
        if 'SVLEN' in record.info:
            val = record.info['SVLEN']
            sv_len = abs(val[0]) if isinstance(val, (tuple, list)) else abs(val)
        elif record.stop:
            sv_len = record.stop - record.pos
        end = pos + sv_len
        
        # 遍历有效样本列
        for sample in active_samples:
            sample_data = record.samples[sample]
            gt = sample_data['GT']
            
            # 如果包含缺失值 '.'
            if gt and (None in gt):
                allele_1, allele_2 = gt
                is_changed = False
                
                # 记录原本是不是全缺失的基因型 (./. 或 .|.)
                was_double_missing = (allele_1 is None and allele_2 is None)
                
                # 检查第一单倍型
                if allele_1 is None:
                    if check_overlap(chrom, pos, end, sample_trees[sample]['hap1'], min_overlap):
                        allele_1 = 0
                        is_changed = True
                        
                # 检查第二单倍型
                if allele_2 is None:
                    if check_overlap(chrom, pos, end, sample_trees[sample]['hap2'], min_overlap):
                        allele_2 = 0
                        is_changed = True
                
                # 根据修正状态更新 GT 和分隔符
                if is_changed:
                    sample_data['GT'] = (allele_1, allele_2)
                    
                    # 如果原本是全缺失，且两边都没能成功纠正（虽然理论上is_changed为True时至少有一边变了，这里做兜底保护）
                    if was_double_missing and (allele_1 is None and allele_2 is None):
                        sample_data.phased = False
                    else:
                        # 只要有任意一边成功被 BED 覆盖并修正为 0，就应当转化为物理相位 |
                        sample_data.phased = True
                        
                    refine_counts[sample] += 1
                else:
                    # 如果没有触发任何修正，且原先是全缺失状态（如 .|.），强行统一降级为未相位的 ./.
                    if allele_1 is None and allele_2 is None:
                        sample_data.phased = False
                        
        out_vcf.write(record)
        
    in_vcf.close()
    out_vcf.close()
    
    print("\n=== Refinement Report ===")
    for sample, count in refine_counts.items():
        print(f"Sample {sample}: Refined {count} alleles to 0.")
    print("Everything done!")

# --- 运行入口 ---
if __name__ == "__main__":
    # 在这里配置你的实际路径
    INPUT_VCF = "../02-sampel_merge/mix_merge.sv.survivor_1oAll.sort.ctRM.truvari_collapse.invRM.uniform.vcf"
    OUTPUT_VCF = "mix_merge.sortsv.ctrm.gtrefine.vcf"
    BED_DIRECTORY = "./bed"  # 存放所有样本两个单倍型 bed 的目录
    
    batch_refine_vcf(
        in_vcf_path = INPUT_VCF,
        out_vcf_path = OUTPUT_VCF,
        bed_dir = BED_DIRECTORY,
        min_overlap = 1.0  # 重叠度阈值
    )
