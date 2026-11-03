# Structural variant callset integration

The workflow integrates SV calls from PAV, SVIM-asm, LGvar, pbsv, and Sniffles2. It combines calls within samples, merges callsets across the cohort, and prepares simple SV and inversion datasets for population analyses.

| Script | Analysis |
| --- | --- |
| [01_caller_merge.sh](01_caller_merge.sh) | Combine per-sample caller VCFs with bcftools and SURVIVOR, harmonize genotypes, and take the union of both merged callsets. |
| [02_sample_merge.sh](02_sample_merge.sh) | Merge per-sample callsets, apply the target-region filter, and collapse redundant calls with Truvari. |
| [consensusSV_from_bcftoolsMerge.py](consensusSV_from_bcftoolsMerge.py) | Derive consensus calls and genotypes from the bcftools merge. |
| [correctGT_for_survivorSV.py](correctGT_for_survivorSV.py) | Harmonize genotypes in the SURVIVOR callset. |
| [unionSVsets.py](unionSVsets.py) | Combine the normalized bcftools and SURVIVOR callsets. |
| [GTrefine.py](GTrefine.py) | Refine missing genotypes using haplotype coverage intervals. |
| [inversion_merge.py](inversion_merge.py) | Cluster inversion calls by overlap and generate population haplotype matrices. |

Set the sample list, caller input paths, and helper directory in the caller merge script. Set the cohort input paths and target-region BED in the sample merge script. The workflow uses bcftools, SURVIVOR, Truvari, bgzip, tabix, and Python 3. Python dependencies include icecream, pysam, intervaltree, pandas, and NumPy.
