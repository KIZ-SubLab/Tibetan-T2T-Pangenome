# Pangenome construction and genotyping

The graph workflow constructs chromosome-level pangenomes from phased haplotype assemblies. The genotyping workflow uses a phased variant panel to genotype short-read samples.

| Workflow | Analysis |
| --- | --- |
| [Minigraph-Cactus](minigraph_cactus/) | Graph construction, followed by filtering with vcfbub and normalization with vcfwave. |
| [PanGenie](pangenie/) | Panel preparation, index construction, per-sample genotyping, and genotype VCF merging. |
