# Tibetan T2T Pangenome

This repository accompanies the Tibetan T2T pangenome study by Yaoxi He and colleagues. The project generated 70 phased haplotype assemblies from 35 Tibetan trios, including 16 near telomere-to-telomere haplotypes. The Tibetan pangenome graph represents 122.05 Mb of sequence beyond the T2T-CHM13 linear reference. A joint graph integrates the Tibetan assemblies with 870 haplotypes from APGp1, HGSVC3, and HPRC2 for population comparisons.

The study examines genomic diversity in Tibetans and investigates how structural variation and archaic ancestry contribute to candidate signals of high-altitude adaptation. The analyses span population differentiation, trait associations, and regulatory annotation.

## Contents

- [Repository structure](#repository-structure)
- [Data resources](#data-resources)
- [Citation](#citation)

## Repository structure

| Module | Analysis |
| --- | --- |
| [Genome assembly](assembly/) | Read preparation, parental k-mer profiling, trio phasing with hifiasm, and nextPolish2 polishing. |
| [Assembly quality control](qc/) | Yak quality values and trio evaluation, Inspector, and QUAST. |
| [Genome annotation](genome_annotation/) | Telomere architecture, telomeric variant repeats, centromeric satellites, and higher-order repeats. |
| [Pangenome construction and genotyping](pangenome/) | Minigraph-Cactus graphs, graph VCF processing, and PanGenie genotyping. |
| [Structural variant discovery and merging](structural_variants/) | SV calling from HiFi reads and phased assemblies, callset integration, and genotype refinement. |
| [Tandem repeat analysis](tandem_repeats/) | Tandem repeat genotyping, population differentiation, and repeat motif annotation for Figure 5. |
| [Selection simulations](selection/) | Coalescent simulations for empirical population branch statistic (PBS) significance. |
| [Archaic introgression](introgression/) | Archaic segment detection, SV genotyping, haplotype analyses, and local genealogies. |
| [Functional annotation and enrichment](functional_annotation/) | TESV regulatory annotation, hypoxia response models, and enrichment of candidate genes for GWAS traits and hypoxia pathways. |

## Data resources

| Resource | Archive or project record |
| --- | --- |
| Study BioProject | [NGDC BioProject PRJCA043737](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA043737) |
| Short-read whole-genome sequencing, PacBio HiFi reads, and ultra-long Oxford Nanopore reads | [Genome Sequence Archive for Human (GSA-Human)](https://ngdc.cncb.ac.cn/gsa-human/) |
| Phased haplotype assemblies | [Genome Warehouse (GWH)](https://ngdc.cncb.ac.cn/gwh/) |
| Variant calls | [Genome Variation Map (GVM)](https://ngdc.cncb.ac.cn/gvm/) |
| Tibetan pangenome graph | [NGDC BioProject PRJCA043737](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA043737) |

Researchers can request controlled access to raw reads and individual-level genome data through the Data Access Committee of the Kunming Institute of Zoology, Chinese Academy of Sciences. Requests include a research proposal and institutional ethics approval.

### External datasets

- **Reference assemblies and comparison panels.** T2T-CHM13 v2.0 provides the primary graph backbone. Comparative haplotypes come from APGp1, HGSVC3, and HPRC2. [T2T-CHM13](https://github.com/marbl/CHM13) · [APGp1](https://github.com/Asian-Pan-Genome/APGp1)
- **Archaic genomes and modern outgroup.** Neanderthal and Denisovan genomes support archaic ancestry comparisons. Yoruba (YRI) samples from the 1000 Genomes Project serve as a modern outgroup. [Archaic genome projects](https://www.eva.mpg.de/genetics/genome-projects/) · [International Genome Sample Resource](https://www.internationalgenome.org/)
- **Regulatory annotation.** ENCODE4 candidate cis-regulatory elements and their target-gene links complement the study's hypoxia time-course data. [ENCODE](https://www.encodeproject.org/)

## Citation

He, Y. et al. *Tibetan near-complete pangenome reveals complex variants underlying high-altitude adaptation*. bioRxiv (2025). [doi:10.64898/2025.12.16.694547](https://doi.org/10.64898/2025.12.16.694547)
