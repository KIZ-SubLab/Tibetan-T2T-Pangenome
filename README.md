# Tibetan-T2T-Pangenome

Analysis code for *Tibetan near-complete pangenome reveals complex variants underlying high-altitude adaptation* by Yaoxi He and colleagues.

## Study overview

The study used 35 Tibetan parent–offspring trios (105 participants) to generate 70 phased haplotype assemblies. These include 16 near-telomere-to-telomere (near-T2T) haplotypes from eight diploid genomes and 54 high-quality haplotypes from 27 diploid genomes. The Tibetan pangenome graph integrates the 70 haplotypes with T2T-CHM13 and contains 122.05 Mb of sequence absent from the linear reference. A separate joint pangenome combines the Tibetan assemblies with 870 haplotypes from APGp1, HGSVC3, and HPRC2 for population comparisons.

The manuscript characterizes small and structural variants, including variation in telomeres, centromeres, and the Y chromosome. It identifies population-differentiated insertions, deletions, and tandem repeats, and uses trait associations and regulatory annotation to assess candidates relevant to high-altitude adaptation. The archaic ancestry analysis recovers the known Denisovan-derived haplotype at `EPAS1` and identifies additional candidate introgressed regions, some of which carry Tibetan-enriched structural variants (TESVs).

## Data sources and access

The manuscript identifies **PRJCA043737** as the study BioProject at the National Genomics Data Center (NGDC). On 3 November 2026, NGDC listed the BioProject as unreleased. [NGDC BioProject PRJCA043737](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA043737)

| Study output | Archive named in the manuscript |
| --- | --- |
| Whole-genome sequencing reads, including short reads, PacBio HiFi, and ultra-long Oxford Nanopore reads | [Genome Sequence Archive for Human (GSA-Human)](https://ngdc.cncb.ac.cn/gsa-human/) |
| Phased haplotype assemblies | [Genome Warehouse (GWH)](https://ngdc.cncb.ac.cn/gwh/) |
| Variant calls | [Genome Variation Map (GVM)](https://ngdc.cncb.ac.cn/gvm/) |

The manuscript also lists the pangenome graph under this BioProject. It describes controlled access to raw sequencing and individual-level genome data for nonprofit scientific research. Requests go to the Data Access Committee of the Kunming Institute of Zoology, Chinese Academy of Sciences, with a research proposal and institutional ethics approval.

### Additional inputs

- **Reference and comparison genomes.** T2T-CHM13 v2.0 is the primary graph backbone. Some analyses use GRCh38 or GRCh37 coordinates. APGp1, HGSVC3, and HPRC2 assemblies provide comparative haplotypes. [T2T-CHM13 reference files](https://github.com/marbl/CHM13)
- **Archaic ancestry inputs.** High-coverage Neanderthal and Denisovan genomes are used for archaic comparisons, with Yoruba (YRI) samples from the 1000 Genomes Project as a modern outgroup. [Archaic genome projects](https://www.eva.mpg.de/genetics/genome-projects/), [International Genome Sample Resource](https://www.internationalgenome.org/)
- **Regulatory and phenotypic inputs.** TESV annotation uses ENCODE4 candidate cis-regulatory elements (cCREs) and cCRE-to-gene links, together with hypoxia time-course data from human umbilical vein endothelial cells. The manuscript also describes short-read genotyping and quantitative traits for association analyses. Separate accessions for the hypoxia assays and phenotype measurements are unspecified in the manuscript. [ENCODE](https://www.encodeproject.org/)

## Script directory

| Analysis | Directory | Contents |
| --- | --- | --- |
| Genome assembly | [assembly/genome_assembly](scripts/assembly/genome_assembly/) | Read QC, parental k-mer sets, hifiasm trio assembly, assembly statistics, and nextPolish2 polishing. |
| Assembly quality control | [assembly/qc](scripts/assembly/qc/) | Yak QV and trio evaluation, Inspector, QUAST, and QC summaries. |
| Pangenome graph | [pangenome/minigraph_cactus](scripts/pangenome/minigraph_cactus/) | Chromosome-level Minigraph-Cactus graph construction and graph VCF processing with vcfbub and vcfwave. |
| Graph-based genotyping | [pangenome/pangenie](scripts/pangenome/pangenie/) | PanGenie panel and index preparation, per-sample genotyping, and VCF merging. |
| Structural variant discovery | [structural_variants/sv_calling](scripts/structural_variants/sv_calling/) | Read-based SV calling with Sniffles2 and pbsv, plus assembly-based calling with SVIM-asm and PAV. |
| Selection simulations | [selection/pbs_simulation](scripts/selection/pbs_simulation/) | Neutral coalescent simulations and empirical population branch statistic (PBS) P values matched by local mutation and recombination rates. |
| Archaic segment detection | [introgression/segments](scripts/introgression/segments/) | SPrime, HMM and hmmix scripts, a consensus overlap filter, and the TRACE Slurm dispatcher. |
| Archaic SV genotyping | [introgression/archaic_sv_genotyping](scripts/introgression/archaic_sv_genotyping/) | CHM13-based input preparation and Paragraph genotyping for archaic genomes. |
| Introgression statistics | [introgression/fd](scripts/introgression/fd/) | Inputs, window and region estimates, validation, and plots for archaic ancestry using `f_D`. |
| SV haplotype heatmaps | [introgression/heatmap](scripts/introgression/heatmap/) | Slurm and shell drivers for feasibility checks, summaries, and heatmaps. |
| Haplotype homozygosity | [introgression/ehh](scripts/introgression/ehh/) | Linkage disequilibrium (LD) proxies, phase clustering, and extended haplotype homozygosity (EHH) plots. |
| Locus plots | [introgression/locuszoom_daf](scripts/introgression/locuszoom_daf/) | Locus views that combine PBS, archaic `f_D`, derived allele frequencies, and regulatory tracks. |
| Genealogy and allele frequency history | [introgression/relate](scripts/introgression/relate/) | Relate local trees, conversion for CLUES2, and allele frequency trajectories. |
| Local coalescence time | [introgression/msmc2](scripts/introgression/msmc2/) | Coordinate conversion, callable masks, MSMC2 fitting and decoding, and local TMRCA plots. |
| TESV functional annotation | [functional_annotation/tesv](scripts/functional_annotation/tesv/) | cCRE annotation, tissue activity and target-gene links, HUVEC population and time models, and matched-SV enrichment. |

## Running the scripts

Start with the module README where provided. Install the listed tools and set the reference assembly, input paths, and sample tables before running jobs. The release workflows include module-specific configurations or usage notes. Check the expected coordinate system for each input, especially when combining CHM13, GRCh38, and GRCh37 files.

Several downstream workflows retain their original Linux and Slurm paths. Set up the paths and scheduler options for your environment. The `fd` and `relate` drivers require a `CONFIG.env` file. The EHH and MSMC2 Slurm wrappers source separate project configurations. The TRACE dispatcher also expects companion SINGER and TRACE jobs, and the SV heatmap driver expects `build_manifest.R` and a `modules/` directory. Those companion files are absent from this script collection.
