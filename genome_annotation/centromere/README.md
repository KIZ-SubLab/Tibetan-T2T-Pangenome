# Centromere annotation and introgression

These workflows characterize centromeric satellites and alpha-satellite
higher-order repeats (HORs), compare Tibetan haplotypes with APGp1, and analyze
candidate centromeric introgression.

## Layout

Copy `config.env.example` to `config.env` and set the input, manifest, and software paths.

```
config.env.example
01_Centromere_features/
  01_repeatmasker_annotation.sh      RepeatMasker + satellite-class extraction
  01b_parse_satellite_classes.py     .out -> BED6 satellite classes
  02_alpha_satellite_hor_pipeline.sh monomer -> STV -> k-mer HOR classification
  03_hor_length_and_stats.py         HOR length by haplotype and chromosome, with Kruskal–Wallis and Dunn tests
  04_satellite_composition.py        satellite bp and proportion per centromere
02_Centromeric_introgression/
  01_build_archaic_kmer_databases.sh 5 archaic DBs: count/merge/filter/difference
  02_build_afr_kmer_database.sh      AFR joint DB + count>=10 filter
  03_arc_kmer_scan.sh                ARCkmerFinder scan + >=11 + HOR intersect
  04_introgression_fraction.py       per-haplotype autosomal fraction
  05_snv_projection.sh               hg38 -> haplotype projection (asm5)
  06_sequence_comparison.sh          20-kb/1-kb windows + identity filter
  06b_mannwhitney.py                 one-sided Mann-Whitney U
  07_phylogeny_pipeline.sh           ortholog -> full-flank MSA -> 50-kb -> MAFFT
  07b_project_windows.py             column-block projection to taxon sequences
  08_beast_analysis.py               BEAST2 XML template, execution, and TreeAnnotator
```

## Software

| tool | version used |
|---|---|
| RepeatMasker | 4.1.6 |
| meryl | 1.4.1 (canonical k-mers, k=15, by default) |
| minimap2 | `-x asm5` preset |
| MAFFT | v7.525 |
| BEAST2 | 2.7.7 |
| TreeAnnotator | 2.7.7 |
| bedtools | any recent version |

## Key parameters

- Archaic k-mer filter: `meryl greater-than 2` = **total count >= 3**
- African k-mer filter: `meryl greater-than 9` = **total count >= 10**
- Merged archaic-specific database: `union-sum` (counts summed)
- Final merged distinct k-mers: 156,790,652
- ARCkmerFinder windows: 2,000 bp, non-overlapping; column 4 = **number of
  k-mer occurrences** matching the archaic-specific database
- Significance: **>= 11** k-mers per 2-kb window, uniform across haplotypes
- Candidate window: >=11 AND >=1 bp overlap with that haplotype's merged HOR intervals
- SNV projection: same-chromosome, MAPQ >=30, query coverage >=0.85, best hit per query
- Sequence comparison: 20-kb windows, 1-kb step, MAPQ >=5, aligned length >=5 kb,
  identity >=85%, highest-identity alignment, one-sided Mann-Whitney U
- Phylogeny: 50-kb windows projected from the full-flank MSA, then **re-aligned
  per window** with MAFFT (`--retree 2 --maxiterate 0 --inputorder`)
- BEAST: HKY + 5 gamma categories + pInv **fixed at 0**; strict clock;
  calibrated Yule; one human-chimp calibration; no gorilla;
  1 chain x 50,000,000 states; log/tree every 1,000; 10% burn-in

## Reference data and manifests

- **APGp1 collaborator**: alpha-satellite monomer annotation (`cenanno.bed`),
  HOR classification (`HORclass.bed`), and the `mon2stv.py` /
  HOR k-mer classifier programs.
- **Public data**: archaic SRA/ENA runs, 1KGP African CRAMs, hg38 and the
  haplotype assemblies.
- Input manifests: `ARCHAIC_MANIFEST`, `AFR_SAMPLE_LIST`,
  `HOR_MANIFEST`, `SNV_SEGMENT_MANIFEST`, `TAXON_MANIFEST`,
  `WINDOW_MANIFEST`, `SEQ_COMPARISON_MANIFEST`.

## Statistical analyses

The HOR statistics workflow groups haplotypes by chromosome. It uses a
Kruskal–Wallis test followed by Dunn pairwise comparisons with
Benjamini–Hochberg adjustment. The historical production analysis used
per-chromosome Welch two-sample tests.

The phylogeny workflow builds alignments from centromere flanks.
`08_beast_analysis.py` prepares a BEAST2 XML template from an aligned FASTA.
Use the original BEAST2 run XML when reproducing the production analysis.

## APGp1 comparison set: 152 vs 160

Verified directly from `7.all_APG_stv/0.list`:

| quantity | value |
|---|---|
| APGp1 haplotypes in the list | **320** |
| APGp1 individuals (320 / 2, Mat+Pat) | **160** |
| Tibetan haplotypes present in that list (C101-C108, Mat+Pat) | **16** (= 8 individuals) |
| haplotypes after removing the 8 Tibetan individuals | **304** |
| comparison individuals (304 / 2) | **152** |

So both numbers count **individuals**, each represented by two haplotype
assemblies: 160 is the full APGp1 set, and 152 is that set after removing
the 8 Tibetan individuals that were already analysed separately.
The selection logic is a plain list subtraction of the 16 Tibetan
haplotype IDs; see `01_Centromere_features/APGp1_comparison_selection.md`.

## Archaic k-mer workflow

| step | script | loop object |
|---|---|---|
| archaic per-run counting (5 sources) | `02_Centromeric_introgression/01_build_archaic_kmer_databases.sh` | **archaic libraries** |
| archaic union-sum per source | same | archaic sources |
| archaic filter `greater-than 2` (= count >= 3) | same | archaic sources |
| African joint database + `greater-than 9` (= count >= 10) | `02_.../02_build_afr_kmer_database.sh` | **661 African genomes** |
| African subtraction (`meryl difference`) | `01_build_archaic_kmer_databases.sh` | archaic sources |
| final archaic `union-sum` | `01_build_archaic_kmer_databases.sh` | archaic sources |
| ARCkmerFinder whole-genome scan | `02_.../03_arc_kmer_scan.sh` | **16 Tibetan haplotypes** |

Archaic libraries are used to build the k-mer databases. Tibetan haplotypes are the inputs to the genome scans.
