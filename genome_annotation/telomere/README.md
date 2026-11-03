# Telomere annotation and variation

These workflows detect telomeres in phased assemblies and characterize telomere
architecture, telomeric variant repeats (TVRs), and subtelomeric TAR1 elements.
Use `config.example.sh` to set the input and output paths for your analyses.

```
01_telomere_detection.sh   Part I   detection + extraction
02_telomere_summary.py     Part I   recovery / length statistics
03_CFB_classification.py   Part I   C / F / B architecture classification
04_Btype_variants.py       Part I   B-type motif and substitution spectrum
05_TVR_calling.py          Part II  orientation-normalised TVR caller + blocks
06_TVR_summary.py          Part II  TVR load, variant-length classes, substitution/indel
07_TVR_population_comparison.py     arm-arm Jaccard + Fisher/BH composition test
08_TAR1_annotation.py               VGP-boundary TAR1 10-kb annotation
```

---

## Part I — telomere annotation and architecture

```
assembly FASTA
  -> 01 detection (VGP FindTelomereWindows, bedtools merge, seqkit subseq)
  -> 02 recovery / length
  -> 03 C/F/B classification
  -> 04 B-type motif spectrum
```

**Detection.** The 1-kbp window / 200-bp step / canonical TTAGGG ≥10 % rule is
implemented inside the VGP `telomere.jar`; `01_telomere_detection.sh` wraps the
official VGP driver so detection parameters stay identical to the published
run (identity 99.9, repeat-density threshold 0.4, terminal 15 kb,
`bedtools merge -d 100`).

**Classification (03).**
- B-type: total motifs ≥ 2000/6 and some non-canonical motif > 5 % of motifs.
- Remaining telomeres are C/F candidates.
- C-type if `max(count(TTAGGG), count(CCCTAA)) * 6 / length > 0.75`, else F-type.
- The 0.75 cut-off is the fixed value chosen after KDE / k-means / GMM
  exploration; the clustering itself is not re-run.

```bash
bash 01_telomere_detection.sh --fasta-dir "$ASSEMBLY_DIR" \
    --sample-list "$SAMPLE_LIST_ALL" --output-dir "$TELOMERE_WORK"

python 02_telomere_summary.py --bed-dir "${TELOMERE_WORK}/telo_bed" \
    --sample-list "$SAMPLE_LIST_ALL" --expected-ends "$EXPECTED_ENDS_PER_HAPLOID" \
    --output telomere_summary.tsv

python 03_CFB_classification.py --fasta-dir "${TELOMERE_WORK}/telo_fasta" \
    --sample-list "$SAMPLE_LIST_ALL" --output telomere_types.tsv

python 04_Btype_variants.py --input telomere_types.tsv --type B \
    --motif-output Btype_motifs.tsv --substitution-output Btype_substitutions.tsv
```

## Part II — TVR and TAR1

```
telomeric FASTA
  -> 05 orientation normalisation + TVR calling + TVR blocks
  -> 06 per-arm TVR load / variant-length classes / substitution vs indel
  -> 07 arm-arm Jaccard + between-group composition test
VGP telomere BED + RepeatMasker TAR1 -> 08 10-kb boundary-proximal annotation
```

**Orientation (05).** p-arm is reverse-complemented, q-arm keeps assembly
orientation, so every arm is reported as `0 bp = subtelomere–telomere boundary`,
increasing towards the distal telomere. Calling scans 2 strands × 6 frame
offsets and labels each 6-bp slot canonical / rotated / ambiguous / tvr, choosing
the frame with `score = 10*n_canonical - 2*n_rotated - 3*n_ambiguous - 0.2*n_tvr`.
Consecutive TVR slots become blocks; blocks are kept with 6–300 bp.

**TVR load (06).** `sum(block length, 6–300 bp) / telomeric interval length`.
Variant-length classes come from the Hamming distance of each TVR 6-mer to its
strand canonical (1 bp / 2 bp / 3–4 bp / 5–6 bp). Indel-associated slots are
those inside a block whose `length != n_hexamers*6`, or slots not contiguous
with the previous TVR slot; no alignment is used.

**Composition test (07).** For each arm × hexamer, two-sided Fisher exact test on
haplotype-level presence/absence, followed by Benjamini–Hochberg FDR. The unit
is "a haplotype carries the hexamer", never the copy number.

```bash
python 05_TVR_calling.py --input "${TELO_FASTA_DIR}/sample.telo.fasta" \
    --hexamer-tsv "${TVR_WORK}/sample_hexamers.tsv" \
    --blocks-tsv  "${TVR_WORK}/sample_tvr_blocks.tsv" \
    --min-tvr "$TVR_MIN_BLOCK" --max-tvr "$TVR_MAX_BLOCK"

python 06_TVR_summary.py --hexamer-dir "$TVR_WORK" --blocks-dir "$TVR_WORK" \
    --sample-list "$SAMPLE_LIST_A" --label A --output tvr_summary_A.tsv

python 07_TVR_population_comparison.py \
    --group-a-dir "$TVR_WORK/a" --group-a-samples "$SAMPLE_LIST_A" --group-a-label A \
    --group-b-dir "$TVR_WORK/b" --group-b-samples "$SAMPLE_LIST_B" --group-b-label B \
    --output-dir comparison/
```

**TAR1 (08).** Boundary = proximal end of the p-arm telomeric array
(`max(end)` of VGP windows) or proximal start of the q-arm array (`min(start)`).
TAR1 is called present when a RepeatMasker TAR1 feature overlaps
`[boundary, boundary+10 kb]` (p) or `[boundary−10 kb, boundary]` (q) — interior
direction only. RepeatMasker output produced on extracted windows is mapped back
to genomic coordinates with `--gff-coordinate-mode windowed` plus the assembly
and extract `.fai` directories; no group-specific code exists. Haplotypes
without a boundary or annotation are excluded from the denominator.
Arm is TAR1+ when positives ≥ `ceil(n_effective/2)`.

```bash
python 08_TAR1_annotation.py --vpg-bed-dir "$VPG_BED_DIR" --gff-dir "$RM_GFF_DIR" \
    --sample-list "$SAMPLE_LIST_A" --group-label A --output-dir tar1_out

# windowed RepeatMasker output
python 08_TAR1_annotation.py --vpg-bed-dir "$VPG_BED_DIR" --gff-dir "$RM_GFF_DIR" \
    --sample-list "$SAMPLE_LIST_B" --group-label B --output-dir tar1_out \
    --gff-coordinate-mode windowed \
    --assembly-fai-dir "$ASSEMBLY_FAI_DIR" --extract-fai-dir "$EXTRACT_FAI_DIR"
```

---

## Requirements

python ≥3.8 (numpy, scipy, biopython optional), bedtools, seqkit, samtools, java,
VGP assembly pipeline (`telomere/telomere_analysis.sh`, `telomere/telomere.jar`).
RepeatMasker is used upstream to produce the TAR1 GFF and is not re-run here.
