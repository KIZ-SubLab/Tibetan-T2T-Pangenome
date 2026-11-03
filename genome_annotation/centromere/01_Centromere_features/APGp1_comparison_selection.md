# APGp1 comparison set selection (152 vs 160)

Evidence: `2.APT_annotation/7.all_APG_stv/0.list` (320 haplotype IDs).

```
320 haplotypes                     = 160 individuals x (Mat + Pat)
16 Tibetan haplotypes in that list = C101-C108, Mat + Pat = 8 individuals
320 - 16 = 304 haplotypes          = 152 individuals
```

* `160` = the full APGp1 individual set (as 320 haplotype assemblies).
* `152` = the comparison set after removing the 8 Tibetan individuals that
  were analysed separately.

Selection rule (list subtraction, no other filtering applied):

```bash
# TIBETAN_ID_PATTERN matches the 8 Tibetan individuals (C101..C108)
grep -v -E "$TIBETAN_ID_PATTERN" 0.list > comparison_152.haplotypes
# 304 lines -> 152 individuals
```

Both counts refer to **individuals**, each represented by two haplotype
assemblies (Mat and Pat). This is a list-level subtraction; no additional
sample QC or population filter was found in the production scripts.
