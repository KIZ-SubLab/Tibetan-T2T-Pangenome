#!/usr/bin/env python3
"""Prepare sequence-resolved CHM13 INS/DEL records for Paragraph."""

import argparse
import collections
import re
import sys

import pysam


DNA = re.compile(r"^[ACGTNacgtn]+$")


def scalar(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument(
        "--keep-filtered",
        action="store_true",
        help="Keep records whose FILTER is not PASS or '.'.",
    )
    args = parser.parse_args()

    reference = pysam.FastaFile(args.reference)
    source = pysam.VariantFile(args.input)
    header = source.header.copy()
    header.add_meta(
        "paragraph_prepare",
        value=(
            "INS/DEL only; DEL REF rebuilt from CHM13 using SVLEN; "
            "INS anchor replaced with CHM13; malformed records removed"
        ),
    )
    output = pysam.VariantFile(args.output, "wz", header=header)

    counts = collections.Counter()
    corrected = collections.Counter()
    seen = set()

    for record in source:
        counts["input"] += 1
        svtype = scalar(record.info.get("SVTYPE"))
        if svtype not in {"INS", "DEL"}:
            counts["skip_non_ins_del"] += 1
            continue

        filters = list(record.filter.keys())
        if not args.keep_filtered and filters and filters != ["PASS"]:
            counts["skip_filter_not_pass"] += 1
            continue

        if len(record.alleles or ()) != 2:
            counts["skip_not_biallelic"] += 1
            continue
        old_ref, old_alt = record.alleles
        if not old_ref or not old_alt or not DNA.fullmatch(old_ref) or not DNA.fullmatch(old_alt):
            counts["skip_non_sequence_allele"] += 1
            continue
        if record.contig not in reference.references:
            counts["skip_contig_absent"] += 1
            continue

        raw_svlen = scalar(record.info.get("SVLEN"))
        try:
            svlen = int(raw_svlen)
        except (TypeError, ValueError):
            svlen = len(old_alt) - len(old_ref)

        if svtype == "INS":
            if len(old_ref) != 1 or len(old_alt) <= 1:
                counts["skip_malformed_ins"] += 1
                continue
            inserted = old_alt[1:]
            if svlen <= 0 or len(inserted) != svlen:
                counts["skip_ins_svlen_mismatch"] += 1
                continue
            anchor = reference.fetch(record.contig, record.pos - 1, record.pos).upper()
            if len(anchor) != 1 or anchor == "N":
                counts["skip_unusable_reference"] += 1
                continue
            new_ref = anchor
            new_alt = anchor + inserted.upper()
            if old_ref.upper() != new_ref or old_alt[0].upper() != anchor:
                corrected["corrected_ins_anchor"] += 1
            record.alleles = (new_ref, new_alt)
            # For a sequence-resolved insertion, REF spans the anchor only.
            record.stop = record.pos

        else:
            deletion_length = abs(svlen)
            if svlen >= 0 or deletion_length < 1:
                counts["skip_invalid_del_svlen"] += 1
                continue
            if len(old_alt) != 1 or len(old_ref) != deletion_length + 1:
                counts["skip_del_svlen_mismatch"] += 1
                continue
            start0 = record.pos - 1
            end0 = start0 + deletion_length + 1
            if end0 > reference.get_reference_length(record.contig):
                counts["skip_del_out_of_bounds"] += 1
                continue
            new_ref = reference.fetch(record.contig, start0, end0).upper()
            if len(new_ref) != deletion_length + 1 or "N" in new_ref:
                counts["skip_unusable_reference"] += 1
                continue
            new_alt = new_ref[0]
            if old_ref.upper() != new_ref or old_alt.upper() != new_alt:
                corrected["corrected_del_sequence"] += 1
            record.alleles = (new_ref, new_alt)
            record.stop = record.pos + deletion_length

        key = (record.contig, record.pos, record.ref, record.alts[0], svtype)
        if key in seen:
            counts["skip_exact_duplicate"] += 1
            continue
        seen.add(key)

        # Final invariants required by Paragraph.
        genome_anchor = reference.fetch(record.contig, record.pos - 1, record.pos).upper()
        if record.ref[0].upper() != genome_anchor or record.alts[0][0].upper() != genome_anchor:
            raise RuntimeError(f"Internal anchor validation failed at {record.contig}:{record.pos}")

        output.write(record)
        counts[f"written_{svtype.lower()}"] += 1
        counts["written_total"] += 1

    output.close()
    source.close()
    reference.close()

    with open(args.stats, "w", encoding="utf-8") as report:
        report.write("metric\tcount\n")
        for key in sorted(counts):
            report.write(f"{key}\t{counts[key]}\n")
        for key in sorted(corrected):
            report.write(f"{key}\t{corrected[key]}\n")

    print(
        f"Wrote {counts['written_total']} records "
        f"(INS={counts['written_ins']}, DEL={counts['written_del']}).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
