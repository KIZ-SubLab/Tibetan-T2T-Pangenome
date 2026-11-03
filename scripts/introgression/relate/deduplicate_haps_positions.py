#!/usr/bin/env python3
"""Retain one deterministic biallelic record per physical coordinate."""

import argparse
from itertools import groupby


p = argparse.ArgumentParser()
p.add_argument("--input", required=True)
p.add_argument("--output", required=True)
p.add_argument("--audit", required=True)
a = p.parse_args()


def parsed(line):
    fields = line.split()
    if len(fields) < 7:
        raise ValueError(f"malformed haps line: {line[:120]!r}")
    position = int(fields[2])
    states = [int(x) for x in fields[5:]]
    if any(x not in (0, 1) for x in states):
        raise ValueError(f"non-binary haplotype state at {position}")
    allele_count = sum(states)
    minor_count = min(allele_count, len(states) - allele_count)
    return position, minor_count, fields


with open(a.input, encoding="utf-8") as source, open(
    a.output, "w", encoding="utf-8"
) as output, open(a.audit, "w", encoding="utf-8") as audit:
    audit.write(
        "position\tkept_id\tkept_ancestral\tkept_derived\tkept_MAC\t"
        "dropped_id\tdropped_ancestral\tdropped_derived\tdropped_MAC\n"
    )
    previous_position = -1
    for position, rows_iter in groupby(source, key=lambda line: int(line.split()[2])):
        if position < previous_position:
            raise ValueError(f"input is not sorted: {position} after {previous_position}")
        previous_position = position
        rows = [(line, *parsed(line)[1:]) for line in rows_iter]
        rows.sort(key=lambda item: (-item[1], item[2][1], item[2][3], item[2][4]))
        kept_line, kept_mac, kept = rows[0]
        output.write(kept_line)
        for _, dropped_mac, dropped in rows[1:]:
            audit.write(
                f"{position}\t{kept[1]}\t{kept[3]}\t{kept[4]}\t{kept_mac}\t"
                f"{dropped[1]}\t{dropped[3]}\t{dropped[4]}\t{dropped_mac}\n"
            )

