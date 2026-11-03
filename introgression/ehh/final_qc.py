#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True)
a = ap.parse_args()
root = Path(a.root)
checks = []


def add(name, observed, expected, passed):
    checks.append({"check": name, "observed": str(observed), "expected": str(expected),
                   "status": "PASS" if passed else "FAIL"})


manifest = read_rows(root / "metadata/group_manifest.tsv")
proxy = read_rows(root / "results/highest_mappable_polarizable_proxy.tsv")
clusters = read_rows(root / "metadata/phase_clusters.tsv")
summary = read_rows(root / "results/all_groups_EHH_summary.tsv")

add("group_manifest_rows", len(manifest), 64, len(manifest) == 64)
gc = Counter(r["group"] for r in manifest)
add("group_counts_DEN_NEA_GHOST", f"{gc['DEN']},{gc['NEA']},{gc['GHOST']}", "3,20,41",
    (gc["DEN"], gc["NEA"], gc["GHOST"]) == (3, 20, 41))
add("unique_original_SVs", len({r["original_sv_id"] for r in manifest}), 63,
    len({r["original_sv_id"] for r in manifest}) == 63)
add("proxy_rows", len(proxy), 63, len(proxy) == 63)
add("proxy_status_ok", sum(r["status"] == "ok" for r in proxy), 63,
    len(proxy) == 63 and all(r["status"] == "ok" for r in proxy))
add("proxy_AN70", sum(r["lead_snp_an"] == "70" for r in proxy), 63,
    len(proxy) == 63 and all(r["lead_snp_an"] == "70" for r in proxy))
add("proxy_selection_rule_records_AN70", sum("AN=70 required" in r["selection_rule"] for r in proxy), 63,
    len(proxy) == 63 and all("AN=70 required" in r["selection_rule"] for r in proxy))
add("phase_clusters", len(clusters), 44, len(clusters) == 44)
phase_qc = list((root / "results").glob("*.phase_QC.tsv"))
phase_vcf = list((root / "vcf/phased").glob("*.phased.vcf.gz"))
phase_tbi = list((root / "vcf/phased").glob("*.phased.vcf.gz.tbi"))
add("phase_QC_files", len(phase_qc), 44, len(phase_qc) == 44 and all(p.stat().st_size > 0 for p in phase_qc))
add("phased_VCF_files", len(phase_vcf), 44, len(phase_vcf) == 44 and all(p.stat().st_size > 0 for p in phase_vcf))
add("phased_VCF_indexes", len(phase_tbi), 44, len(phase_tbi) == 44 and all(p.stat().st_size > 0 for p in phase_tbi))

add("EHH_summary_rows", len(summary), 64, len(summary) == 64)
sg = Counter(r["group"] for r in summary)
add("EHH_summary_group_counts", f"{sg['DEN']},{sg['NEA']},{sg['GHOST']}", "3,20,41",
    (sg["DEN"], sg["NEA"], sg["GHOST"]) == (3, 20, 41))
status = Counter(r["status"] for r in summary)
add("two_group_EHH_cases", status["OK_TWO_GROUP_AA_DA_EHH"], 63,
    status["OK_TWO_GROUP_AA_DA_EHH"] == 63)
add("single_group_EHH_cases", status["INSUFFICIENT_FOR_TWO_GROUP_EHH_SINGLE_GROUP_MAJOR_ALLELE"], 1,
    status["INSUFFICIENT_FOR_TWO_GROUP_EHH_SINGLE_GROUP_MAJOR_ALLELE"] == 1)
add("samples_per_case_35", sum(r.get("sample_count") == "35" for r in summary), 64,
    len(summary) == 64 and all(r.get("sample_count") == "35" for r in summary))
add("haplotypes_per_case_70", sum(r.get("haplotype_count") == "70" for r in summary), 64,
    len(summary) == 64 and all(r.get("haplotype_count") == "70" for r in summary))

case_dirs = []
for group in ("DEN", "NEA", "GHOST"):
    case_dirs.extend([p for p in (root / group).iterdir() if p.is_dir()])
add("independent_case_directories", len(case_dirs), 64, len(case_dirs) == 64)
for label, pattern in (("case_metadata", "case_metadata.tsv"), ("EHH_source", "*.EHH_source.tsv"),
                       ("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")):
    found = []
    for d in case_dirs:
        found.extend(d.glob(pattern))
    add(f"nonempty_{label}_files", sum(p.stat().st_size > 0 for p in found), 64,
        len(found) == 64 and all(p.stat().st_size > 0 for p in found))

core_ok, bounded_ok = 0, 0
for row in summary:
    source = root / row["group"] / row["case_id"] / f"{row['case_id']}.GRCh38.EHH_source.tsv"
    src = read_rows(source)
    core = [x for x in src if x["distance_from_core_bp"] == "0"]
    distances = [int(x["distance_from_core_bp"]) for x in src]
    if distances and min(distances) >= -300000 and max(distances) <= 300000:
        bounded_ok += 1
    if len(core) == 1:
        if row["status"] == "OK_TWO_GROUP_AA_DA_EHH":
            core_ok += float(core[0]["EHH_AA"]) == 1.0 and float(core[0]["EHH_DA"]) == 1.0
        else:
            core_ok += float(core[0]["EHH_major"]) == 1.0
add("core_EHH_equals_1", core_ok, 64, core_ok == 64)
add("EHH_distances_within_300kb", bounded_ok, 64, bounded_ok == 64)

out = root / "FINAL_QC.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["check", "observed", "expected", "status"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(checks)
failed = [x for x in checks if x["status"] != "PASS"]
print(f"checks={len(checks)} failed={len(failed)}")
if failed:
    for x in failed:
        print(x)
    raise SystemExit(1)

