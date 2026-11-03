#!/usr/bin/env python3
"""One-command PBS + archaic fD + regulatory locus figure workflow (Python only)."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


WORK_ROOT = Path("/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/locuszoom_cms_fd_daf")
FD_BASE = Path("/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/EPAS1/fd_analysis")
GHOST_FD_BASE = Path("/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/Ghost_DISC1/fd_analysis")
FD_ENV = Path("/work/home/songjiaxu/anaconda3/envs/fd_trace")
PBS_DEFAULT = WORK_ROOT / "snp.all.sites.pbs.PBS_gt0.tsv.gz"
KG_ROOT = Path("/work/share/acc2159rgq/1kGPx30/1KGP_Change_Format/1kGP_SNV")
FAM = Path("/work/home/songjiaxu/jianxin_cms/caculate_results/cms_corr_4_parament/daf_SGPL1_single/relate_time_395/1KGP/1000G_high_coverage.2504.fam")
TBN_VCF = Path("/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/vcf_hg38/TBN35.grch38.biSNV.vcf.gz")
REFERENCE_CACHE = Path("/work/home/songjiaxu/jianxin_cms/caculate_results/cms_corr_4_parament/daf_SGPL1_single/relate_time_395/96/north_locuszoom_cms/_reference_cache")
ARCHAIC_ROOTS = {
    "den": Path("/work/home/songjiaxu/introgression/fastq/den/SNPad_GRCh38/all_chr_compressed"),
    "nea": Path("/work/home/songjiaxu/introgression/fastq/nea/Altai/snpad_GRCh38/all_chr_compressed"),
}
ARCHAIC_META = {
    "den": {"code": "DEN", "label": "Denisovan", "sample": "Denisova"},
    "nea": {"code": "NEA", "label": "Altai Neanderthal", "sample": "Altai"},
    "ghost": {"code": "GHOST", "label": "Ghost candidate", "sample": "not directly sampled"},
}


def run(command: list[str], stdout: Path | None = None, stderr: Path | None = None):
    print("+", " ".join(map(str, command)), flush=True)
    out_handle = stdout.open("w") if stdout else None
    err_handle = stderr.open("w") if stderr else None
    try:
        subprocess.run(command, check=True, stdout=out_handle, stderr=err_handle)
    finally:
        if out_handle:
            out_handle.close()
        if err_handle:
            err_handle.close()


def require_file(path: Path, description: str, indexed: bool = False):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing {description}: {path}")
    if indexed and not any(Path(str(path) + suffix).is_file() for suffix in [".tbi", ".csi"]):
        raise FileNotFoundError(f"Missing index (.tbi or .csi) for {description}: {path}")


def parse_chrom(value: str) -> tuple[str, int]:
    # str.removeprefix() requires Python >= 3.9, while the HPCC login
    # environment may still expose Python 3.8 as `python3`.
    text = value.strip().lower()
    if text.startswith("chr"):
        text = text[3:]
    if not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2])", text):
        raise argparse.ArgumentTypeError("chromosome must be an autosome from 1 to 22")
    number = int(text)
    return f"chr{number}", number


def regional_pbs_max(path: Path, chrom: str, start: int, end: int) -> tuple[int, float]:
    opener = gzip.open if path.suffix == ".gz" else open
    count, maximum = 0, float("-inf")
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["CHR", "POS", "LOG_R", "LOG_MU", "PBS"]
        if reader.fieldnames != expected:
            raise RuntimeError(f"Unexpected PBS header: {reader.fieldnames}; expected {expected}")
        for row in reader:
            if row["CHR"] != chrom:
                continue
            pos = int(row["POS"])
            if start <= pos <= end:
                value = float(row["PBS"])
                count += 1
                maximum = max(maximum, value)
    if count == 0:
        raise RuntimeError(f"No PBS rows in {chrom}:{start}-{end}")
    return count, maximum


def ensure_ancestral_fasta(chrom_no: int) -> Path:
    ancestor_dir = FD_BASE / "inputs/homo_sapiens_ancestor_GRCh38"
    target = ancestor_dir / f"homo_sapiens_ancestor_{chrom_no}.fa"
    if target.is_file() and target.stat().st_size > 0:
        return target
    archive = FD_BASE / "inputs/homo_sapiens_ancestor_GRCh38.tar.gz"
    require_file(archive, "ancestral FASTA archive")
    ancestor_dir.mkdir(parents=True, exist_ok=True)
    member = f"homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom_no}.fa"
    print(f"Extracting {member} from ancestral FASTA archive...", flush=True)
    run(["tar", "-xzf", str(archive), "-C", str(FD_BASE / "inputs"), member])
    require_file(target, f"chr{chrom_no} ancestral FASTA")
    return target


def write_fd_config(path: Path, *, run_root: Path, chrom: str, chrom_no: int,
                    start: int, end: int, archaic: str, ancestor: Path):
    meta = ARCHAIC_META[archaic]
    archaic_vcf = ARCHAIC_ROOTS[archaic] / f"chr{chrom_no}.vcf.gz"
    kg_vcf = KG_ROOT / f"1kGP_high_coverage.chr{chrom_no}.filtered.SNV_phased_panel.vcf.gz"
    bed = run_root / "fd_analysis/region.bed"
    config = run_root / "fd_analysis/config.env"
    bed.write_text(f"{chrom}\t{start - 1}\t{end}\n")
    values = {
        "BASE": FD_BASE,
        "ENV_DIR": FD_ENV,
        "REGION_LABEL": f"{chrom}_{start}_{end}",
        "ANALYSIS_TAG": f"{chrom}_{start}_{end}_{meta['code']}",
        "ARCHAIC_CODE": meta["code"],
        "ARCHAIC_LABEL": meta["label"],
        "ARCHAIC_SAMPLES": meta["sample"],
        "BED": bed,
        "KG_VCF": kg_vcf,
        "FAM": FAM,
        "TBN_VCF": TBN_VCF,
        "ARCHAIC_VCF": archaic_vcf,
        "ANCESTRAL_FASTA": ancestor,
        "OUTPUT_DIR": run_root / "fd_analysis/run",
        "WINDOW_SNPS": 100,
        "STEP_SNPS": 50,
        "WINDOW_JACKKNIFE_BLOCK_SNPS": 10,
        "JACKKNIFE_BP": 50000,
    }
    config.write_text(
        "\n".join(f"{key}={shlex.quote(str(value))}" for key, value in values.items()) + "\n")
    return config, archaic_vcf, kg_vcf


def write_ghost_config(*, run_root: Path, chrom: str, chrom_no: int,
                       start: int, end: int, ancestor: Path):
    den_vcf = ARCHAIC_ROOTS["den"] / f"chr{chrom_no}.vcf.gz"
    nea_vcf = ARCHAIC_ROOTS["nea"] / f"chr{chrom_no}.vcf.gz"
    kg_vcf = KG_ROOT / f"1kGP_high_coverage.chr{chrom_no}.filtered.SNV_phased_panel.vcf.gz"
    bed = run_root / "fd_analysis/region.bed"
    config = run_root / "fd_analysis/config.env"
    bed.write_text(f"{chrom}\t{start - 1}\t{end}\n")
    values = {
        "BASE": GHOST_FD_BASE,
        "ENV_DIR": FD_ENV,
        "EPAS_FD_BASE": FD_BASE,
        "REGION_LABEL": f"{chrom}_{start}_{end}_ghost_candidate",
        "ANALYSIS_TAG": f"{chrom}_{start}_{end}_GHOST",
        "BED": bed,
        "KG_VCF": kg_vcf,
        "FAM": FAM,
        "TBN_VCF": TBN_VCF,
        "DEN_VCF": den_vcf,
        "NEA_VCF": nea_vcf,
        "ANCESTRAL_FASTA": ancestor,
        "OUTPUT_DIR": run_root / "fd_analysis/run",
        "WINDOW_SNPS": 100,
        "STEP_SNPS": 50,
        "WINDOW_JACKKNIFE_BLOCK_SNPS": 10,
        "REGION_JACKKNIFE_BP": 50000,
    }
    config.write_text(
        "\n".join(f"{key}={shlex.quote(str(value))}" for key, value in values.items()) + "\n")
    return config, den_vcf, nea_vcf, kg_vcf


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Run regional DEN/NEA or ghost-candidate A/B/C fD tests and build a PBS/fD/cCRE/ENCODE/GENCODE locus figure.")
    parser.add_argument("--chrom", required=True, help="Autosome, e.g. 2 or chr2")
    parser.add_argument("--start", required=True, type=int, help="GRCh38 1-based inclusive start")
    parser.add_argument("--end", required=True, type=int, help="GRCh38 1-based inclusive end")
    parser.add_argument("--archaic", required=True, choices=["den", "nea", "ghost"], help="Archaic source or ghost-candidate mode")
    parser.add_argument("--pbs-ymax", type=float, help="Manual PBS y-axis maximum; omit for automatic scaling")
    parser.add_argument("--pbs-file", type=Path, default=PBS_DEFAULT)
    parser.add_argument("--point-size", type=float, default=15.0)
    parser.add_argument("--run-name", help="Unique result directory name under the project root")
    args = parser.parse_args()

    chrom, chrom_no = parse_chrom(args.chrom)
    if args.start < 1 or args.end <= args.start:
        parser.error("require 1 <= start < end")
    if args.pbs_ymax is not None and args.pbs_ymax <= 0:
        parser.error("--pbs-ymax must be positive")
    if args.point_size <= 0:
        parser.error("--point-size must be positive")
    meta = ARCHAIC_META[args.archaic]
    run_name = args.run_name or f"{chrom}_{args.start}_{args.end}_{meta['code']}_PBS_fD"
    if not re.fullmatch(r"[A-Za-z0-9._-]+", run_name):
        parser.error("--run-name may contain only letters, numbers, dot, underscore, and hyphen")
    run_root = WORK_ROOT / run_name
    if run_root.exists():
        raise FileExistsError(f"Refusing to overwrite existing run directory: {run_root}")

    require_file(args.pbs_file, "PBS table")
    require_file(FAM, "1KGP population FAM")
    require_file(TBN_VCF, "TBN VCF", indexed=True)
    require_file(FD_BASE / "scripts/run_archaic_fd.sh", "fd workflow entry script")
    if args.archaic == "ghost":
        require_file(GHOST_FD_BASE / "scripts/run_ghost_fd.sh", "ghost A/B/C fd workflow entry script")
    require_file(REFERENCE_CACHE / "gencode.v50.basic.annotation.gtf.gz", "GENCODE GTF")
    require_file(REFERENCE_CACHE / "encodeCcreCombined_ENCODE3.bb", "ENCODE3 cCRE bigBed")
    pbs_count, pbs_max = regional_pbs_max(args.pbs_file, chrom, args.start, args.end)
    if args.pbs_ymax is not None and args.pbs_ymax < pbs_max:
        parser.error(
            f"--pbs-ymax {args.pbs_ymax:g} is below the current regional PBS maximum {pbs_max:g}; "
            "choose a larger value so no source point is clipped")

    ancestor = ensure_ancestral_fasta(chrom_no)
    run_root.mkdir(parents=True)
    for folder in ["data", "figures", "metadata", "qa", "refs", "scripts", "fd_analysis"]:
        (run_root / folder).mkdir()
    (run_root / "refs/gencode.v50.basic.annotation.gtf.gz").symlink_to(
        REFERENCE_CACHE / "gencode.v50.basic.annotation.gtf.gz")
    (run_root / "refs/encodeCcreCombined_ENCODE3.bb").symlink_to(
        REFERENCE_CACHE / "encodeCcreCombined_ENCODE3.bb")

    script_dir = Path(__file__).resolve().parent
    helper_names = ["01_prepare_pbs_fd_tables.py", "02_plot_pbs_fd_locus.py",
                    "03_final_qa.py", "fetch_ucsc_layered_h3k27ac.py", Path(__file__).name]
    for name in helper_names:
        shutil.copy2(script_dir / name, run_root / "scripts" / name)

    if args.archaic == "ghost":
        config, den_vcf, nea_vcf, kg_vcf = write_ghost_config(
            run_root=run_root, chrom=chrom, chrom_no=chrom_no,
            start=args.start, end=args.end, ancestor=ancestor)
        require_file(den_vcf, f"Denisovan chr{chrom_no} VCF", indexed=True)
        require_file(nea_vcf, f"Altai Neanderthal chr{chrom_no} VCF", indexed=True)
        require_file(kg_vcf, f"1KGP chr{chrom_no} VCF", indexed=True)
        run(["bash", str(GHOST_FD_BASE / "scripts/run_ghost_fd.sh"), str(config)])
    else:
        config, archaic_vcf, kg_vcf = write_fd_config(
            run_root / "fd_analysis/config.env", run_root=run_root, chrom=chrom, chrom_no=chrom_no,
            start=args.start, end=args.end, archaic=args.archaic, ancestor=ancestor)
        require_file(archaic_vcf, f"{meta['label']} chr{chrom_no} VCF", indexed=True)
        require_file(kg_vcf, f"1KGP chr{chrom_no} VCF", indexed=True)
        run(["bash", str(FD_BASE / "scripts/run_archaic_fd.sh"), str(config)])
    fd_windows = run_root / "fd_analysis/run/results/windows_100SNP_step50.tsv"
    require_file(fd_windows, "newly generated fD windows")
    run([
        sys.executable, str(run_root / "scripts/01_prepare_pbs_fd_tables.py"), str(run_root),
        "--pbs-file", str(args.pbs_file), "--fd-windows", str(fd_windows),
        "--gtf", str(run_root / "refs/gencode.v50.basic.annotation.gtf.gz"),
        "--archaic-label", meta["label"], "--archaic-code", meta["code"], "--chrom", chrom,
        "--start", str(args.start), "--end", str(args.end),
    ], run_root / "metadata/prepare.stdout.log", run_root / "metadata/prepare.stderr.log")
    run([
        sys.executable, str(run_root / "scripts/fetch_ucsc_layered_h3k27ac.py"), str(run_root),
        "--chrom", chrom, "--start", str(args.start), "--end", str(args.end),
    ], run_root / "metadata/h3k27ac.stdout.log", run_root / "metadata/h3k27ac.stderr.log")
    plot_command = [
        sys.executable, str(run_root / "scripts/02_plot_pbs_fd_locus.py"), str(run_root),
        "--chrom", chrom, "--start", str(args.start), "--end", str(args.end),
        "--archaic-code", meta["code"], "--archaic-label", meta["label"],
        "--point-size", str(args.point_size),
    ]
    qa_command = [
        sys.executable, str(run_root / "scripts/03_final_qa.py"), str(run_root),
        "--chrom", chrom, "--start", str(args.start), "--end", str(args.end),
        "--archaic-code", meta["code"], "--archaic-label", meta["label"],
    ]
    if args.pbs_ymax is not None:
        plot_command += ["--pbs-ymax", str(args.pbs_ymax)]
        qa_command += ["--pbs-ymax", str(args.pbs_ymax)]
    run(plot_command, run_root / "metadata/plot.stdout.log", run_root / "metadata/plot.stderr.log")
    run(qa_command, run_root / "qa/qa.stdout.log", run_root / "qa/qa.stderr.log")

    provenance = {
        "command": " ".join(sys.argv),
        "region": f"{chrom}:{args.start}-{args.end}",
        "archaic": meta,
        "ghost_interpretation": (
            "A/B/C sampled-archaic affinity tests; no directly sampled ghost fD"
            if args.archaic == "ghost" else None
        ),
        "pbs_rows_preflight": pbs_count,
        "pbs_max_preflight": pbs_max,
        "pbs_ymax_requested": args.pbs_ymax,
        "inputs": {str(path): sha256(path) for path in [args.pbs_file, config, run_root / "fd_analysis/region.bed"]},
    }
    (run_root / "metadata/integrated_run.json").write_text(json.dumps(provenance, indent=2) + "\n")
    y_text = str(args.pbs_ymax) if args.pbs_ymax is not None else "automatic"
    (run_root / "README.zh-CN.md").write_text(
        f"# {chrom}:{args.start:,}-{args.end:,} PBS–{meta['label']} fD 整合图\n\n"
        f"本目录由一体化 Python 入口生成。古人类来源：{meta['label']}；PBS y 轴上限：{y_text}。\n\n"
        "流程依次运行区域 fd 分析、PBS/fd 表整理、cCRE/H3K27ac/GENCODE 注释、四格式绘图和自动 QA。\n")
    print(f"DONE: {run_root}")


if __name__ == "__main__":
    main()
