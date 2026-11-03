#!/usr/bin/env bash
# Clean the panel VCF for the PanGenie preparation pipeline:
#   - remove symbolic records other than <DEL>
#   - resolve symbolic <DEL> alleles against the reference
#   - convert ./. to .|.
#   - sort and check REF alleles
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${HERE}/config.sh"

IN=${1:-${PANEL_VCF}}
OUTDIR="${WORKDIR}/input"
OUT="${OUTDIR}/panel.sorted.vcf"
mkdir -p "${OUTDIR}"

command -v bcftools >/dev/null
command -v samtools >/dev/null

TMP=$(mktemp -d -p "${OUTDIR}" .tmp.XXXXXX)
trap 'rm -rf "${TMP}"' EXIT

echo "[1/4] removing symbolic records other than <DEL>"
bcftools view -e 'ALT~"<" && ALT!="<DEL>"' -Ov \
    -o "${TMP}/no_symbolic.vcf" "${IN}"

echo "[2/4] resolving symbolic DEL alleles"
python3 - "${TMP}/no_symbolic.vcf" "${REFERENCE}" "${TMP}/seq_resolved.vcf" <<'PY'
import subprocess
import sys

in_vcf, ref_fa, out_vcf = sys.argv[1:4]


def get_ref(chrom, start, end):
    res = subprocess.check_output(
        ["samtools", "faidx", ref_fa, f"{chrom}:{start}-{end}"],
        stderr=subprocess.DEVNULL,
    )
    return "".join(res.decode().splitlines()[1:]).upper()


with open(in_vcf) as fin, open(out_vcf, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[4] == "<DEL>":
            end = None
            for item in fields[7].split(";"):
                if item.startswith("END="):
                    end = int(item[4:])
                    break
            pos = int(fields[1])
            if end is not None and end > pos:
                ref = get_ref(fields[0], pos, end)
                fields[3] = ref
                fields[4] = ref[0]
        fout.write("\t".join(fields) + "\n")
PY

echo "[3/4] converting ./. to .|."
sed 's/\t\.\/\./\t.|./g' "${TMP}/seq_resolved.vcf" > "${TMP}/phased.vcf"

echo "[4/4] sorting and checking REF alleles"
bcftools sort -Ov -o "${OUT}" "${TMP}/phased.vcf"
bcftools norm --check-ref e --fasta-ref "${REFERENCE}" "${OUT}" \
    >/dev/null 2> "${OUTDIR}/check_ref.log" || true

echo "output: ${OUT}"
echo "records: $(grep -vc '^#' "${OUT}")"
