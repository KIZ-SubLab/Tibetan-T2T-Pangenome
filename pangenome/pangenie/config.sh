#!/usr/bin/env bash
# Edit the paths below before running the pipeline.
# Values can also be overridden from the environment.

PANGENIE_SIF=${PANGENIE_SIF:-/path/to/pangenie.sif}
PANGENIE_REPO=${PANGENIE_REPO:-/path/to/pangenie}
REFERENCE=${REFERENCE:-/path/to/reference.fa}
PANEL_VCF=${PANEL_VCF:-/path/to/panel.phased.vcf.gz}
WORKDIR=${WORKDIR:-/path/to/pangenie_run}
SAMPLE_SHEET=${SAMPLE_SHEET:-/path/to/samples.tsv}

THREADS=${THREADS:-16}
SINGULARITY=${SINGULARITY:-singularity}

# Keep tool temporary files inside the working directory.
TMPDIR="${WORKDIR}/tmp"
export TMPDIR
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"
mkdir -p "${TMPDIR}"

# Run PanGenie from the Singularity/Apptainer image.
# If PanGenie is installed natively, replace this with:
#   pangenie_exec() { "$@"; }
pangenie_exec() {
    "${SINGULARITY}" exec "${PANGENIE_SIF}" "$@"
}
