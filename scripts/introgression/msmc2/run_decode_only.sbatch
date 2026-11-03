#!/usr/bin/env bash
#SBATCH --job-name=TBN2_decode_only
#SBATCH --account=acc2159rgq
#SBATCH --partition=wzhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/TMRCA/logs/%x.%j.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/TMRCA/logs/%x.%j.err

set -euo pipefail
BASE=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/TMRCA
source "$BASE/config/config.sh"
bash "$BASE/scripts/04b_decode_only.sh"
