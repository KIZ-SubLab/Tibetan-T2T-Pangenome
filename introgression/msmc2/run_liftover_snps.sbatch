#!/usr/bin/env bash
#SBATCH --job-name=T2T012_snp_lift
#SBATCH --account=acc2159rgq
#SBATCH --partition=wzhcnormal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/TMRCA/logs/%x.%j.out
#SBATCH --error=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/TMRCA/logs/%x.%j.err

set -euo pipefail
BASE=/work/home/songjiaxu/introgression/data_HMM_intro/GRCh38/output/rm1kb-na-segments/merged_fragments_v3/PBS_sv-snp/TMRCA_dating/TMRCA
source "$BASE/config/config.sh"
bash "$BASE/scripts/01_liftover_snps_gatk.sh"
