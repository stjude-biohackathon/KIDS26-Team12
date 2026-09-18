#!/bin/bash
# =============================================================================
# run_pipeline.sh  —  run the pipeline on a plain server (no scheduler)
#   SMOKE=1  bash run_pipeline.sh   # ~10-20 min: 20k probes, 3 small cancers
#   SUBSET=1 bash run_pipeline.sh   # prototype: 6 non-CNS cancers, all probes
#   bash run_pipeline.sh            # full development cohort
#   CONTROL=permute_within_cancer SUBSET=1 bash run_pipeline.sh   # tissue-identity control
#   CONTROL=permute SUBSET=1 bash run_pipeline.sh                 # pure label-permutation null
# Progress: tail -f logs/<name>.log     Re-run the same command to resume.
# For the St. Jude HPC (LSF) use submit_hpc.sh instead.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
PROCESSED=${PROCESSED:-/home/KIDS26/DATA/kids26_team12_share/data/processed}
CORES=${CORES:-8}
SUBSET_CANCERS=${SUBSET_CANCERS:-BRCA,UCEC,BLCA,STAD,LUSC,HNSC}

if [[ "${SMOKE:-0}" == "1" ]]; then
  NAME=smoke;  PREP_EXTRA="--max_rows=20000 --top_k=3000 --cancers=ACC,CHOL,KICH,UVM"; FIT_EXTRA="--max_features=500 --alphas=0.5,1"
elif [[ "${SUBSET:-0}" == "1" ]]; then
  NAME=subset; PREP_EXTRA="--top_k=50000 --cancers=$SUBSET_CANCERS"; FIT_EXTRA=""
else
  NAME=full;   PREP_EXTRA="--top_k=50000"; FIT_EXTRA=""
fi
CONTROL=${CONTROL:-none}; INNER=${INNER:-cancer}
FEAT=features_$NAME; OUT=results/${NAME}_${CONTROL}
FIT_EXTRA="$FIT_EXTRA --control=$CONTROL --inner=$INNER"
mkdir -p logs

nohup bash -c "
  set -e
  if [[ ! -f $FEAT/X_dev.rds ]]; then
    Rscript 01_prepare_features.R --beta=$PROCESSED/beta.tsv --samples=$PROCESSED/master_samples.tsv \
      --provenance=$PROCESSED/beta.provenance.json --out=$FEAT $PREP_EXTRA
  fi
  Rscript 02_train_loco.R --features=$FEAT --out=$OUT --cores=$CORES $FIT_EXTRA
" > logs/${NAME}_${CONTROL}.log 2>&1 &
echo "Started '$NAME' (PID $!). Watch with:  tail -f logs/${NAME}_${CONTROL}.log"
