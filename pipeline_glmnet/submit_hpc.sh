#!/bin/bash
# =============================================================================
# submit_hpc.sh  —  submit the pipeline to the St. Jude HPC (LSF)
#   MODE=smoke  bash submit_hpc.sh   # quick test
#   MODE=subset bash submit_hpc.sh   # prototype on 6 non-CNS cancers
#   MODE=full   bash submit_hpc.sh   # full development cohort (overnight)
# Two jobs: (1) feature prep, (2) LOCO training that starts when (1) succeeds.
# Check before first use:  module avail R   |   bqueues   |   ls $PROCESSED
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
WORK=$(pwd)

# ---------------- CONFIG (edit for your account) ----------------------------------
PROCESSED=${PROCESSED:-/lustre_scratch/shared_scratch/kids26_team12_share/data/processed}
R_MODULE=${R_MODULE:-R/4.5.0}
QUEUE=${QUEUE:-standard}
CORES=${CORES:-8}
SUBSET_CANCERS=${SUBSET_CANCERS:-BRCA,UCEC,BLCA,STAD,LUSC,HNSC}
MODE=${MODE:-subset}
CONTROL=${CONTROL:-none}     # none | permute | permute_within_cancer
INNER=${INNER:-cancer}       # cancer (nested LOCO) | random
# ----------------------------------------------------------------------------------

case "$MODE" in
  smoke)  PREP_EXTRA="--max_rows=20000 --top_k=3000 --cancers=ACC,CHOL,KICH,UVM"; FIT_EXTRA="--max_features=500 --alphas=0.5,1"
          PREP_MEM=8000;  PREP_W=1:00; FIT_MEM=16000; FIT_W=1:00 ;;
  subset) PREP_EXTRA="--top_k=50000 --cancers=$SUBSET_CANCERS"; FIT_EXTRA=""
          PREP_MEM=24000; PREP_W=6:00; FIT_MEM=48000; FIT_W=8:00 ;;
  full)   PREP_EXTRA="--top_k=50000"; FIT_EXTRA=""
          PREP_MEM=48000; PREP_W=8:00; FIT_MEM=64000; FIT_W=16:00 ;;
  *) echo "MODE must be smoke, subset or full"; exit 1 ;;
esac
FEAT=$WORK/features_$MODE; OUT=$WORK/results/${MODE}_${CONTROL}; TAG="hrd_${MODE}_${CONTROL}_$(date +%m%d%H%M)"
FIT_EXTRA="$FIT_EXTRA --control=$CONTROL --inner=$INNER"
mkdir -p "$WORK/logs"
LOAD="module load $R_MODULE"
# LSF rusage[mem] is PER SLOT at St. Jude and multiplies by -n (docs/23), so the
# per-slot request below is the total divided by CORES.
FIT_MEM_SLOT=$(( FIT_MEM / CORES ))

if [[ ! -f "$FEAT/X_dev.rds" ]]; then
  bsub -q "$QUEUE" -J "${TAG}_prep" -n 1 -R "rusage[mem=$PREP_MEM]" -W "$PREP_W" \
       -o "$WORK/logs/${TAG}_prep.%J.log" \
       "$LOAD && Rscript $WORK/01_prepare_features.R --beta=$PROCESSED/beta.tsv \
          --samples=$PROCESSED/master_samples.tsv --provenance=$PROCESSED/beta.provenance.json \
          --out=$FEAT $PREP_EXTRA"
  DEP=(-w "done(${TAG}_prep)")
else
  echo "Features exist in $FEAT - skipping prep."; DEP=()
fi

bsub -q "$QUEUE" -J "${TAG}_train" -n "$CORES" -R "span[hosts=1]" -R "rusage[mem=$FIT_MEM_SLOT]" -W "$FIT_W" \
     ${DEP[@]+"${DEP[@]}"} -o "$WORK/logs/${TAG}_train.%J.log" \
     "$LOAD && Rscript $WORK/02_train_loco.R --features=$FEAT --out=$OUT --cores=$CORES $FIT_EXTRA"

echo "Submitted $TAG. Monitor: bjobs -J '${TAG}*'   Logs: $WORK/logs/"
echo "Results will be in $OUT (macro_metrics.txt, loco_metrics.tsv, frozen_nonCNS.rds)"
