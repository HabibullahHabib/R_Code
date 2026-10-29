#!/bin/bash
set -euo pipefail
export PATH=/workspace/.r-runtime/bin:$PATH
export R_LIBS_SITE=/workspace/.r-runtime/usr/lib/R/site-library
export R_LIBS_USER=/workspace/.r-runtime/library
export LD_LIBRARY_PATH=/workspace/.r-runtime/usr/lib/x86_64-linux-gnu
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export XDG_CACHE_HOME=/workspace/.r-runtime/cache
mkdir -p "$XDG_CACHE_HOME"
export SIM_WORKERS=${SIM_WORKERS:-4}
mode=${1:-full}
exec flock -n "/workspace/paper-feedback-10-treatments/${mode}.lock" Rscript --vanilla /workspace/paper-feedback-10-treatments/run_all.R "/workspace/paper-feedback-10-treatments" "$mode"
