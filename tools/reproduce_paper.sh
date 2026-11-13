#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Isolate the checkout's source version from older user-library installations.
paper_library=$(mktemp -d /tmp/capnet-paper-library.XXXXXX)
export R_LIBS_USER="$paper_library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PAPER_WORKERS="${PAPER_WORKERS:-4}"
R CMD INSTALL --no-docs --no-html --no-multiarch --library="$paper_library" .
Rscript tools/run_paper_experiments.R all
Rscript tools/check_paper_results.R
Rscript tools/render_paper_vignettes.R
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd manuscript/main.tex
