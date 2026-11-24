#!/bin/sh
# From any working directory; pass --hash to reauthenticate in a clean archive.
set -eu
here=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
taskset -c "${CPU:-5}" nice -n 19 "${PYTHON:-python3}" -Bu "$here/../replay.py" --caps 18 "$@"
taskset -c "${CPU:-5}" nice -n 19 "${PYTHON:-python3}" -Bu "$here/replay.py" --host 175 --tip 5 --stock 8 "$@"
