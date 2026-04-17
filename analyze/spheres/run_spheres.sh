#!/bin/bash
#SBATCH --job-name=spheres-compare
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=7500M
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user=jl4415@princeton.edu
#SBATCH --mail-type=END,FAIL

module purge
module load openmpi/gcc/4.1.2

cd /home/jl4415/mini-ramses || exit 1

export SLURM_CPU_BIND=cores
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
SCRIPT_DIR="${REPO_ROOT}/analyze/spheres"

PYTHON_BIN="${PYTHON_BIN:-python3}"
SOLVER="${SOLVER:-both}"            # both | fmm | mg
FMM_EXEC="${FMM_EXEC:-${REPO_ROOT}/fmm_spheres}"
MG_EXEC="${MG_EXEC:-${REPO_ROOT}/mg_spheres}"
LEVEL_START="${LEVEL_START:-6}"
LEVEL_END="${LEVEL_END:-10}"
NPROCS="${NPROCS:-${SLURM_NTASKS:-1}}"
RUN_ROOT="${RUN_ROOT:-${REPO_ROOT}/runs/spheres_compare}"
NML_TEMPLATE="${NML_TEMPLATE:-${REPO_ROOT}/namelist/spheres.nml}"
NML_DIR="${NML_DIR:-}"
ANALYZE_NOUT="${ANALYZE_NOUT:-2}"
DO_ANALYZE="${DO_ANALYZE:-1}"
CONV_CSV="${CONV_CSV:-}"
CONV_PLOT="${CONV_PLOT:-}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --solver MODE        Solver mode: both|fmm|mg (default: both)
  --fmm-exec PATH      FMM executable (default: ${REPO_ROOT}/fmm_spheres)
  --mg-exec PATH       MG executable (default: ${REPO_ROOT}/mg_spheres)
  --level-start N      First level (default: 6)
  --level-end N        Last level (default: 10)
  --nprocs N           MPI ranks (default: SLURM_NTASKS or 1)
  --run-root PATH      Root output directory (default: runs/spheres_compare)
  --nout N             Output index for analysis (default: 2)
  --no-analyze         Skip analysis/plot stage
  --help               Show this message

Environment variables are also supported:
  SOLVER, FMM_EXEC, MG_EXEC, LEVEL_START, LEVEL_END, NPROCS,
  RUN_ROOT, NML_TEMPLATE, NML_DIR, ANALYZE_NOUT, DO_ANALYZE,
  CONV_CSV, CONV_PLOT, PYTHON_BIN
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --solver)
      SOLVER="$2"
      shift 2
      ;;
    --fmm-exec)
      FMM_EXEC="$2"
      shift 2
      ;;
    --mg-exec)
      MG_EXEC="$2"
      shift 2
      ;;
    --level-start)
      LEVEL_START="$2"
      shift 2
      ;;
    --level-end)
      LEVEL_END="$2"
      shift 2
      ;;
    --nprocs)
      NPROCS="$2"
      shift 2
      ;;
    --run-root)
      RUN_ROOT="$2"
      shift 2
      ;;
    --nout)
      ANALYZE_NOUT="$2"
      shift 2
      ;;
    --no-analyze)
      DO_ANALYZE=0
      shift 1
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1
      ;;
  esac
done

case "$SOLVER" in
  both|fmm|mg) ;;
  *)
    echo "Invalid --solver '$SOLVER' (must be both|fmm|mg)" >&2
    exit 1
    ;;
esac

if (( LEVEL_END < LEVEL_START )); then
  echo "Invalid level range: LEVEL_END < LEVEL_START" >&2
  exit 1
fi

if [[ -z "$NML_DIR" ]]; then
  NML_DIR="${RUN_ROOT}/namelists"
fi
if [[ -z "$CONV_CSV" ]]; then
  CONV_CSV="${RUN_ROOT}/spheres_convergence.csv"
fi
if [[ -z "$CONV_PLOT" ]]; then
  CONV_PLOT="${RUN_ROOT}/spheres_convergence.png"
fi

mkdir -p "$RUN_ROOT" "$NML_DIR"

"$PYTHON_BIN" "${SCRIPT_DIR}/generate_spheres_namelists.py" \
  --template "$NML_TEMPLATE" \
  --output-dir "$NML_DIR" \
  --level-start "$LEVEL_START" \
  --level-end "$LEVEL_END"

if command -v srun >/dev/null 2>&1; then
  LAUNCHER=(srun -n "$NPROCS")
elif command -v mpirun >/dev/null 2>&1 && (( NPROCS > 1 )); then
  LAUNCHER=(mpirun -np "$NPROCS")
else
  LAUNCHER=()
fi

run_case() {
  local case_name="$1"
  local executable="$2"

  if [[ ! -x "$executable" ]]; then
    echo "Executable not found or not executable for case '$case_name': $executable" >&2
    exit 1
  fi

  for level in $(seq "$LEVEL_START" "$LEVEL_END"); do
    local level_tag
    level_tag=$(printf "%02d" "$level")

    local src_nml="${NML_DIR}/spheres_l${level_tag}.nml"
    local run_dir="${RUN_ROOT}/${case_name}/level_${level_tag}"
    local runtime_nml_name="spheres_runtime_l${level_tag}.nml"
    local runtime_nml="${run_dir}/${runtime_nml_name}"

    if [[ ! -f "$src_nml" ]]; then
      echo "Missing generated namelist: $src_nml" >&2
      exit 1
    fi

    mkdir -p "$run_dir"
    rm -rf "${run_dir}/output_"* "${run_dir}/time_"*.txt

    # If namelist has initfile='namelist/...', rewrite to absolute repo path.
    sed -E "s#(initfile\([0-9]+\)\s*=\s*')namelist/#\1${REPO_ROOT}/namelist/#g" \
      "$src_nml" > "$runtime_nml"

    case_label=$(echo "$case_name" | tr '[:lower:]' '[:upper:]')
    echo "=== Case ${case_label} | level ${level} | run dir: ${run_dir} ==="
    (
      cd "$run_dir"
      if [[ ${#LAUNCHER[@]} -gt 0 ]]; then
        "${LAUNCHER[@]}" "$executable" "$runtime_nml_name"
      else
        "$executable" "$runtime_nml_name"
      fi
    )
  done
}

#if [[ "$SOLVER" == "both" || "$SOLVER" == "fmm" ]]; then
#  run_case "fmm" "$FMM_EXEC"
#fi

#if [[ "$SOLVER" == "both" || "$SOLVER" == "mg" ]]; then
#  run_case "mg" "$MG_EXEC"
#fi

if [[ "$DO_ANALYZE" == "1" ]]; then
  ANALYZE_CMD=(
    "$PYTHON_BIN" "${SCRIPT_DIR}/spheres_convergence.py"
    --level-min "$LEVEL_START"
    --level-max "$LEVEL_END"
    --nout "$ANALYZE_NOUT"
    --csv "$CONV_CSV"
    --plot "$CONV_PLOT"
  )

  if [[ "$SOLVER" == "both" ]]; then
    ANALYZE_CMD+=(--case "FMM=${RUN_ROOT}/fmm" --case "MG=${RUN_ROOT}/mg")
  elif [[ "$SOLVER" == "fmm" ]]; then
    ANALYZE_CMD+=(--case "FMM=${RUN_ROOT}/fmm")
  else
    ANALYZE_CMD+=(--case "MG=${RUN_ROOT}/mg")
  fi

  "${ANALYZE_CMD[@]}"
fi

echo "Finished spheres run: solver=${SOLVER}, levels=${LEVEL_START}..${LEVEL_END}"
