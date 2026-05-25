#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Run full BG6 binary analysis + movie pipeline.

Usage:
  analyze/binary/run_all_bg6.sh [options]

Options:
  --repo-root PATH         Repository root (default: inferred from script path)
  --binary-root PATH       Binary run input root (default: <repo>/binary_test)
  --output-root PATH       Analysis/movie output root (default: <repo>/analyze/binary/results)
  --levels "L1 L2 ..."      Levels list (default: "bg6_l6 bg6_l7 bg6_l8")
  --skip-analyses          Skip all analysis scripts
  --skip-movies            Skip GIF/MP4 movie generation
  --skip-six-panel         Skip combined 6-panel movie generation
  --skip-mp4               Generate GIF/frames only (skip MP4 conversion)
  --help                   Show this help

Environment overrides:
  REPO_ROOT, BINARY_ROOT, OUTPUT_ROOT, LEVELS
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/../.." && pwd)"
REPO_ROOT="${REPO_ROOT:-$REPO_ROOT_DEFAULT}"
BINARY_ROOT="${BINARY_ROOT:-$REPO_ROOT/binary_test}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$REPO_ROOT/analyze/binary/results}"
LEVELS_STR="${LEVELS:-bg6_l6 bg6_l7 bg6_l8}"
SKIP_ANALYSES=0
SKIP_MOVIES=0
SKIP_SIX_PANEL=0
SKIP_MP4=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-root)
      REPO_ROOT="$2"; shift 2 ;;
    --binary-root)
      BINARY_ROOT="$2"; shift 2 ;;
    --output-root)
      OUTPUT_ROOT="$2"; shift 2 ;;
    --levels)
      LEVELS_STR="$2"; shift 2 ;;
    --skip-analyses)
      SKIP_ANALYSES=1; shift ;;
    --skip-movies)
      SKIP_MOVIES=1; shift ;;
    --skip-six-panel)
      SKIP_SIX_PANEL=1; shift ;;
    --skip-mp4)
      SKIP_MP4=1; shift ;;
    --help|-h)
      usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2 ;;
  esac
done

IFS=' ' read -r -a LEVELS <<< "$LEVELS_STR"

LOG_DIR="$OUTPUT_ROOT/logs"
RUNS_ROOT="$OUTPUT_ROOT/bg6_runs"
AGG_DIR="$OUTPUT_ROOT/analysis_fixed_dt_bg6"
MOVIES_DIR="$OUTPUT_ROOT/movies_bg6"

mkdir -p "$LOG_DIR" "$RUNS_ROOT" "$AGG_DIR"

echo "[run_all_bg6] repo_root=$REPO_ROOT"
echo "[run_all_bg6] binary_input_root=$BINARY_ROOT"
echo "[run_all_bg6] output_root=$OUTPUT_ROOT"
echo "[run_all_bg6] levels=${LEVELS[*]}"

if [[ $SKIP_ANALYSES -eq 0 ]]; then
  echo "[run_all_bg6] Running per-level analyses"
  mkdir -p "$REPO_ROOT/amr_test/runs"

  for level in "${LEVELS[@]}"; do
    tag="equal_mass_binary_amr_${level}"
    mg_dir="$BINARY_ROOT/mg/$tag"
    fmm_dir="$BINARY_ROOT/fmm/$tag"
    runroot="$RUNS_ROOT/$level"

    if [[ ! -d "$mg_dir" || ! -d "$fmm_dir" ]]; then
      echo "Missing run directories for $level:" >&2
      echo "  $mg_dir" >&2
      echo "  $fmm_dir" >&2
      exit 1
    fi

    mkdir -p "$runroot"
    ln -sfn "$mg_dir" "$runroot/mg"
    ln -sfn "$fmm_dir" "$runroot/fmm"

    ln -sfn "$runroot" "$REPO_ROOT/amr_test/runs/latest"

    python3 "$REPO_ROOT/analyze/binary/binary_mg_fmm_compare.py" \
      > "$LOG_DIR/analysis_mg_fmm_${level}.log" 2>&1

    python3 "$REPO_ROOT/analyze/binary/binary_vs_analytic.py" "$runroot" --output-dir "$runroot/analysis" \
      > "$LOG_DIR/analysis_vs_analytic_${level}.log" 2>&1

    echo "  done analyses for $level"
  done

  echo "[run_all_bg6] Running cross-level analyses"
  depth_args=()
  lmax_args=()
  for level in "${LEVELS[@]}"; do
    lev_num="${level##*_l}"
    runroot="$RUNS_ROOT/$level"
    depth_args+=(--run "$lev_num" "$runroot")
    lmax_args+=(--run "$lev_num" "$runroot")
  done

  python3 "$REPO_ROOT/analyze/binary/binary_amr_depth_compare.py" \
    "${depth_args[@]}" \
    --output-dir "$AGG_DIR" \
    > "$LOG_DIR/analysis_amr_depth_bg6.log" 2>&1

  python3 "$REPO_ROOT/analyze/binary/binary_lmax_compilation.py" \
    "${lmax_args[@]}" \
    --output "$AGG_DIR/binary_lmax_compilation_bg6.png" \
    > "$LOG_DIR/analysis_lmax_compilation_bg6.log" 2>&1

  echo "[run_all_bg6] Building fixed_dt-style summary tables"
  python3 - "$RUNS_ROOT" "$AGG_DIR" "${LEVELS[@]}" <<'PY'
import sys
from pathlib import Path
import pandas as pd

runs_root = Path(sys.argv[1])
agg_dir = Path(sys.argv[2])
levels = sys.argv[3:]

rows = []
delta_rows = []
for level in levels:
    runroot = runs_root / level
    summary = pd.read_csv(runroot / "analysis" / "binary_solver_summary.csv")
    deltas = pd.read_csv(runroot / "analysis" / "binary_solver_deltas.csv")
    va = pd.read_csv(runroot / "analysis" / "binary_vs_analytic_summary.csv")

    for solver in ("MG", "FMM"):
        s = summary.loc[summary["solver"] == solver].iloc[0]
        v = va.loc[va["solver"] == solver].iloc[0]
        rows.append(
            {
                "solver": solver,
                "level": level.upper(),
                "t_end": float(s["time_max"]),
                "max_abs_energy_rel_drift": float(s["max_abs_energy_rel_drift"]),
                "max_abs_lmag_rel_drift": float(s["max_abs_lmag_rel_drift"]),
                "max_com_drift": float(s["max_com_drift"]),
                "max_abs_sep_rel_dev": float(s["max_abs_sep_rel_dev"]),
                "max_abs_phase_err": float(v["max_phase_err"]),
                "final_phase_err": float(v["final_phase_err"]),
                "max_body_err": float(v["max_body_err"]),
            }
        )

    t_end = float(deltas["time"].iloc[-1])
    delta_rows.append(
        {
            "level": level.upper(),
            "t_end": t_end,
            "max_abs_phase_diff": float(deltas["phase_diff"].abs().max()),
            "final_phase_diff": float(abs(deltas["phase_diff"].iloc[-1])),
            "max_abs_sep_diff": float(deltas["sep_diff"].abs().max()),
            "max_abs_com_r_diff": float(deltas["com_r_diff"].abs().max()),
            "max_abs_energy_drift_diff": float(deltas["energy_drift_diff"].abs().max()),
            "max_abs_lmag_drift_diff": float(deltas["lmag_drift_diff"].abs().max()),
        }
    )

summary_df = pd.DataFrame(rows)
deltas_df = pd.DataFrame(delta_rows)
summary_df.to_csv(agg_dir / "binary_fixed_dt_summary.csv", index=False)
deltas_df.to_csv(agg_dir / "binary_fixed_dt_solver_deltas_summary.csv", index=False)
print(agg_dir / "binary_fixed_dt_summary.csv")
print(agg_dir / "binary_fixed_dt_solver_deltas_summary.csv")
PY
fi

if [[ $SKIP_MOVIES -eq 0 ]]; then
  if [[ $SKIP_MP4 -eq 0 ]]; then
    command -v ffmpeg >/dev/null 2>&1 || {
      echo "ffmpeg not found. Re-run with --skip-mp4 or install ffmpeg." >&2
      exit 1
    }
  fi

  echo "[run_all_bg6] Building GIF movies"
  python3 "$REPO_ROOT/analyze/binary/binary_bg6_grid_movie.py" \
    --levels "${LEVELS[@]}" \
    --binary-root "$BINARY_ROOT" \
    --global-limits \
    --output-dir "$MOVIES_DIR" \
    > "$LOG_DIR/movie_bg6.log" 2>&1

  if [[ $SKIP_SIX_PANEL -eq 0 ]]; then
    echo "[run_all_bg6] Building 6-panel GIF movie"
    SIX_PANEL_DIR="$MOVIES_DIR/six_panel"
    level_tag="$(IFS=_; echo "${LEVELS[*]}")"
    level_tag="${level_tag// /_}"

    python3 "$REPO_ROOT/analyze/binary/binary_bg6_six_panel_movie.py" \
      --levels "${LEVELS[@]}" \
      --binary-root "$BINARY_ROOT" \
      --output-dir "$SIX_PANEL_DIR" \
      > "$LOG_DIR/movie_bg6_six_panel.log" 2>&1

    if [[ $SKIP_MP4 -eq 0 ]]; then
      ffmpeg -y -framerate 12 -i "$SIX_PANEL_DIR/frames/frame_%04d.png" \
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
        -c:v libx264 -pix_fmt yuv420p "$SIX_PANEL_DIR/binary_${level_tag}_mg_vs_fmm_grid_6panel.mp4" \
        > "$LOG_DIR/movie_bg6_six_panel_ffmpeg.log" 2>&1
    fi
  fi

  if [[ $SKIP_MP4 -eq 0 ]]; then
    echo "[run_all_bg6] Building MP4 movies"
    for level in "${LEVELS[@]}"; do
      ffmpeg -y -framerate 12 -i "$MOVIES_DIR/$level/frames/frame_%04d.png" \
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
        -c:v libx264 -pix_fmt yuv420p "$MOVIES_DIR/$level/binary_${level}_mg_vs_fmm_grid.mp4" \
        > "$LOG_DIR/movie_${level}_ffmpeg.log" 2>&1

      cp -f "$MOVIES_DIR/$level/binary_${level}_mg_vs_fmm_grid.gif" "$MOVIES_DIR/binary_${level}_mg_vs_fmm_grid.gif"
      cp -f "$MOVIES_DIR/$level/binary_${level}_mg_vs_fmm_grid.mp4" "$MOVIES_DIR/binary_${level}_mg_vs_fmm_grid.mp4"
    done
  fi
fi

echo "[run_all_bg6] Done"
echo "  analyses: $AGG_DIR"
echo "  movies:   $MOVIES_DIR"
echo "  logs:     $LOG_DIR"
