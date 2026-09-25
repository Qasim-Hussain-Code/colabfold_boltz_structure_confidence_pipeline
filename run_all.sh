#!/usr/bin/env bash
# =============================================================================
#  run_all.sh - the whole pipeline, in order
# =============================================================================
#  Every stage is idempotent: a finished stage prints that it is skipping and
#  returns. A run that died at stage 5 is resumed by running this again, and
#  --from restarts from a named stage after a code change.
#
#  Order, and what each stage costs on the machine this was developed on, which
#  is sixteen threads with seven gigabytes visible to the kernel and no usable
#  accelerator. The measured figures are in logs/*.resources.tsv and the README
#  quotes them from there rather than from this comment.
#
#      00_configure        seconds. Run twice: once before anything is
#                          installed, and again with --calibrate afterwards,
#                          which measures this machine rather than assuming it.
#      01_install          four environments. Measured at 97 MB, 2.5 GB, 514 MB
#                          and 1.8 GB of host disk, and about nine minutes in
#                          total when nothing is cached.
#      02_build_holdout    a few hundred requests to a public archive, cached,
#                          so a second run makes none.
#      03_fetch_references one structure per target, and one ligand per target
#                          in the docking subset.
#      04_run_msa          serial by design, not by necessity. See the note in
#                          that script.
#      05_predict          the long one. On this processor a 140-residue target
#                          took 167 seconds and a 179-residue target 258, for a
#                          single model with no alignment. Multiply by the arms,
#                          the targets and the model count, and read the
#                          measured distribution out of results/timing.tsv
#                          afterwards.
#      06_score_structures seconds per target.
#      07_validate_geometry seconds per target.
#      08_predict_pedv     the application arm, run on its own because its
#                          construct is longer than anything in the held-out
#                          set.
#      09_dock             hands the receptors to the previous stage of this
#                          roadmap, unmodified, and collects its tables.
#      10_analyse          seconds.
#      11_figures          seconds.
#      12_report           one render.
#
#  Usage:
#      bash run_all.sh                          # everything, resuming
#      bash run_all.sh --from 05_predict        # restart at prediction
#      bash run_all.sh --arm af2_msa_notmpl     # one arm only
#      bash run_all.sh --smoke                  # a few targets, one arm
#      bash run_all.sh --help
#
#  Options:
#      --from STAGE    start at this stage name
#      --only STAGE    run just this stage
#      --arm NAME      restrict prediction and scoring to one arm
#      --limit N       cap the targets per stage
#      --smoke         three targets, one arm: proves the wiring
#      --force         pass --force to every stage it applies to
#      --skip-pedv     leave out the application arm
#      --skip-docking  leave out the handoff to the docking pipeline
#      -h, --help      this text
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${ROOT}/scripts"

STAGES=(00_configure 01_install 02_build_holdout 03_fetch_references 04_run_msa
        05_predict 06_score_structures 07_validate_geometry 08_predict_pedv
        09_dock 10_analyse 11_figures 12_report)

FROM=""; ONLY=""; ARM=""; LIMIT=""; SMOKE=0; FORCE=0; SKIP_PEDV=0; SKIP_DOCK=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)         FROM="$2"; shift 2 ;;
        --only)         ONLY="$2"; shift 2 ;;
        --arm)          ARM="$2"; shift 2 ;;
        --limit)        LIMIT="$2"; shift 2 ;;
        --smoke)        SMOKE=1; shift ;;
        --force)        FORCE=1; shift ;;
        --skip-pedv)    SKIP_PEDV=1; shift ;;
        --skip-docking) SKIP_DOCK=1; shift ;;
        -h|--help)      sed -n '2,62p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

stage_index() {
    local want="$1" i
    for i in "${!STAGES[@]}"; do
        [[ "${STAGES[$i]}" == "$want" ]] && { echo "$i"; return; }
    done
    echo "-1"
}
START=0
if [[ -n "$FROM" ]]; then
    START="$(stage_index "$FROM")"
    (( START >= 0 )) || { echo "[error] unknown stage: ${FROM}" >&2; exit 1; }
fi
should_run() {
    local name="$1"
    if [[ -n "$ONLY" ]]; then [[ "$ONLY" == "$name" ]] && return 0 || return 1; fi
    local i; i="$(stage_index "$name")"
    (( i >= START ))
}
banner() {
    echo
    echo "=============================================================="
    echo " $1"
    echo "=============================================================="
}

FORCE_FLAG=(); (( FORCE == 1 )) && FORCE_FLAG=(--force)
if (( SMOKE == 1 )); then
    [[ -n "$LIMIT" ]] || LIMIT=3
    [[ -n "$ARM" ]] || ARM="af2_nomsa"
    echo "[run_all] smoke run: ${LIMIT} targets, arm ${ARM}"
fi
LIMIT_FLAG=(); [[ -n "$LIMIT" ]] && LIMIT_FLAG=(--limit "$LIMIT")

# -----------------------------------------------------------------------------
# 00. Only when project.conf is absent. Re-running it would overwrite a
# deliberately edited budget on a resumed run, and would throw away a measured
# memory curve in favour of a projection.
# -----------------------------------------------------------------------------
if should_run 00_configure; then
    if [[ ! -f "${ROOT}/project.conf" || "$ONLY" == "00_configure" ]]; then
        banner "00_configure"
        bash "${SCRIPTS}/00_configure.sh" --yes
    else
        echo "[00_configure] project.conf exists; skipping (--only 00_configure to redo)."
    fi
fi
[[ -f "${ROOT}/project.conf" ]] || {
    echo "[error] project.conf missing. Run: bash scripts/00_configure.sh --yes" >&2; exit 1; }
# shellcheck source=/dev/null
source "${ROOT}/project.conf"

PY="$(dirname "$(dirname "$CONDA_SH")")"
PY="$(dirname "$PY")/envs/${CONDA_ENV_ANALYSIS}/bin/python"

if should_run 01_install; then
    banner "01_install"
    bash "${SCRIPTS}/01_install.sh" "${FORCE_FLAG[@]}"
    # The memory ceiling in project.conf is a projection until this runs. It
    # predicts three short targets and refits the curve from what this machine
    # actually did.
    if [[ ! -s "${RESULTS_DIR}/environment/memory_curve.tsv" ]]; then
        banner "00_configure --calibrate"
        bash "${SCRIPTS}/00_configure.sh" --calibrate || \
            echo "[run_all] calibration failed; the ceiling stays a projection and says so"
    fi
fi

if should_run 02_build_holdout; then
    banner "02_build_holdout_set"
    "$PY" "${SCRIPTS}/02_build_holdout_set.py" --config "${ROOT}/project.conf" \
        "${FORCE_FLAG[@]}"
fi
if should_run 03_fetch_references; then
    banner "03_fetch_references"
    bash "${SCRIPTS}/03_fetch_references.sh" "${FORCE_FLAG[@]}" "${LIMIT_FLAG[@]}"
fi
if should_run 04_run_msa; then
    banner "04_run_msa"
    bash "${SCRIPTS}/04_run_msa.sh" "${FORCE_FLAG[@]}" "${LIMIT_FLAG[@]}"
fi

# -----------------------------------------------------------------------------
# 05. One invocation per arm. The arms that need an alignment are skipped with
# a message when none is present, rather than failing, so a run without the
# alignment stage still produces the single-sequence arm.
# -----------------------------------------------------------------------------
if should_run 05_predict; then
    banner "05_predict"
    IFS=',' read -r -a ARM_LIST <<< "${ARM:-$ARMS}"
    for arm in "${ARM_LIST[@]}"; do
        [[ -n "$arm" ]] || continue
        case "$arm" in
            null_template|null_unrelated) continue ;;   # built by the scoring stage
        esac
        bash "${SCRIPTS}/05_predict.sh" --arm "$arm" "${FORCE_FLAG[@]}" "${LIMIT_FLAG[@]}" || \
            echo "[run_all] arm ${arm} failed; continuing with the rest"
    done
    # The five-model protocol on a subset, so that what model selection is
    # worth is measured rather than assumed.
    if [[ -z "$ARM" && "${AF2_SUBSET_N_TARGETS:-0}" != "0" ]]; then
        banner "05_predict, five models on a subset"
        bash "${SCRIPTS}/05_predict.sh" --arm af2_msa_notmpl \
            --num-models "${AF2_NUM_MODELS_SUBSET}" \
            --limit "${AF2_SUBSET_N_TARGETS}" --force || \
            echo "[run_all] the five-model subset failed; continuing"
    fi
    # Repeat seeds on one target, because run-to-run spread is part of the
    # measurement and is almost never reported.
    banner "05_predict --seed-variance"
    bash "${SCRIPTS}/05_predict.sh" --arm "${ARM:-af2_msa_notmpl}" --seed-variance \
        --limit "${SEED_VARIANCE_N_TARGETS:-1}" || \
        echo "[run_all] the seed experiment failed; continuing"
fi

if should_run 06_score_structures; then
    # The two floors first. They are copied structures rather than
    # predictions, so they cost a download each and no inference, and the
    # scoring stage treats them exactly as it treats a prediction.
    banner "06a_build_null_floors"
    "$PY" "${SCRIPTS}/06a_build_null_floors.py" --config "${ROOT}/project.conf"         "${FORCE_FLAG[@]}" "${LIMIT_FLAG[@]}" ||         echo "[run_all] the floors failed; continuing without them"
    banner "06_score_structures"
    "$PY" "${SCRIPTS}/06_score_structures.py" --config "${ROOT}/project.conf" \
        ${ARM:+--arm "$ARM"} "${FORCE_FLAG[@]}"
fi
if should_run 07_validate_geometry; then
    banner "07_validate_geometry"
    "$PY" "${SCRIPTS}/07_validate_geometry.py" --config "${ROOT}/project.conf" \
        ${ARM:+--arm "$ARM"} "${FORCE_FLAG[@]}"
fi
if should_run 08_predict_pedv && (( SKIP_PEDV == 0 )) && (( SMOKE == 0 )); then
    banner "08_predict_pedv"
    bash "${SCRIPTS}/08_predict_pedv.sh" "${FORCE_FLAG[@]}" || \
        echo "[run_all] the application arm failed; continuing"
fi
if should_run 09_dock && (( SKIP_DOCK == 0 )) && (( SMOKE == 0 )); then
    banner "09_dock_into_predictions"
    for arm in crystal predicted; do
        bash "${SCRIPTS}/09_dock_into_predictions.sh" --arm "$arm" "${FORCE_FLAG[@]}" || \
            echo "[run_all] the ${arm} docking arm failed; continuing"
    done
fi
if should_run 10_analyse; then
    banner "10_analyse"
    "$PY" "${SCRIPTS}/10_analyse.py" --config "${ROOT}/project.conf"
fi
if should_run 11_figures; then
    banner "11_figures"
    "$PY" "${SCRIPTS}/11_figures.py" --config "${ROOT}/project.conf"
fi
if should_run 12_report; then
    banner "12_report"
    QUARTO="$(dirname "$PY")/quarto"
    if [[ -x "$QUARTO" ]] || command -v quarto >/dev/null 2>&1; then
        [[ -x "$QUARTO" ]] || QUARTO="$(command -v quarto)"
        mkdir -p "${RESULTS_DIR}/report"
        # Rendered from the repository root, with the config named explicitly.
        # The report locates project.conf by walking up from the working
        # directory, and a report that cannot find its own config depending on
        # how it was invoked is not reproducible.
        ( cd "$ROOT" && CSC_CONFIG="${ROOT}/project.conf" \
            "$QUARTO" render "scripts/12_report.qmd" --to html \
                --output-dir "${RESULTS_DIR}/report" ) || \
            echo "[12_report] the render failed; the tables in results/ are unaffected"
    else
        echo "[12_report] quarto not on PATH; skipping."
        echo "            Every number it would show is already in results/."
    fi
fi

banner "done"
echo "Tables:  ${RESULTS_DIR}"
echo "Figures: ${FIGURES_DIR}"
echo "Logs:    ${LOG_DIR}"
if ls "${LOG_DIR}"/*.resources.tsv >/dev/null 2>&1; then
    echo
    echo "Per-stage elapsed time, peak resident memory and disk:"
    { printf 'stage\telapsed_s\tpeak_rss_mb\tdata_growth_mb\tdata_total_mb\tcache_mb\tpeak_disk_mb\thost_start_mb\thost_min_mb\tjobs\tnote\n'
      for f in "${LOG_DIR}"/*.resources.tsv; do
          [[ -f "$f" ]] && tail -n +2 "$f"
      done
    } > "${LOG_DIR}/summary.tsv"
    column -t -s "$(printf '\t')" "${LOG_DIR}/summary.tsv" | cut -c1-140 | sed 's/^/  /'
fi
