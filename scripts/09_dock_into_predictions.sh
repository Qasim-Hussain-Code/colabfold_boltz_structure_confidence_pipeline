#!/usr/bin/env bash
# =============================================================================
#  09_dock_into_predictions.sh - hand the predicted receptors to stage 1
# =============================================================================
#  This stage runs no docking code of its own. It clones the previous stage of
#  this roadmap at a pinned commit, prepares a manifest in the format that
#  pipeline reads, calls its scripts unmodified, and collects its tables. The
#  point of the comparison is that the docking is identical and only the
#  receptor differs, which is only true if the same code runs.
#
#  Two arms, and the second is what makes the first mean anything.
#
#      predicted   dock into the predicted receptor
#      crystal     dock the same ligand into the experimental receptor of the
#                  same target, with the same protocol
#
#  Without the crystal arm the predicted number is uninterpretable: this set is
#  not the set the previous stage reported on, so its published figure is not a
#  fair comparator. The crystal arm is that comparator, measured here, on these
#  targets. The published figure is quoted alongside as context and the README
#  says plainly which is which.
#
#  Three fixes were needed in that repository before this could work, and all
#  three were made there rather than worked around here, as the brief for this
#  stage requires. The first: its scoring stage looked for reference ligands
#  only under the two dataset names it shipped with, so a run on any other
#  dataset was skipped silently, with no message and no row. The second: when
#  the protonation step falls back to its second tool, the receptor file it
#  writes carries header records that the docking program rejects, which turned
#  a protonation fallback into a docking failure for every affected receptor.
#  The third: ligand preparation ran the whole set through one worker pool with
#  no per-item handling, so a ligand the charge model could not parameterise
#  ended the preparation of every other ligand with a traceback and no table.
#  One heme in this set does exactly that. The pinned commit carries all three.
#
#  A coordinate frame problem this stage has to solve, and it is the one most
#  likely to produce a wrong answer quietly. That pipeline measures pose
#  accuracy in place, with no superposition, and centres its search box on the
#  reference ligand. A predicted structure comes out in its own arbitrary
#  frame. Docking into it with a box placed on the crystal ligand would put the
#  box somewhere in solvent. So each predicted receptor is superposed onto its
#  experimental structure first, and the superposition is recorded per target.
#
#  Usage:
#      bash scripts/09_dock_into_predictions.sh --arm predicted
#      bash scripts/09_dock_into_predictions.sh --arm crystal
#
#  Options:
#      --arm NAME   predicted or crystal
#      --limit N    stop after N targets
#      --jobs N     dock N targets at once, default DOCK_JOBS in
#                   project.conf or 1
#      --force      ignore the stage stamp and redo
#      -h, --help   this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
csc_load_conf

ARM="predicted"; LIMIT=0; FORCE=0; DOCK_JOBS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arm)     ARM="$2"; shift 2 ;;
        --limit)   LIMIT="$2"; shift 2 ;;
        --jobs)    DOCK_JOBS="$2"; shift 2 ;;
        --force)   FORCE=1; shift ;;
        -h|--help) sed -n '2,48p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done
case "$ARM" in
    predicted|crystal) ;;
    *) echo "[error] --arm must be predicted or crystal" >&2; exit 1 ;;
esac

# How many targets are docked at once. This is not a setting of the docking
# protocol and changing it cannot change a result: each search runs with one
# CPU and a fixed seed, so it is deterministic whatever else is running. It
# matters because the ligands here are nucleotides and cofactors rather than
# the drug-like ligands the docking pipeline was benchmarked on, and one of
# them took ten minutes where that pipeline averaged eighteen seconds. Run
# one at a time on this set, a single arm is most of a day.
[[ -n "$DOCK_JOBS" ]] || DOCK_JOBS="${DOCKING_JOBS:-1}"
[[ "$DOCK_JOBS" =~ ^[0-9]+$ ]] || DOCK_JOBS=1

STAGE="09_dock_${ARM}"
csc_skip_if_done "$STAGE" "$FORCE" && exit 0

PY_ANALYSIS="$(csc_env_bin "$CONDA_ENV_ANALYSIS" python)"
SUBSET="${CONFIG_DIR}/docking_subset.tsv"
[[ -s "$SUBSET" ]] || { echo "[error] ${SUBSET} missing; run 02_build_holdout_set.py" >&2; exit 3; }

csc_stage_start "$STAGE"

# -----------------------------------------------------------------------------
# 1. The previous stage, at the pinned commit.
# -----------------------------------------------------------------------------
# shellcheck disable=SC2153
if [[ ! -d "${STAGE1_DIR}/.git" ]]; then
    echo "[${STAGE}] cloning the docking pipeline at ${STAGE1_COMMIT:0:12}"
    mkdir -p "$(dirname "$STAGE1_DIR")"
    csc_run "clone_stage1" git clone --quiet "$STAGE1_REPO" "$STAGE1_DIR"
fi
( cd "$STAGE1_DIR" && git fetch --quiet --all && git checkout --quiet "$STAGE1_COMMIT" )
GOT="$(cd "$STAGE1_DIR" && git rev-parse HEAD)"
if [[ "$GOT" != "$STAGE1_COMMIT" ]]; then
    echo "[error] the clone is at ${GOT}, not the pinned ${STAGE1_COMMIT}" >&2
    exit 3
fi
echo "[${STAGE}] docking pipeline at $(cd "$STAGE1_DIR" && git rev-parse --short HEAD)"

# That repository reads its settings from a file at its own root. It is written
# here rather than by its own configuration script, because its script sizes a
# disk budget for its own datasets and would refuse a budget appropriate to
# this one. Every path is literal: it reads the file without shell expansion.
S1_DATA="${DATA_DIR}/stage1_data"
S1_RESULTS="${DATA_DIR}/stage1_results"
mkdir -p "$S1_DATA" "$S1_RESULTS" "${S1_RESULTS}/runs"
cat > "${STAGE1_DIR}/project.conf" <<CONF
THREADS=${THREADS}
RAM_GB=${RAM_GB}
DISK_GB=${DISK_GB}
JOBS=1
ARMS=${STAGE1_ARMS}
METHODS=${STAGE1_METHODS}
DATASETS=${ARM}
SEED=${SEED}
SEED_REPLICATES="${SEED}"
BOX_SIZE=25
EXHAUSTIVENESS=8
TOP_N_POSES=5
GNINA_CNN_SCORING=rescore
GNINA_CNN_MODEL=crossdock_default2018
RMSD_PASS=2.0
RMSD_STRICT=1.0
XDOCK_MAX_RESOLUTION=2.5
XDOCK_MIN_SEQ_IDENTITY=0.95
XDOCK_MAX_PARTNERS=1
POCKET_METHOD=auto
REPO_DIR=${STAGE1_DIR}
CONFIG_DIR=${STAGE1_DIR}/config
SCRIPTS_DIR=${STAGE1_DIR}/scripts
DATA_DIR=${S1_DATA}
RESULTS_DIR=${S1_RESULTS}
FIGURES_DIR=${S1_RESULTS}/figures
LOG_DIR=${S1_RESULTS}/logs
CONDA_SH=${CONDA_SH}
CONDA_ENV_NAME=vgb_bench
CONDA_ENV_VINA=vgb_vina
CONDA_ENV_SMINA=vgb_smina
CONDA_ENV_GNINA_RT=vgb_gnina
GNINA_BIN=${S1_DATA}/tools/gnina
CONF
# Its results directory is pointed outside its own tree on purpose. That
# repository tracks its own results, and letting this stage write into them
# would mix two experiments in one set of tables.

# -----------------------------------------------------------------------------
# 2. Receptors in the docking pipeline's format, in the crystal frame.
# -----------------------------------------------------------------------------
csc_run "prepare_inputs" "$PY_ANALYSIS" "${SCRIPT_DIR}/lib_dock_handoff.py" \
    --config "${REPO_DIR}/project.conf" --arm "$ARM" \
    --stage1-config "${STAGE1_DIR}/project.conf" \
    ${LIMIT:+--limit "$LIMIT"}

MANIFEST="${STAGE1_DIR}/config/dataset_${ARM}.tsv"
[[ -s "$MANIFEST" ]] || { echo "[error] no manifest was written for arm ${ARM}" >&2; exit 4; }
N_TARGETS="$(( $(wc -l < "$MANIFEST") - 1 ))"
echo "[${STAGE}] ${N_TARGETS} receptors prepared for the ${ARM} arm"
(( N_TARGETS > 0 )) || { echo "[error] the manifest is empty" >&2; exit 4; }

# -----------------------------------------------------------------------------
# 3. That pipeline's own stages, unmodified.
# -----------------------------------------------------------------------------
set +u
# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate vgb_bench
set -u

# Its preparation stage writes its outcome tables fresh on each call, so the
# two arms would overwrite each other's record. Each call's tables are copied
# aside immediately and merged when both arms have run.
run_stage1() {
    local what="$1"; shift
    csc_run "stage1_${what}" python "${STAGE1_DIR}/scripts/${what}" "$@"
}

run_stage1 04_prepare.py --config "${STAGE1_DIR}/project.conf" \
    --dataset "$ARM" --what receptors --jobs 1
for f in receptor_prep ligand_prep; do
    [[ -f "${S1_RESULTS}/preparation/${f}.tsv" ]] && \
        cp "${S1_RESULTS}/preparation/${f}.tsv" "${S1_RESULTS}/preparation/${f}_${ARM}_receptors.tsv"
done
run_stage1 04_prepare.py --config "${STAGE1_DIR}/project.conf" \
    --dataset "$ARM" --what ligands --jobs 1
for f in receptor_prep ligand_prep; do
    [[ -f "${S1_RESULTS}/preparation/${f}.tsv" ]] && \
        cp "${S1_RESULTS}/preparation/${f}.tsv" "${S1_RESULTS}/preparation/${f}_${ARM}_ligands.tsv"
done

# Its box stage writes one file for the whole run and skips when that file
# exists, so the second arm needs --force. Pocket detection is off: the arm
# that needs it is not run here.
run_stage1 05_define_boxes.py --config "${STAGE1_DIR}/project.conf" \
    --dataset "$ARM" --no-detected --jobs 1 --force

for method in $(csc_split "$STAGE1_METHODS"); do
    csc_run "stage1_dock_${method}" bash "${STAGE1_DIR}/scripts/06_dock.sh" \
        --arm A2_genconf_refbox --method "$method" --dataset "$ARM" --jobs "$DOCK_JOBS"
done
# The floor, run once with the first method, exactly as that pipeline runs it.
csc_run "stage1_dock_null" bash "${STAGE1_DIR}/scripts/06_dock.sh" \
    --arm A0_null --method "$(csc_split "$STAGE1_METHODS" | head -1)" \
    --dataset "$ARM" --jobs "$DOCK_JOBS" || echo "[${STAGE}] the floor arm failed; continuing"

csc_run "stage1_score" python "${STAGE1_DIR}/scripts/07_score_poses.py" \
    --config "${STAGE1_DIR}/project.conf" --jobs 1 --force

# -----------------------------------------------------------------------------
# 4. Collect, and compute the rates the same way that pipeline does.
# -----------------------------------------------------------------------------
csc_run "collect_results" "$PY_ANALYSIS" "${SCRIPT_DIR}/lib_dock_handoff.py" \
    --config "${REPO_DIR}/project.conf" --arm "$ARM" \
    --stage1-config "${STAGE1_DIR}/project.conf" --collect

csc_stage_end "arm=${ARM} targets=${N_TARGETS}"
# A run that covered part of the set must not mark the stage finished. The
# prediction stage had this exact fault: a single-target test stamped the arm,
# and the full run that followed skipped it and reported success having done
# nothing.
if (( LIMIT > 0 )); then
    echo "[${STAGE}] this run covered ${N_TARGETS} targets under --limit, so the"
    echo "           stage is not marked as finished. Re-run without --limit."
else
    csc_mark_done "$STAGE"
fi
