#!/usr/bin/env bash
# =============================================================================
#  08_predict_pedv.sh - the application arm
# =============================================================================
#  A single-conformation model returns one structure and one set of confidence
#  numbers. This arm asks what that output does with a domain that the
#  experimental record shows in more than one place, and whether anything the
#  model reports tells you that it moves.
#
#  The question is not which deposited structure is correct. Several entries of
#  this protein exist and they place the N-terminal domain differently relative
#  to the body. This stage does not adjudicate between them. It measures where
#  the prediction puts that domain, against each deposited arrangement, and
#  reports which one it resembles, and it stops there.
#
#  What is measured rather than assumed
#  ------------------------------------
#  The domain boundary is taken from the deposited coordinates rather than from
#  a number quoted in a paper: 08a_pedv_analysis.py reads the entries, finds
#  the domain by its own contacts, and writes the residue range it used. The
#  displacement between deposited arrangements is measured here too, by
#  superposing on the body and reporting how far the domain moves. If the
#  deposited models turn out to agree more closely than the literature implies,
#  that is what the table will say.
#
#  The construct, and what it costs
#  --------------------------------
#  Predicting the assembled trimer is out of reach on this machine and would
#  answer a different question anyway. Predicting the domain alone cannot
#  answer this one at all: a domain on its own has no body to sit against. So
#  the construct is the domain plus enough of the body to define the interface,
#  as one chain. The length settled on is written into the manifest with the
#  reasoning, and the README says what a single-chain construct cannot tell
#  you: nothing about how three copies pack against each other.
#
#  Which arm is a prediction here
#  ------------------------------
#  What runs here is AlphaFold2 twice, once with an alignment and once without.
#  Its recorded cutoff is 30 April 2018 and every deposited entry of this
#  protein was released in 2019 or later, so with templates off both arms are
#  predicting rather than recalling, and the pair says how much of whatever
#  the model produces comes from the alignment.
#
#  The comparison that is missing, and why. Boltz-2's recorded cutoff is 1 June
#  2023, which is after all six of these entries, so those coordinates were
#  available to it during training and the same construct run through it would
#  be an uneven and more interesting comparison. It is not run here: the weights
#  and their run do not fit in the disk this machine has left, which is recorded
#  rather than presented as a choice. The README says the same, because an arm
#  that was planned and not run is a limitation and not an omission to be
#  quietly dropped.
#
#  Usage:
#      bash scripts/08_predict_pedv.sh
#      bash scripts/08_predict_pedv.sh --entries 6U7K,6VV5
#
#  Options:
#      --entries LIST  comma-separated accessions to measure against
#      --no-predict    measure the deposited entries and stop
#      --force         ignore the stage stamp and redo
#      -h, --help      this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
csc_load_conf

STAGE=08_predict_pedv
ENTRIES="${PEDV_ENTRIES:-}"
NO_PREDICT=0; FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --entries)    ENTRIES="$2"; shift 2 ;;
        --no-predict) NO_PREDICT=1; shift ;;
        --force)      FORCE=1; shift ;;
        -h|--help)    sed -n '2,54p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$ENTRIES" ]] || ENTRIES="6U7K,6VV5,7W6M,7W73,7Y6S,7Y6T"

csc_skip_if_done "$STAGE" "$FORCE" && exit 0

PY_ANALYSIS="$(csc_env_bin "$CONDA_ENV_ANALYSIS" python)"
[[ -n "$PY_ANALYSIS" ]] || { echo "[error] ${CONDA_ENV_ANALYSIS} not installed; run 01_install.sh" >&2; exit 3; }

csc_stage_start "$STAGE"

PEDV_DIR="${DATA_DIR}/pedv"
mkdir -p "$PEDV_DIR"

# --- the deposited entries, measured -----------------------------------------
# Every number this arm reports about the experimental record comes from the
# coordinates, fetched here and checksummed like any other input.
for entry in $(csc_split "$ENTRIES"); do
    dest="${PEDV_DIR}/${entry}.cif.gz"
    if [[ -s "$dest" ]]; then continue; fi
    echo "[${STAGE}] fetching ${entry}"
    csc_run "fetch_${entry}" curl -sS --fail --max-time 300 --retry 3 \
        -o "$dest" "https://files.rcsb.org/download/${entry}.cif.gz" || \
        { echo "[${STAGE}] ${entry} could not be fetched"; rm -f "$dest"; }
done

# How long a construct this machine can actually predict from an alignment,
# fitted from the peak memory every finished prediction recorded rather than
# from the calibration curve in project.conf. That curve was measured on
# single-sequence runs, which cost far less memory at the same length: 2.9 GB
# against 4.1 GB at 179 residues here, and the gap widens with length. Using
# it for an alignment run lets a construct through that the kernel then kills.
#
# With too little data to fit, this returns 0 and the cap in project.conf is
# the only one that applies.
FITS_IN_MEMORY="$("$PY_ANALYSIS" "${SCRIPT_DIR}/lib_memory.py" --config "${REPO_DIR}/project.conf" --arms af2_msa --budget-mb "$(( RAM_GB * 1024 ))" --quiet 2>/dev/null || echo 0)"
[[ "$FITS_IN_MEMORY" =~ ^[0-9]+$ ]] || FITS_IN_MEMORY=0
if (( FITS_IN_MEMORY > 0 )); then
    echo "[${STAGE}] the measured memory curve allows ${FITS_IN_MEMORY} residues"
else
    echo "[${STAGE}] not enough finished predictions to fit a memory curve;"
    echo "           only the cap in project.conf applies"
fi

csc_run "measure_deposited" "$PY_ANALYSIS" "${SCRIPT_DIR}/08a_pedv_analysis.py" \
    --config "${REPO_DIR}/project.conf" --entries "$ENTRIES" --measure     --max-construct "$FITS_IN_MEMORY"

if (( NO_PREDICT == 1 )); then
    csc_stage_end "measured only"
    exit 0
fi

# --- the construct ------------------------------------------------------------
# Written by the analysis script from the residue range it derived, so the
# sequence predicted and the range measured cannot drift apart.
CONSTRUCT="${PEDV_DIR}/construct.fasta"
[[ -s "$CONSTRUCT" ]] || { echo "[error] no construct was written; see the measurement step" >&2; exit 4; }
LEN="$(awk 'NR>1 {n += length($0)} END {print n}' "$CONSTRUCT")"
echo "[${STAGE}] construct is ${LEN} residues"

# The same ceiling every other target is held to.
PEAK="$(awk -v L="$LEN" -v b="${MEM_BASE_MB}" -v q="${MEM_QUAD_MB_PER_KRES2}" \
    'BEGIN { printf "%d", b + q * (L/1000.0)^2 }')"
BUDGET=$(( RAM_GB * 1024 ))
if (( PEAK > BUDGET )); then
    SHORT="$(awk -v a="$(( PEAK - BUDGET ))" 'BEGIN{printf "%.1f", a/1024}')"
    echo "[error] the construct projects to ${PEAK} MB against ${BUDGET} MB available," >&2
    echo "        short by ${SHORT} GB. Shorten it or raise the memory budget." >&2
    csc_record_exclusion "pedv_construct" "$STAGE" "" \
        "projected peak ${PEAK} MB exceeds the ${BUDGET} MB available at ${LEN} residues" \
        "the construct length is derived by 08a_pedv_analysis.py"
    csc_stage_end "refused for memory"
    exit 4
fi

# --- the alignment, once, from the same server as every other target ---------
MSA="${PEDV_DIR}/construct.a3m"
if [[ ! -s "$MSA" ]]; then
    WORK="$(csc_scratch "${STAGE}_msa")"
    mkdir -p "$WORK"
    csc_run "msa_pedv" "$(csc_env_bin "$CONDA_ENV_COLABFOLD" colabfold_batch)" \
        --msa-only --msa-mode "$MSA_MODE" --host-url "https://api.colabfold.com" \
        "$CONSTRUCT" "$WORK" || echo "[${STAGE}] the alignment could not be fetched"
    found="$(find "$WORK" -name '*.a3m' -size +0 | head -1)"
    [[ -n "$found" ]] && cp "$found" "$MSA"
    rm -rf "$WORK"
fi
[[ -s "$MSA" ]] && echo "[${STAGE}] alignment depth $(grep -c '^>' "$MSA")"

# --- predict ------------------------------------------------------------------
for arm in af2_msa_notmpl af2_nomsa; do
    out="${PEDV_DIR}/${arm}"
    if [[ -d "$out" && -n "$(find "$out" -name '*_scores_rank_001_*.json' 2>/dev/null)" && $FORCE -eq 0 ]]; then
        echo "[${STAGE}] ${arm} already predicted; skipping"
        continue
    fi
    mkdir -p "$out"
    args=(--model-type "$AF2_MODEL_TYPE" --num-models 1 --num-recycle "$AF2_NUM_RECYCLE"
          --random-seed "$SEED" --rank "$AF2_RANK" --data "${CACHE_DIR}/colabfold")
    input="$CONSTRUCT"
    case "$arm" in
        af2_nomsa) args+=(--msa-mode single_sequence) ;;
        *) [[ -s "$MSA" ]] && input="$MSA" ;;
    esac
    csc_run "predict_pedv_${arm}" \
        "$(csc_env_bin "$CONDA_ENV_COLABFOLD" colabfold_batch)" \
        "${args[@]}" "$input" "$out" || \
        echo "[${STAGE}] ${arm} failed; continuing"
done

# --- where the prediction put the domain --------------------------------------
csc_run "compare_prediction" "$PY_ANALYSIS" "${SCRIPT_DIR}/08a_pedv_analysis.py" \
    --config "${REPO_DIR}/project.conf" --entries "$ENTRIES" --compare

csc_stage_end "construct=${LEN} residues"
csc_mark_done "$STAGE"
