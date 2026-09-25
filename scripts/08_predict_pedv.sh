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
#  One of the two models compared in this repository was trained on structures
#  released before a date that precedes every deposited entry of this protein,
#  and the other on structures released after several of them. With templates
#  off, only the first is making a prediction in any useful sense; for the
#  second these coordinates were available during training. Both are run and
#  the table says which is which, because the comparison is more interesting
#  for being uneven.
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

csc_run "measure_deposited" "$PY_ANALYSIS" "${SCRIPT_DIR}/08a_pedv_analysis.py" \
    --config "${REPO_DIR}/project.conf" --entries "$ENTRIES" --measure

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
