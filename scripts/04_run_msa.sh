#!/usr/bin/env bash
# =============================================================================
#  04_run_msa.sh - fetch one alignment per target from the remote server
# =============================================================================
#  This stage is serial, and the reason is courtesy rather than memory.
#
#  The genetic databases a local search needs run to roughly a terabyte. They
#  do not fit on this machine and never will, so the alignments come from the
#  public server the prediction tool ships with. That server is free, shared,
#  and maintained by a small group. Its own instructions ask that queries be
#  made serially from a single address and not from several machines at once.
#  So this stage ignores JOBS, submits one query at a time, and waits between
#  them. A benchmark that finishes an hour sooner by hammering a free service
#  is not a benchmark worth running.
#
#  Everything is cached and every rerun skips what is already on disk. The
#  alignments are committed to the repository, compressed, because they are the
#  only artefact that makes a prediction reproducible: the server's databases
#  are versioned independently of this code, so the same query next year may
#  return a different alignment and therefore a different structure. Anyone
#  rerunning this repository from the committed alignments gets the inputs this
#  run used rather than whatever the server holds on the day.
#
#  What is recorded per target: the depth actually returned, the query date,
#  the server address, the database versions the server publishes, and the
#  queue depth at submission. The depth matters because it is the single
#  strongest predictor of how well the prediction will do, and a target whose
#  alignment came back shallow is a different kind of target from one whose
#  alignment is deep.
#
#  Usage:
#      bash scripts/04_run_msa.sh
#      bash scripts/04_run_msa.sh --limit 5
#      bash scripts/04_run_msa.sh --targets 1ABC_1,2DEF_1
#
#  Options:
#      --limit N       stop after N targets
#      --targets LIST  comma-separated target ids
#      --delay S       seconds between queries (default 5)
#      --force         refetch alignments already on disk
#      -h, --help      this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
csc_load_conf

STAGE=04_run_msa
LIMIT=0; TARGETS=""; DELAY=5; FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --limit)   LIMIT="$2"; shift 2 ;;
        --targets) TARGETS="$2"; shift 2 ;;
        --delay)   DELAY="$2"; shift 2 ;;
        --force)   FORCE=1; shift ;;
        -h|--help) sed -n '2,42p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

csc_skip_if_done "$STAGE" "$FORCE" && exit 0

PY_ANALYSIS="$(csc_env_bin "$CONDA_ENV_ANALYSIS" python)"
CF_BIN="$(csc_env_bin "$CONDA_ENV_COLABFOLD" colabfold_batch)"
[[ -n "$CF_BIN" ]] || { echo "[error] colabfold_batch not installed; run 01_install.sh" >&2; exit 3; }
MANIFEST="${CONFIG_DIR}/targets.tsv"
[[ -s "$MANIFEST" ]] || { echo "[error] ${MANIFEST} missing; run 02_build_holdout_set.py" >&2; exit 3; }

MSA_DIR="${DATA_DIR}/msa"
KEEP_DIR="${RESULTS_DIR}/msa"
mkdir -p "$MSA_DIR" "$KEEP_DIR"
DEPTHS="${RESULTS_DIR}/msa_depth.tsv"
[[ -s "$DEPTHS" ]] || printf 'target_id\tsequence_length\tdepth\tbytes\tserver\tqueue_at_submission\telapsed_s\tstatus\treason\trecorded\n' > "$DEPTHS"

WORK="$(csc_scratch "$STAGE")"
mkdir -p "$WORK"
cleanup() {
    local rc=$?
    rm -rf "$WORK"
    (( rc != 0 )) && echo "[${STAGE}] exited with status ${rc}; no stamp written" >&2
    return 0
}
trap cleanup EXIT INT TERM

csc_stage_start "$STAGE"

SERVER="https://api.colabfold.com"
# The server's own database versions, recorded once per run. They are published
# on a page rather than in the response, so this is the version of the record
# rather than of the answer, and the README says so.
DB_NOTE="$(curl -sS --max-time 20 "${SERVER}/queue" 2>/dev/null || echo "")"
{
    printf 'key\tvalue\trecorded\n'
    printf 'server\t%s\t%s\n' "$SERVER" "$(date -Iseconds)"
    printf 'queue_at_start\t%s\t%s\n' "${DB_NOTE:-unavailable}" "$(date -Iseconds)"
    printf 'query_date\t%s\t%s\n' "$(date -I)" "$(date -Iseconds)"
    printf 'msa_mode\t%s\t%s\n' "$MSA_MODE" "$(date -Iseconds)"
    printf 'policy\t%s\t%s\n' \
        "serial queries from one address, ${DELAY} s apart, as the server's instructions ask" \
        "$(date -Iseconds)"
} > "${RESULTS_DIR}/msa_server.tsv"

pick_targets() {
    if [[ -n "$TARGETS" ]]; then csc_split "$TARGETS"; return; fi
    awk -F'\t' '
        NR==1 { for (i = 1; i <= NF; i++) { if ($i == "sequence_length") len = i; if ($i == "target_id") id = i } next }
        len && id { print $len "\t" $id }' "$MANIFEST" | sort -n | cut -f2
}
mapfile -t TARGET_LIST < <(pick_targets)
(( LIMIT > 0 )) && TARGET_LIST=("${TARGET_LIST[@]:0:LIMIT}")

echo "[${STAGE}] ${#TARGET_LIST[@]} targets, serial, ${DELAY} s apart"
N_OK=0; N_SKIP=0; N_FAIL=0
for target in "${TARGET_LIST[@]}"; do
    [[ -n "$target" ]] || continue
    dest="${MSA_DIR}/${target}.a3m"
    kept="${KEEP_DIR}/${target}.a3m.gz"
    if [[ -s "$kept" && $FORCE -eq 0 ]]; then
        # Restore from the committed copy rather than asking the server again.
        # This is what makes the repository reproducible from a clean clone
        # without a single query.
        [[ -s "$dest" ]] || gunzip -c "$kept" > "$dest"
        echo "  ${target}: already have an alignment; skipping"
        N_SKIP=$(( N_SKIP + 1 ))
        continue
    fi

    queue="$(curl -sS --max-time 20 "${SERVER}/queue" 2>/dev/null || echo "")"
    out_dir="${WORK}/${target}"
    mkdir -p "$out_dir"
    "$PY_ANALYSIS" "${SCRIPT_DIR}/lib_predict.py" prepare \
        --config "${REPO_DIR}/project.conf" --target "$target" \
        --format fasta --out "${out_dir}/${target}.fasta"

    t0="$(date +%s)"
    rc=0
    # --msa-only stops after the alignment and predicts nothing, so the cost
    # here is the server's rather than this machine's.
    csc_run "msa_${target}" "$CF_BIN" --msa-only --msa-mode "$MSA_MODE" \
        --host-url "$SERVER" "${out_dir}/${target}.fasta" "$out_dir" || rc=$?
    t1="$(date +%s)"

    found="$(find "$out_dir" -name '*.a3m' -size +0 | head -1)"
    if (( rc != 0 )) || [[ -z "$found" ]]; then
        echo "  ${target}: no alignment returned (status ${rc})"
        printf '%s\t%s\t\t\t%s\t%s\t%d\tfailed\tthe server returned no alignment (status %s)\t%s\n' \
            "$target" "" "$SERVER" "${queue:-}" "$(( t1 - t0 ))" "$rc" "$(date -Iseconds)" >> "$DEPTHS"
        csc_record_exclusion "$target" "$STAGE" "" "no alignment returned by the server" "status ${rc}"
        N_FAIL=$(( N_FAIL + 1 ))
        rm -rf "$out_dir"
        sleep "$DELAY"
        continue
    fi

    cp "$found" "$dest"
    depth="$(grep -c '^>' "$dest" || echo 0)"
    bytes="$(stat -c %s "$dest")"
    # Committed compressed. These are text and compress by roughly ten to one.
    gzip -9 -c "$dest" > "$kept"
    len="$(awk -F'\t' -v t="$target" '
        NR==1 { for (i=1;i<=NF;i++) { if ($i=="sequence_length") l=i; if ($i=="target_id") id=i } next }
        $id == t { print $l; exit }' "$MANIFEST")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%d\tok\t\t%s\n' \
        "$target" "$len" "$depth" "$bytes" "$SERVER" "${queue:-}" "$(( t1 - t0 ))" "$(date -Iseconds)" >> "$DEPTHS"
    echo "  ${target}: depth ${depth}, $(( bytes / 1024 )) kB, $(( t1 - t0 )) s"
    N_OK=$(( N_OK + 1 ))
    rm -rf "$out_dir"

    # Between queries, not after the last one.
    if [[ "$target" != "${TARGET_LIST[-1]}" ]]; then sleep "$DELAY"; fi
done

# A size gate on what gets committed. Deep alignments of well studied families
# can be tens of megabytes; the repository has a 50 MB ceiling per file and
# anything approaching it is listed rather than tracked.
BIG="$(find "$KEEP_DIR" -name '*.a3m.gz' -size +40M -printf '%f %s\n' || true)"
if [[ -n "$BIG" ]]; then
    echo "[${STAGE}] alignments too large to track:"
    echo "$BIG" | while read -r f s; do
        echo "  ${f}: $(( s / 1048576 )) MB"
        csc_record_exclusion "${f%%.a3m.gz}" "$STAGE" "" \
            "compressed alignment is $(( s / 1048576 )) MB, above the 40 MB tracking limit" \
            "the alignment is on disk and the prediction used it; only the committed copy is omitted"
    done
fi

csc_stage_end "ok=${N_OK} skipped=${N_SKIP} failed=${N_FAIL}"
echo "[${STAGE}] ${N_OK} fetched, ${N_SKIP} already present, ${N_FAIL} failed"
csc_mark_done "$STAGE"
