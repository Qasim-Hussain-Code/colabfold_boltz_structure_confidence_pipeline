#!/usr/bin/env bash
# =============================================================================
#  05_predict.sh - one arm, one model, one target set
# =============================================================================
#  One invocation covers one arm. run_all.sh loops over the arms named in
#  project.conf. Keeping the loop outside means a failed arm can be rerun on
#  its own without redoing the rest, which matters when an arm takes a day.
#
#  Arms
#      af2_msa_notmpl  AlphaFold2 through ColabFold on the cached alignment,
#                      templates off. The primary arm.
#      af2_msa_tmpl    the same with templates on. Every template hit is
#                      recorded with its identity and its release date, because
#                      the template search runs against a current archive with
#                      no identity filter and can hand a held-out target its own
#                      structure. If it does, that target's prediction is a copy
#                      and the arm is measuring homology modelling.
#      af2_nomsa       the same weights and code with single-sequence mode. The
#                      only variable between this and the primary arm is the
#                      alignment, which is what makes the ablation clean.
#      boltz2_msa      Boltz-2 on the same alignment.
#
#  The null floors are built by 06_score_structures.py rather than here,
#  because copying a structure is not a prediction and pretending otherwise by
#  running it through this stage would invite the timing tables to include it.
#
#  Weights. One set on disk at a time: fetch, record the size and checksum,
#  run every target, delete. The cache is outside the repository and, on a
#  machine where the Linux filesystem is a virtual disk, outside that too: a
#  file written inside the virtual disk costs host disk permanently, while one
#  written to the mounted host drive is released when deleted. --keep-weights
#  overrides for a machine with room.
#
#  Memory. Every target is checked against the ceiling in project.conf before
#  it is run. A target above it is refused with its projected peak and the
#  shortfall in gigabytes, and the refusal is a row in results/excluded.tsv.
#  The alternative is an out-of-memory kill two hours in, which leaves nothing
#  to read afterwards.
#
#  Usage:
#      bash scripts/05_predict.sh --arm af2_msa_notmpl
#      bash scripts/05_predict.sh --arm boltz2_msa --limit 5
#      bash scripts/05_predict.sh --calibrate --out results/environment/memory_curve.tsv
#      bash scripts/05_predict.sh --arm af2_msa_notmpl --seed-variance
#
#  Options:
#      --arm NAME        one of the arms above
#      --limit N         stop after N targets
#      --targets LIST    comma-separated target ids, instead of the manifest
#      --num-models N    override the model count for this run
#      --seed-variance   rerun one target once per seed in SEED_REPLICATES
#      --calibrate       predict three short targets of increasing length with
#                        no alignment, record peak memory and elapsed time, and
#                        write the table 00_configure.sh --calibrate fits
#      --out FILE        where --calibrate writes its table
#      --keep-weights    do not delete the weights when the arm finishes
#      --force           ignore the stage stamp and redo
#      -h, --help        this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
csc_load_conf

ARM=""; LIMIT=0; TARGETS=""; NUM_MODELS=""; SEED_VAR=0; CALIBRATE=0
CAL_OUT=""; KEEP_WEIGHTS=0; FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arm)           ARM="$2"; shift 2 ;;
        --limit)         LIMIT="$2"; shift 2 ;;
        --targets)       TARGETS="$2"; shift 2 ;;
        --num-models)    NUM_MODELS="$2"; shift 2 ;;
        --seed-variance) SEED_VAR=1; shift ;;
        --calibrate)     CALIBRATE=1; shift ;;
        --out)           CAL_OUT="$2"; shift 2 ;;
        --keep-weights)  KEEP_WEIGHTS=1; shift ;;
        --force)         FORCE=1; shift ;;
        -h|--help)       sed -n '2,62p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

if (( CALIBRATE == 1 )); then
    ARM="af2_nomsa"
    [[ -n "$CAL_OUT" ]] || CAL_OUT="${RESULTS_DIR}/environment/memory_curve.tsv"
elif [[ -z "$ARM" ]]; then
    echo "[error] --arm is required" >&2; exit 1
fi

case "$ARM" in
    af2_msa_notmpl|af2_msa_tmpl|af2_nomsa) MODEL=alphafold2_ptm; SOURCE=colabfold ;;
    boltz2_msa)                            MODEL=boltz2;         SOURCE=boltz2 ;;
    *) echo "[error] ${ARM} is not an arm this stage runs" >&2; exit 1 ;;
esac

STAGE="05_predict_${ARM}"
(( SEED_VAR == 1 )) && STAGE="05_predict_seedvariance_${ARM}"
(( CALIBRATE == 1 )) && STAGE="05_predict_calibrate"
(( CALIBRATE == 0 )) && { csc_skip_if_done "$STAGE" "$FORCE" && exit 0; }

PY_ANALYSIS="$(csc_require_bin "$CONDA_ENV_ANALYSIS" python)" || exit 3
# The inference tool is resolved once, here, rather than inside the loop. A
# lookup that fails inside the loop produces a bare status of 127 on the first
# target, which says nothing about what was missing or where it was sought.
case "$SOURCE" in
    colabfold) PREDICT_BIN="$(csc_require_bin "$CONDA_ENV_COLABFOLD" colabfold_batch)" || exit 3 ;;
    boltz2)    PREDICT_BIN="$(csc_require_bin "$CONDA_ENV_BOLTZ" boltz)" || exit 3 ;;
esac
MANIFEST="${CONFIG_DIR}/targets.tsv"
[[ -s "$MANIFEST" ]] || { echo "[error] ${MANIFEST} missing; run 02_build_holdout_set.py" >&2; exit 3; }

# -----------------------------------------------------------------------------
# Weights: fetch, record, and delete when the arm is done.
# -----------------------------------------------------------------------------
# The parameter archive is fetched by the tool's own downloader rather than by
# a URL written here, so that the file and the code that reads it cannot come
# apart. config/sources.tsv records which archive that is, and the size and
# checksum of every file it leaves behind are recorded below.
AF2_DIR="${CACHE_DIR}/colabfold"
BOLTZ_DIR="${CACHE_DIR}/boltz"
WEIGHTS_TSV="${RESULTS_DIR}/environment/weights.tsv"
mkdir -p "$(dirname "$WEIGHTS_TSV")"
[[ -s "$WEIGHTS_TSV" ]] || printf 'model\tfile\tbytes\tsha256\tsource\tnote\trecorded\n' > "$WEIGHTS_TSV"

record_weight_file() {  # record_weight_file <model> <path> <note>
    local model="$1" f="$2" note="$3"
    [[ -f "$f" ]] || return 0
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$model" "$(basename "$f")" \
        "$(stat -c %s "$f")" "$(sha256sum "$f" | awk '{print $1}')" \
        "fetched by 05_predict.sh" "$note" "$(date -Iseconds)" >> "$WEIGHTS_TSV"
}

fetch_af2_weights() {
    if [[ -f "${AF2_DIR}/params/download_finished.txt" ]]; then
        echo "[${STAGE}] AlphaFold2 parameters already in the cache"
        return 0
    fi
    csc_gate_disk 3700 "the AlphaFold2 parameters" || return 4
    echo "[${STAGE}] fetching the AlphaFold2 parameters, about 3.5 GB"
    mkdir -p "${AF2_DIR}/params"
    # The download function is called with an explicit model type and directory
    # rather than through the module's own entry point. That entry point, given
    # no arguments, fetches the multimer parameters as well as these and puts
    # both in a default cache directory of its own choosing. On this machine
    # that meant 3.8 GB of weights for a model this benchmark does not run,
    # written inside the virtual disk where deleting them does not give the
    # space back. It got 1.6 GB in before it was stopped.
    # The downloader draws a progress bar on standard error, one line per
    # update, which writes tens of thousands of lines into a log that is meant
    # to be read afterwards. The bar is suppressed where the library honours
    # the setting and does no harm where it does not.
    export TQDM_DISABLE=1
    csc_run "fetch_af2_weights" "$(csc_env_bin "$CONDA_ENV_COLABFOLD" python)" -c \
        "from pathlib import Path; from colabfold.download import download_alphafold_params; download_alphafold_params('${AF2_MODEL_TYPE}', Path('${AF2_DIR}'))" \
        || { unset TQDM_DISABLE; echo "[error] the parameter download failed" >&2; return 4; }
    unset TQDM_DISABLE
    local f
    for f in "${AF2_DIR}"/params/*ptm*.npz; do
        record_weight_file alphafold2_ptm "$f" "one of five parameter sets for this model type"
    done
    # The archive holds ten parameter sets: five for this model type and five
    # for the original release, which this benchmark never loads. They are
    # extracted together and the unused half is 1.8 GB. Deleting it leaves the
    # success marker in place, so nothing re-fetches, and the five that remain
    # are the ones every arm here reads.
    local removed=0 bytes=0
    for f in "${AF2_DIR}"/params/params_model_*.npz; do
        case "$f" in *_ptm.npz) continue ;; esac
        [[ -f "$f" ]] || continue
        bytes=$(( bytes + $(stat -c %s "$f") ))
        rm -f "$f"
        removed=$(( removed + 1 ))
    done
    if (( removed > 0 )); then
        echo "[${STAGE}] removed ${removed} parameter sets this benchmark does not use, $(( bytes / 1048576 )) MB"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' alphafold2_ptm "params_model_1..5.npz" "$bytes" "" \
            "deleted by 05_predict.sh" \
            "the archive also carries five parameter sets for the original model type, which no arm here runs" \
            "$(date -Iseconds)" >> "$WEIGHTS_TSV"
    fi
    # The licence file inside the archive, recorded because it does not agree
    # with the statement in the model's current repository and the README says
    # so rather than choosing one of them.
    if [[ -f "${AF2_DIR}/params/LICENSE" ]]; then
        cp "${AF2_DIR}/params/LICENSE" "${RESULTS_DIR}/environment/alphafold_params_LICENSE.txt"
        head -3 "${AF2_DIR}/params/LICENSE" | tr '\n' ' ' \
            | sed "s|^|alphafold_params_license\t|" >> /dev/null
    fi
}

fetch_boltz_weights() {
    mkdir -p "$BOLTZ_DIR"
    if [[ -f "${BOLTZ_DIR}/boltz2_conf.ckpt" ]]; then
        echo "[${STAGE}] Boltz-2 weights already in the cache"
        return 0
    fi
    # The peak is the moment the molecule archive and its extracted contents
    # are both on disk, which is larger than what is left afterwards.
    csc_gate_disk 6100 "the Boltz-2 weights" || return 4
    # The affinity checkpoint is fetched by the tool whenever the file is
    # absent, whether or not affinity is requested, and nothing here requests
    # it. The check is for existence alone, so an empty placeholder skips a two
    # gigabyte download of a model that would never be loaded. This is a
    # deliberate pre-population of the cache, it is recorded in the weights
    # table as such, and it is the difference between this arm fitting on this
    # machine and not.
    if [[ ! -f "${BOLTZ_DIR}/boltz2_aff.ckpt" ]]; then
        : > "${BOLTZ_DIR}/boltz2_aff.ckpt"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' boltz2 boltz2_aff.ckpt 0 "" \
            "placeholder written by 05_predict.sh" \
            "empty on purpose: the affinity model is never loaded here and the downloader checks only that the file exists" \
            "$(date -Iseconds)" >> "$WEIGHTS_TSV"
    fi
    echo "[${STAGE}] fetching the Boltz-2 weights, about 4.4 GB"
    # The tool downloads on its first run rather than offering a download
    # command, so a trivial prediction is used to trigger it. The prediction
    # itself is thrown away.
    local probe="${DATA_DIR}/work/boltz_weight_probe"
    rm -rf "$probe"; mkdir -p "$probe"
    printf 'version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: MKVLIS\n      msa: empty\n' \
        > "${probe}/probe.yaml"
    csc_run "fetch_boltz_weights" "$(csc_env_bin "$CONDA_ENV_BOLTZ" boltz)" predict \
        "${probe}/probe.yaml" --out_dir "$probe" --cache "$BOLTZ_DIR" \
        --accelerator cpu --diffusion_samples 1 --sampling_steps 10 \
        --recycling_steps 0 --num_workers 0 --seed 1 || true
    rm -rf "$probe"
    [[ -f "${BOLTZ_DIR}/boltz2_conf.ckpt" ]] || { echo "[error] the Boltz-2 weights did not arrive" >&2; return 4; }
    record_weight_file boltz2 "${BOLTZ_DIR}/boltz2_conf.ckpt" "structure and confidence model"
    # The molecule archive is kept next to its own extracted contents, which
    # doubles its cost for nothing. Emptying the archive keeps the existence
    # check satisfied without keeping the bytes.
    if [[ -s "${BOLTZ_DIR}/mols.tar" && -d "${BOLTZ_DIR}/mols" ]]; then
        record_weight_file boltz2 "${BOLTZ_DIR}/mols.tar" "component dictionary archive, emptied after extraction"
        : > "${BOLTZ_DIR}/mols.tar"
        echo "[${STAGE}] emptied the extracted molecule archive"
    fi
}

delete_weights() {
    (( KEEP_WEIGHTS == 1 )) && { echo "[${STAGE}] keeping the weights (--keep-weights)"; return 0; }
    case "$MODEL" in
        alphafold2_ptm) rm -rf "${AF2_DIR}/params"; echo "[${STAGE}] deleted the AlphaFold2 parameters" ;;
        boltz2)         rm -rf "${BOLTZ_DIR}"; echo "[${STAGE}] deleted the Boltz-2 weights" ;;
    esac
}

# -----------------------------------------------------------------------------
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

# Which targets, and in which order. Shortest first, so that a run cut short by
# time still covers the cheap end of the range rather than stopping part way
# through one long target.
pick_targets() {
    if [[ -n "$TARGETS" ]]; then csc_split "$TARGETS"; return; fi
    # Column 5 is the sequence length and column 1 the identifier. The columns
    # are found by name rather than by position, because a manifest that gains
    # a column would otherwise silently sort by something else: the first
    # version of this sorted on the pre-cutoff identity column and produced a
    # plausible-looking order that had nothing to do with length.
    awk -F'\t' '
        NR==1 { for (i = 1; i <= NF; i++) { if ($i == "sequence_length") len = i; if ($i == "target_id") id = i } next }
        len && id { print $len "\t" $id }' "$MANIFEST" | sort -n | cut -f2
}
mapfile -t TARGET_LIST < <(pick_targets)
(( LIMIT > 0 )) && TARGET_LIST=("${TARGET_LIST[@]:0:LIMIT}")

if (( CALIBRATE == 1 )); then
    # Three targets spread across the length range of the set, so the fit has
    # something to fit rather than three points on top of each other.
    mapfile -t TARGET_LIST < <(
        pick_targets | awk '{a[NR]=$0} END {print a[1]; print a[int((NR+1)/2)]; print a[NR]}'
    )
    echo "[${STAGE}] calibrating on ${TARGET_LIST[*]}"
fi

seq_len_of() {
    awk -F'\t' -v t="$1" '
        NR==1 { for (i = 1; i <= NF; i++) { if ($i == "sequence_length") len = i; if ($i == "target_id") id = i } next }
        $id == t { print $len; exit }' "$MANIFEST"
}

# The memory gate. Refuse rather than let the kernel kill the run.
projected_peak_mb() {
    awk -v L="$1" -v b="${MEM_BASE_MB}" -v q="${MEM_QUAD_MB_PER_KRES2}" \
        'BEGIN { printf "%d", b + q * (L/1000.0)^2 }'
}

# SEED, SEED_REPLICATES and the protocol settings come from project.conf,
# which csc_load_conf sourced above and which fails loudly if any of them is
# missing. The scripts directory is taken from this file's own location rather
# than from the config, so a stale path in project.conf cannot send this stage
# to a different copy of the library it was written against.
# shellcheck disable=SC2153
SEEDS=("$SEED")
# shellcheck disable=SC2153
(( SEED_VAR == 1 )) && read -r -a SEEDS <<< "$SEED_REPLICATES"
N_MODELS="${NUM_MODELS:-$AF2_NUM_MODELS}"
(( CALIBRATE == 1 )) && N_MODELS=1

case "$MODEL" in
    alphafold2_ptm) fetch_af2_weights || exit 4 ;;
    boltz2)         fetch_boltz_weights || exit 4 ;;
esac

# The calibration table is started only once the weights are in hand. An
# earlier version wrote its header before the disk gate ran, so a refusal for
# want of four hundred megabytes destroyed a set of measurements that had cost
# half an hour to take. They were recovered from the per-command log, which is
# the other reason that log exists, and the header moved here.
if (( CALIBRATE == 1 )) && [[ ! -s "$CAL_OUT" ]]; then
    printf 'target_id\tsequence_length\tpeak_rss_mb\telapsed_s\tarm\trecorded\n' > "$CAL_OUT"
fi

N_OK=0; N_SKIP=0; N_REFUSED=0; N_FAIL=0
for target in "${TARGET_LIST[@]}"; do
    [[ -n "$target" ]] || continue
    len="$(seq_len_of "$target")"
    [[ -n "$len" ]] || { echo "  ${target}: not in the manifest, skipped"; continue; }

    peak="$(projected_peak_mb "$len")"
    budget=$(( RAM_GB * 1024 ))
    if (( peak > budget )); then
        short_gb="$(awk -v a="$(( peak - budget ))" 'BEGIN{printf "%.1f", a/1024}')"
        echo "  ${target}: refused, ${len} residues projects to ${peak} MB against ${budget} MB, short by ${short_gb} GB"
        # Written through the shared helper rather than appended by hand, so
        # the table gets its header on first use and the column order cannot
        # drift away from every other writer of this file.
        csc_record_exclusion "$target" "05_predict" "$ARM" \
            "projected peak ${peak} MB exceeds the ${budget} MB available at ${len} residues" \
            "short by ${short_gb} GB; the ceiling is MAX_SEQ_LEN in project.conf"
        N_REFUSED=$(( N_REFUSED + 1 ))
        continue
    fi

    for seed in "${SEEDS[@]}"; do
        tag="$target"
        (( SEED_VAR == 1 )) && tag="${target}__seed${seed}"
        out_dir="${WORK}/${tag}"
        done_marker="${RESULTS_DIR}/confidence/${target}__${ARM}.json.gz"
        (( SEED_VAR == 1 )) && done_marker="${RESULTS_DIR}/confidence/${tag}__${ARM}.json.gz"
        if [[ -f "$done_marker" && $FORCE -eq 0 ]]; then
            echo "  ${tag}: already predicted; skipping"
            N_SKIP=$(( N_SKIP + 1 ))
            continue
        fi
        mkdir -p "$out_dir"
        msa_file="${DATA_DIR}/msa/${target}.a3m"

        rc=0
        if [[ "$SOURCE" == "colabfold" ]]; then
            input="${out_dir}/${target}.fasta"
            cf_args=(--model-type "$AF2_MODEL_TYPE" --num-models "$N_MODELS"
                     --num-recycle "$AF2_NUM_RECYCLE" --random-seed "$seed"
                     --rank "$AF2_RANK" --data "$AF2_DIR")
            case "$ARM" in
                af2_nomsa)
                    # FASTA input, never a cached alignment: single-sequence
                    # mode reads an alignment if it is handed one.
                    "$PY_ANALYSIS" "${SCRIPT_DIR}/lib_predict.py" prepare \
                        --config "${REPO_DIR}/project.conf" --target "$target" \
                        --format fasta --out "$input"
                    cf_args+=(--msa-mode single_sequence)
                    msa_file=""
                    ;;
                af2_msa_tmpl)
                    input="$msa_file"
                    cf_args+=(--templates)
                    ;;
                *)
                    input="$msa_file"
                    ;;
            esac
            if [[ "$ARM" != "af2_nomsa" && ! -s "$input" ]]; then
                echo "  ${tag}: no alignment at ${input}; run 04_run_msa.sh first"
                N_FAIL=$(( N_FAIL + 1 ))
                continue
            fi
            csc_run "predict_${tag}" "$PREDICT_BIN" \
                "${cf_args[@]}" "$input" "$out_dir" || rc=$?
        else
            input="${out_dir}/${target}.yaml"
            msa_arg=""
            [[ -s "$msa_file" ]] && msa_arg="$msa_file"
            "$PY_ANALYSIS" "${SCRIPT_DIR}/lib_predict.py" prepare \
                --config "${REPO_DIR}/project.conf" --target "$target" \
                --format yaml --msa "$msa_arg" --out "$input"
            csc_run "predict_${tag}" "$PREDICT_BIN" predict \
                "$input" --out_dir "$out_dir" --cache "$BOLTZ_DIR" \
                --accelerator cpu --seed "$seed" \
                --recycling_steps "$BOLTZ_RECYCLING_STEPS" \
                --sampling_steps "$BOLTZ_SAMPLING_STEPS" \
                --diffusion_samples "$BOLTZ_DIFFUSION_SAMPLES" \
                --step_scale "$BOLTZ_STEP_SCALE" \
                --output_format mmcif --num_workers 0 || rc=$?
        fi

        # The elapsed time and peak memory this target actually used, read back
        # from the per-command table csc_run just appended to.
        read -r elapsed peak_kb < <(
            awk -F'\t' -v l="predict_${tag}" '$2==l {e=$3; p=$4} END {print e, p}' \
                "${LOG_DIR}/${STAGE}.commands.tsv")

        collect_dir="$out_dir"
        [[ "$SOURCE" == "boltz2" ]] && collect_dir="$(find "$out_dir" -maxdepth 1 -type d -name 'boltz_results_*' | head -1)"
        [[ -n "$collect_dir" ]] || collect_dir="$out_dir"
        if ! "$PY_ANALYSIS" "${SCRIPT_DIR}/lib_predict.py" collect \
                --config "${REPO_DIR}/project.conf" --target "$target" --arm "$ARM" \
                --model "$MODEL" --seed "$seed" --source "$SOURCE" \
                --out-dir "$collect_dir" --elapsed "${elapsed:-}" \
                --peak-rss-kb "${peak_kb:-}" --msa-file "$msa_file" \
                --n-models "$N_MODELS"; then
            N_FAIL=$(( N_FAIL + 1 ))
        else
            N_OK=$(( N_OK + 1 ))
        fi

        if (( CALIBRATE == 1 )); then
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$target" "$len" \
                "$(awk -v k="${peak_kb:-0}" 'BEGIN{printf "%.0f", k/1024}')" \
                "${elapsed:-}" "$ARM" "$(date -Iseconds)" >> "$CAL_OUT"
        fi
        rm -rf "$out_dir"
    done
done

delete_weights
csc_stage_end "arm=${ARM} ok=${N_OK} skipped=${N_SKIP} refused=${N_REFUSED} failed=${N_FAIL}"
echo "[${STAGE}] ${N_OK} predicted, ${N_SKIP} already present, ${N_REFUSED} refused for memory, ${N_FAIL} failed"
(( CALIBRATE == 1 )) && { echo "[${STAGE}] wrote ${CAL_OUT}"; exit 0; }
# A run that covered only part of the set must not mark the stage finished.
# A single-target test did exactly that here, and the full run that followed
# skipped the arm and reported success without predicting anything.
if (( LIMIT > 0 )) || [[ -n "$TARGETS" ]]; then
    echo "[${STAGE}] this run covered part of the set, so the stage is not marked"
    echo "           as finished. Re-run without --limit or --targets to complete it."
else
    csc_mark_done "$STAGE"
fi
