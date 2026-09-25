#!/usr/bin/env bash
# =============================================================================
#  lib_common.sh - sourced by every stage. Not executable on its own.
# =============================================================================
#  Five jobs:
#    1. find and source project.conf, failing loudly if 00_configure.sh has not
#       run, because a stage that silently falls back to shell defaults for
#       THREADS or SEED produces results nobody can reproduce
#    2. activate one of the four conda environments without a login shell
#    3. record elapsed time, peak resident memory and disk use per command and
#       per stage into logs/, so the README quotes measurements
#    4. gate on disk: the data directory, the weight cache, and under WSL the
#       Windows drive that actually holds the Linux filesystem
#    5. stage stamps, so re-running a finished stage skips and says so
#
#  The layout follows stage 1 (vina_gnina_pose_benchmark_pipeline) function for
#  function, with the prefix changed from vgb_ to csc_. Two things are new.
#  Stage 1 could not see the Windows host drive from inside WSL and said so in
#  its data/README.md; this library measures it. And stage 1 had one conda
#  environment to activate, where this stage has four, because ColabFold's JAX
#  pin and Boltz's PyTorch pin cannot share one.
# =============================================================================

# shellcheck disable=SC2148
# (this file is sourced; the shebang is for editors and shellcheck)

CSC_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CSC_REPO_DIR="$(cd "${CSC_LIB_DIR}/.." && pwd)"

# ---- 1. project.conf --------------------------------------------------------
csc_load_conf() {
    local conf="${CSC_REPO_DIR}/project.conf"
    if [[ ! -r "$conf" ]]; then
        echo "[error] project.conf not found. Run scripts/00_configure.sh first." >&2
        echo "        Every stage reads its thread count, seed, memory ceiling and" >&2
        echo "        paths from there; nothing here has a silent default." >&2
        return 3
    fi
    # shellcheck source=/dev/null
    source "$conf"
    local v
    for v in THREADS RAM_GB DISK_GB JOBS SEED DATA_DIR CACHE_DIR RESULTS_DIR \
             LOG_DIR CONFIG_DIR SCRIPTS_DIR FIGURES_DIR CONDA_SH; do
        [[ -n "${!v:-}" ]] || { echo "[error] ${v} missing from project.conf; re-run 00_configure.sh" >&2; return 3; }
    done
    mkdir -p "$DATA_DIR" "$CACHE_DIR" "$RESULTS_DIR" "$FIGURES_DIR" "$LOG_DIR" "$CONFIG_DIR"
}

# ---- 2. conda ---------------------------------------------------------------
# conda.sh is sourced by hand because a non-interactive script has no login
# shell. `conda run` was the alternative; stage 1 found that it buffers stdout
# until the child exits, which makes an hour-long prediction look hung.
csc_conda_init() {
    [[ -n "${CONDA_SH:-}" && -r "$CONDA_SH" ]] || {
        echo "[error] CONDA_SH='${CONDA_SH:-}' is not readable; re-run 00_configure.sh" >&2; return 3; }
    set +u
    # shellcheck source=/dev/null
    source "$CONDA_SH"
    set -u
}

# csc_activate <env_name>
csc_activate() {
    csc_conda_init || return 3
    set +u
    conda activate "$1" || { set -u; echo "[error] conda env '${1}' not found. Run scripts/01_install.sh." >&2; return 3; }
    set -u
}

# csc_conda_root - the conda installation prefix, from CONDA_SH
csc_conda_root() { dirname "$(dirname "$(dirname "$CONDA_SH")")"; }

# csc_env_prefix <env_name> - absolute prefix of an environment, or empty
csc_env_prefix() {
    local p
    p="$(csc_conda_root)/envs/${1}"
    [[ -d "$p" ]] && echo "$p" || echo ""
}

# csc_env_bin <env_name> <binary>
# Absolute path to a binary inside an environment. The prediction stages call
# ColabFold and Boltz by absolute path rather than activating their
# environments, for the reason stage 1 called Vina and smina that way: the
# orchestration python lives in the analysis environment, and activating a
# second environment mid-loop would put its libstdc++ and its numpy first.
csc_env_bin() {
    local prefix; prefix="$(csc_env_prefix "$1")"
    [[ -n "$prefix" && -x "${prefix}/bin/${2}" ]] && { echo "${prefix}/bin/${2}"; return 0; }
    echo ""; return 1
}

# csc_require_bin <env_name> <binary>
# The same, but refuses rather than returning an empty string.
#
# The empty string is the dangerous case. Passed as the command to run, it
# produces a bare "command not found" and a status of 127, with nothing in the
# log to say which tool was missing or where it was looked for. That happened
# here once, on a run where the filesystem holding the environments was under
# heavy write pressure, and cost an arm that had already paid for its weights.
csc_require_bin() {
    local p; p="$(csc_env_bin "$1" "$2" || true)"
    if [[ -z "$p" ]]; then
        echo "[error] ${2} not found in the conda environment '${1}'." >&2
        echo "        Looked in $(csc_env_prefix "$1")/bin. Run scripts/01_install.sh," >&2
        echo "        or check that CONDA_SH in project.conf points at the right conda." >&2
        return 3
    fi
    echo "$p"
}

# ---- 3. measurement ---------------------------------------------------------
# Peak RSS comes from GNU time's %M, in kilobytes. The shell builtin `time`
# cannot report memory, so a missing /usr/bin/time is reported as NA rather
# than replaced by a guess.
csc_time_bin() {
    local t
    for t in /usr/bin/time /bin/time; do
        [[ -x "$t" ]] && { echo "$t"; return; }
    done
    echo ""
}

# csc_dir_kb <dir> - apparent size in kilobytes, 0 if absent
csc_dir_kb() {
    [[ -d "$1" ]] || { echo 0; return; }
    du -sk "$1" 2>/dev/null | awk '{print $1}' || echo 0
}

# csc_free_mb <path> - free megabytes on the filesystem holding <path>
csc_free_mb() {
    df -Pm "$1" 2>/dev/null | awk 'NR==2 {print $4}' || echo 0
}

# csc_is_wsl - true under WSL1 or WSL2
csc_is_wsl() { grep -qi microsoft /proc/version 2>/dev/null; }

# csc_host_free_mb - free megabytes on the Windows drive behind the Linux
# filesystem, or empty when not under WSL.
#
# This is the number stage 1 could not see. Under WSL2 the Linux root is a
# thin-provisioned ext4.vhdx on a Windows drive. `df /` inside the virtual
# machine reported 894 GB free while the host drive fell below 2 GB, and on
# the day this stage was started the host drive had 30 MB left. The vhdx grows
# one-for-one with every byte written inside it and does not shrink when a
# file is deleted, so the only honest gate under WSL is the host drive's own
# free space. HOST_DRIVE_MOUNT in project.conf names it; /mnt/c by default.
csc_host_free_mb() {
    csc_is_wsl || { echo ""; return; }
    local m="${HOST_DRIVE_MOUNT:-/mnt/c}"
    [[ -d "$m" ]] || { echo ""; return; }
    csc_free_mb "$m"
}

# csc_stage_start <stage>
csc_stage_start() {
    CSC_STAGE="$1"
    CSC_T0="$(date +%s)"
    CSC_DISK0="$(csc_dir_kb "$DATA_DIR")"
    CSC_CACHE0="$(csc_dir_kb "$CACHE_DIR")"
    CSC_HOST0="$(csc_host_free_mb)"
    CSC_PEAK_RSS_KB=0
    CSC_PEAK_DISK_KB=$(( CSC_DISK0 + CSC_CACHE0 ))
    CSC_MIN_HOST_MB="${CSC_HOST0:-}"
    echo "[${CSC_STAGE}] started $(date -Iseconds)"
}

# csc_sample_disk - update the stage's peak disk and minimum host free space.
# Called after every command. A weight set fetched and deleted inside one
# command would be missed, which is why fetch, run and delete are always
# separate csc_run calls in 05_predict.sh.
csc_sample_disk() {
    local d c h
    d="$(csc_dir_kb "$DATA_DIR")"
    c="$(csc_dir_kb "$CACHE_DIR")"
    (( d + c > CSC_PEAK_DISK_KB )) && CSC_PEAK_DISK_KB=$(( d + c ))
    h="$(csc_host_free_mb)"
    if [[ -n "$h" ]]; then
        if [[ -z "${CSC_MIN_HOST_MB:-}" ]] || (( h < CSC_MIN_HOST_MB )); then
            CSC_MIN_HOST_MB="$h"
        fi
    fi
}

# csc_run <label> <command...>
# Runs a command under GNU time, appends its elapsed seconds and peak RSS to the
# stage's per-command table, and keeps the largest RSS seen so the stage's peak
# is the peak of its children rather than of the shell.
csc_run() {
    local label="$1"; shift
    local tb; tb="$(csc_time_bin)"
    local tf; tf="$(mktemp)"
    local rc=0 elapsed=NA peak=NA
    if [[ -n "$tb" ]]; then
        "$tb" -f '%e %M' -o "$tf" "$@" || rc=$?
        # GNU time prints "Command exited with non-zero status N" on its own
        # line before the format line when the child fails; take the last line.
        read -r elapsed peak < <(tail -n 1 "$tf") || { elapsed=NA; peak=NA; }
        if [[ "$peak" =~ ^[0-9]+$ ]] && (( peak > CSC_PEAK_RSS_KB )); then
            CSC_PEAK_RSS_KB="$peak"
        fi
    else
        local t0 t1; t0="$(date +%s)"
        "$@" || rc=$?
        t1="$(date +%s)"
        elapsed=$(( t1 - t0 ))
    fi
    local out="${LOG_DIR}/${CSC_STAGE}.commands.tsv"
    [[ -s "$out" ]] || printf 'stage\tlabel\telapsed_s\tpeak_rss_kb\texit_status\trecorded\n' > "$out"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$CSC_STAGE" "$label" "$elapsed" "$peak" "$rc" \
        "$(date -Iseconds)" >> "$out"
    rm -f "$tf"
    csc_sample_disk
    return $rc
}

# csc_stage_end [note]
csc_stage_end() {
    local note="${1:-ok}"
    local t1 disk1 cache1 elapsed peak
    t1="$(date +%s)"
    disk1="$(csc_dir_kb "$DATA_DIR")"
    cache1="$(csc_dir_kb "$CACHE_DIR")"
    csc_sample_disk
    elapsed=$(( t1 - CSC_T0 ))
    peak="$CSC_PEAK_RSS_KB"
    local out="${LOG_DIR}/${CSC_STAGE}.resources.tsv"
    printf 'stage\telapsed_s\tpeak_rss_mb\tdata_dir_growth_mb\tdata_dir_total_mb\tcache_dir_total_mb\tpeak_disk_mb\thost_free_start_mb\thost_free_min_mb\tjobs\tnote\n' > "$out"
    printf '%s\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' \
        "$CSC_STAGE" "$elapsed" \
        "$( [[ "$peak" == "0" ]] && echo NA || echo $(( peak / 1024 )) )" \
        "$(( (disk1 - CSC_DISK0) / 1024 ))" "$(( disk1 / 1024 ))" "$(( cache1 / 1024 ))" \
        "$(( CSC_PEAK_DISK_KB / 1024 ))" "${CSC_HOST0:-NA}" "${CSC_MIN_HOST_MB:-NA}" \
        "${JOBS}" "$note" >> "$out"
    echo "[${CSC_STAGE}] finished in ${elapsed} s; data dir $(( disk1 / 1024 )) MB, cache $(( cache1 / 1024 )) MB, peak disk $(( CSC_PEAK_DISK_KB / 1024 )) MB"
}

# ---- 4. disk gate -----------------------------------------------------------
# csc_gate_disk <projected_additional_mb> <what>
# Refuses, with the shortfall in gigabytes, if the projected addition will not
# fit in the DISK_GB budget, on the filesystems it lands on, or on the Windows
# host drive under WSL. HOST_RESERVE_MB is kept free on the host drive for the
# WSL swap file, which also lives there and grows under memory pressure.
csc_gate_disk() {
    local need_mb="$1" what="$2"
    local used_mb=$(( ( $(csc_dir_kb "$DATA_DIR") + $(csc_dir_kb "$CACHE_DIR") ) / 1024 ))
    local budget_mb=$(( DISK_GB * 1024 ))
    local fail=0
    if (( used_mb + need_mb > budget_mb )); then
        printf '[error] %s needs about %d MB more; %d MB already used of the %d GB budget.\n' \
            "$what" "$need_mb" "$used_mb" "$DISK_GB" >&2
        printf '        Shortfall %.1f GB. Raise --disk in 00_configure.sh or cut the set.\n' \
            "$(awk -v a="$(( used_mb + need_mb - budget_mb ))" 'BEGIN{print a/1024}')" >&2
        fail=1
    fi
    local host; host="$(csc_host_free_mb)"
    local reserve="${HOST_RESERVE_MB:-3072}"
    if [[ -n "$host" ]] && (( need_mb + reserve > host )); then
        printf '[error] %s needs about %d MB and the Windows drive behind WSL has %d MB free,\n' \
            "$what" "$need_mb" "$host" >&2
        printf '        of which %d MB is held back for the WSL swap file. Shortfall %.1f GB.\n' \
            "$reserve" "$(awk -v a="$(( need_mb + reserve - host ))" 'BEGIN{print a/1024}')" >&2
        echo "        df inside WSL cannot see this; see data/README.md." >&2
        fail=1
    fi
    return $fail
}

# ---- 5. idempotency ---------------------------------------------------------
# One stamp file per stage. A stage checks its stamp and skips, which makes
# run_all.sh restartable after a failure without redoing hours of prediction.
# --force on any stage removes its own stamp.
csc_stamp() { echo "${LOG_DIR}/${1}.done"; }
csc_is_done() { [[ -f "$(csc_stamp "$1")" ]]; }
csc_mark_done() { date -Iseconds > "$(csc_stamp "$1")"; }
csc_skip_if_done() {  # csc_skip_if_done <stage> <force 0|1>
    local stage="$1" force="${2:-0}"
    if [[ "$force" == "1" ]]; then rm -f "$(csc_stamp "$stage")"; return 1; fi
    if csc_is_done "$stage"; then
        echo "[${stage}] already completed on $(cat "$(csc_stamp "$stage")"); skipping (--force to redo)."
        return 0
    fi
    return 1
}

# csc_record_exclusion <target> <stage> <arm> <reason> [detail]
# One row per dropped target, in the same table and the same column order the
# Python stages write. Silent exclusions are how benchmark numbers get
# inflated, so every refusal lands here whichever language refused it.
csc_record_exclusion() {
    local f="${RESULTS_DIR}/excluded.tsv"
    mkdir -p "$RESULTS_DIR"
    [[ -s "$f" ]] || printf 'target_id\tstage\tarm\treason\tdetail\trecorded\n' > "$f"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "${5:-}" "$(date -Iseconds)" >> "$f"
}

# ---- misc -------------------------------------------------------------------
csc_need() {
    command -v "$1" >/dev/null 2>&1 || { echo "[error] ${1} not found on PATH" >&2; return 3; }
}
# Comma-separated list to newline-separated, dropping blanks.
csc_split() { tr ',' '\n' <<<"$1" | sed '/^$/d'; }
# Is <needle> in the comma-separated <haystack>?
csc_in_list() { case ",${2}," in *",${1},"*) return 0 ;; *) return 1 ;; esac; }

# csc_scratch <name> - a per-run scratch directory, in shared memory when there
# is room. Stage 1 measured about 3.5 MB of permanent host-drive growth per
# docking run when scratch lived inside the vhdx, because the virtual disk
# never gives space back. /dev/shm never touches it. CSC_WORK_DIR overrides.
csc_scratch() {
    local name="$1" shm_mb
    if [[ -n "${CSC_WORK_DIR:-}" ]]; then echo "${CSC_WORK_DIR}/${name}"; return; fi
    shm_mb="$(csc_free_mb /dev/shm)"
    if [[ -w /dev/shm && -n "$shm_mb" ]] && (( shm_mb > 1024 )); then
        echo "/dev/shm/csc_work/${name}"
    else
        echo "${DATA_DIR}/work/${name}"
    fi
}
