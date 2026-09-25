#!/usr/bin/env bash
# =============================================================================
#  00_configure.sh - measure the machine, project the cost, refuse if it will
#  not fit, write project.conf
# =============================================================================
#  Every later script sources project.conf rather than hardcoding a thread
#  count, a seed, a sequence-length ceiling or a path. Re-run this script to
#  change any of them; it overwrites project.conf and nothing else.
#
#  Two passes, and the difference between them matters.
#
#  The first pass runs before anything is installed. It can measure threads,
#  memory, free disk and the absence of a GPU, but it cannot measure how much
#  memory a prediction takes, because there is nothing to run. So it writes a
#  provisional ceiling derived from published attention scaling and says in its
#  own output that the figure is a projection rather than evidence.
#
#  The second pass is --calibrate, run by run_all.sh after 01_install.sh. It
#  predicts three short targets of increasing length, records peak resident
#  memory and elapsed time for each, fits a quadratic in sequence length, and
#  rewrites MAX_SEQ_LEN and the timing projection from the fit. Stage 1 made
#  the same distinction between a projection from defaults and a gate that has
#  measured the machine, and the gate that measured is the one to trust.
#
#  Why a quadratic. Attention over residues is quadratic in the number of
#  residues, and the Evoformer pair representation is an L by L tensor, so the
#  dominant term in both memory and time grows as L squared. Three points is
#  the minimum that can fit a + b*L^2 with a residual, and the residual is
#  reported so a reader can see whether the fit is worth anything. It is not a
#  law, it is a local approximation over the range measured, and extrapolating
#  it four times past the longest measured target is the least defensible
#  number this script produces.
#
#  Usage:
#      bash scripts/00_configure.sh --threads 16 --ram 7 --disk 25 \
#          --cache-dir /mnt/c/csc_cache --data-dir ~/csc_data --yes
#      bash scripts/00_configure.sh --calibrate          # after 01_install.sh
#
#  Options:
#      --threads N      CPU threads the machine may use (default: detected)
#      --ram GB         RAM in gigabytes the analysis may use (default:
#                       detected, minus a gigabyte for the operating system)
#      --disk GB        budget for DATA_DIR plus CACHE_DIR (default: measured
#                       free space, minus a margin; under WSL the Windows host
#                       drive is measured instead, see data/README.md)
#      --cache-dir DIR  model weights. Gigabytes, fetched and deleted one set
#                       at a time. Keep this outside the repository and, under
#                       WSL, outside the Linux filesystem: a file written
#                       inside the virtual disk consumes host disk permanently,
#                       while one written to /mnt/c is released when deleted.
#      --data-dir DIR   alignments, predictions, scores. Small by comparison.
#      --jobs N         concurrent prediction processes (default 1). Serial is
#                       not politeness here, it is arithmetic: CPU inference
#                       already uses every core, so two jobs halve the cores
#                       each gets and double the memory. The MSA stage is
#                       serial whatever this says, for a different reason given
#                       in 04_run_msa.sh.
#      --hours N        wall-clock budget for the whole run. The script
#                       projects the total from the measured curve and refuses
#                       without --yes if the projection exceeds it.
#      --arms LIST      comma-separated subset of the arms in ARMS below
#      --models LIST    comma-separated subset of alphafold2_ptm,boltz2
#      --max-targets N  cap on the held-out set size (0 means the set size the
#                       filters produce)
#      --calibrate      second pass: measure the memory and time curve on this
#                       machine and rewrite the ceiling. Needs 01_install.sh.
#      --yes, -y        accept every default and proceed past the time gate
#      -h, --help       this text
# =============================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
CONF="${REPO_DIR}/project.conf"
CONF_TMP="${CONF}.tmp.$$"

# A killed run must not leave a half-written project.conf behind, because the
# next script would source it and inherit a truncated variable list.
cleanup() { rm -f "$CONF_TMP"; }
trap cleanup EXIT INT TERM

THREADS=""; RAM_GB=""; DISK_GB=""; JOBS=""; CACHE_DIR=""; DATA_DIR=""
HOURS=""; ASSUME_YES=0; CALIBRATE=0; MAX_TARGETS=""
DEFAULT_ARMS="af2_msa_notmpl,af2_msa_tmpl,af2_nomsa,boltz2_msa,null_template,null_unrelated"
DEFAULT_MODELS="alphafold2_ptm,boltz2"
ARMS="$DEFAULT_ARMS"
MODELS="$DEFAULT_MODELS"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)     THREADS="$2"; shift 2 ;;
        --ram)         RAM_GB="$2"; shift 2 ;;
        --disk)        DISK_GB="$2"; shift 2 ;;
        --jobs)        JOBS="$2"; shift 2 ;;
        --cache-dir)   CACHE_DIR="$2"; shift 2 ;;
        --data-dir)    DATA_DIR="$2"; shift 2 ;;
        --hours)       HOURS="$2"; shift 2 ;;
        --arms)        ARMS="$2"; shift 2 ;;
        --models)      MODELS="$2"; shift 2 ;;
        --max-targets) MAX_TARGETS="$2"; shift 2 ;;
        --calibrate)   CALIBRATE=1; shift ;;
        --yes|-y)      ASSUME_YES=1; shift ;;
        -h|--help)     sed -n '2,74p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------------
# 1. Measure the machine.
# -----------------------------------------------------------------------------
detect_threads() {
    nproc 2>/dev/null && return
    getconf _NPROCESSORS_ONLN 2>/dev/null && return
    sysctl -n hw.ncpu 2>/dev/null && return
    echo 2
}
detect_ram_gb() {
    # MemTotal, not MemAvailable. Under WSL2 the kernel is handed a fraction of
    # the host's RAM, so this reports what the analysis can actually use rather
    # than what the laptop has on the motherboard, and the two differ here.
    if [[ -r /proc/meminfo ]]; then
        awk '/^MemTotal/ {printf "%d", $2 / 1048576}' /proc/meminfo; return
    fi
    if sysctl -n hw.memsize >/dev/null 2>&1; then
        echo $(( $(sysctl -n hw.memsize) / 1073741824 )); return
    fi
    echo 4
}
detect_free_gb() { df -Pk "$1" 2>/dev/null | awk 'NR==2 {printf "%d", $4 / 1048576}' || echo 0; }
detect_free_mb() { df -Pm "$1" 2>/dev/null | awk 'NR==2 {printf "%d", $4}' || echo 0; }
is_wsl() { grep -qi microsoft /proc/version 2>/dev/null; }

# GPU. Every statement this repository makes about hardware is a measurement,
# so the absence of a usable accelerator is established here rather than
# assumed. nvidia-smi absent is not proof on its own: stage 1 found a binary
# linked against CUDA that ran CPU-only on this machine. The prediction stages
# record what the frameworks themselves report at run time.
detect_gpu() {
    local found="none"
    if command -v nvidia-smi >/dev/null 2>&1; then
        local n
        n="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || true)"
        [[ -n "$n" ]] && found="$n"
    fi
    if [[ "$found" == "none" && -d /proc/driver/nvidia ]]; then found="nvidia driver present, nvidia-smi absent"; fi
    echo "$found"
}

DET_THREADS="$(detect_threads)"
DET_RAM_GB="$(detect_ram_gb)"
DET_GPU="$(detect_gpu)"

[[ -n "$DATA_DIR"  ]] || DATA_DIR="${REPO_DIR}/data"
[[ -n "$CACHE_DIR" ]] || CACHE_DIR="${DATA_DIR}/cache"
mkdir -p "$DATA_DIR" "$CACHE_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"
CACHE_DIR="$(cd "$CACHE_DIR" && pwd)"

DET_FREE_DATA="$(detect_free_gb "$DATA_DIR")"
DET_FREE_CACHE="$(detect_free_gb "$CACHE_DIR")"
HOST_MOUNT="/mnt/c"
HOST_FREE_GB=""
HOST_FREE_MB=""
# Held back for the WSL swap file, which lives on the same Windows drive and
# grows under memory pressure. A prediction that pushes the machine into swap
# while the drive is full takes the whole virtual machine down with it.
HOST_RESERVE_MB=3072
if is_wsl && [[ -d "$HOST_MOUNT" ]]; then
    HOST_FREE_GB="$(detect_free_gb "$HOST_MOUNT")"
    HOST_FREE_MB="$(detect_free_mb "$HOST_MOUNT")"
fi

echo "=============================================================="
echo " structure prediction confidence benchmark - configuration"
echo "=============================================================="
echo "  CPU threads visible      : ${DET_THREADS}"
echo "  RAM visible to kernel    : ${DET_RAM_GB} GB"
echo "  GPU                      : ${DET_GPU}"
echo "  free disk at data dir    : ${DET_FREE_DATA} GB  (${DATA_DIR})"
echo "  free disk at cache dir   : ${DET_FREE_CACHE} GB  (${CACHE_DIR})"
if [[ -n "$HOST_FREE_MB" ]]; then
    echo "  free disk on host drive  : ${HOST_FREE_MB} MB  (${HOST_MOUNT})"
    echo "  [note] running under WSL. df inside the virtual machine does not see"
    echo "         the Windows drive that holds the Linux filesystem. A byte"
    echo "         written inside the virtual disk costs a byte of host disk"
    echo "         permanently; a byte written under ${HOST_MOUNT} is released"
    echo "         when the file is deleted. Put --cache-dir under ${HOST_MOUNT}."
fi
echo

ask() {  # ask <prompt> <default> <varname>
    local prompt="$1" default="$2" __var="$3" reply
    if [[ $ASSUME_YES -eq 1 || ! -t 0 ]]; then printf -v "$__var" '%s' "$default"; return; fi
    read -r -p "  ${prompt} [${default}]: " reply || reply=""
    printf -v "$__var" '%s' "${reply:-$default}"
}

# One gigabyte held back for the operating system, the shell and the python
# parent. A prediction that swaps is measuring the swap device.
DEFAULT_RAM=$(( DET_RAM_GB - 1 )); (( DEFAULT_RAM < 1 )) && DEFAULT_RAM=1
DEFAULT_DISK="${HOST_FREE_GB:-$DET_FREE_CACHE}"
DEFAULT_DISK=$(( DEFAULT_DISK - 4 )); (( DEFAULT_DISK < 1 )) && DEFAULT_DISK=1
[[ -n "$THREADS" ]] || ask "CPU threads to use"            "$DET_THREADS" THREADS
[[ -n "$RAM_GB"  ]] || ask "RAM to use (GB)"               "$DEFAULT_RAM" RAM_GB
[[ -n "$DISK_GB" ]] || ask "disk budget, data plus cache (GB)" "$DEFAULT_DISK" DISK_GB
[[ -n "$JOBS"    ]] || ask "concurrent prediction jobs"    "1"            JOBS
[[ -n "$HOURS"   ]] || ask "wall clock budget (hours)"     "72"           HOURS
[[ -n "$MAX_TARGETS" ]] || MAX_TARGETS=0

for pair in "THREADS:$THREADS" "RAM_GB:$RAM_GB" "DISK_GB:$DISK_GB" "JOBS:$JOBS" \
            "HOURS:$HOURS" "MAX_TARGETS:$MAX_TARGETS"; do
    name="${pair%%:*}"; val="${pair#*:}"
    [[ "$val" =~ ^[0-9]+$ ]] || { echo "[error] ${name} must be a non-negative integer, got '${val}'" >&2; exit 1; }
done
(( THREADS >= 1 )) || { echo "[error] need at least 1 thread" >&2; exit 1; }
(( JOBS >= 1 )) || JOBS=1
if (( JOBS > 1 )); then
    echo "[warn] --jobs ${JOBS}. On CPU the inference already spreads over every"
    echo "       core, so concurrent jobs divide the same cores and multiply the"
    echo "       memory. The per-target timings in the README come from the"
    echo "       serial subset and are marked with their job count."
fi

# -----------------------------------------------------------------------------
# 2. The sequence-length ceiling.
#
# Provisional coefficients, used only until --calibrate replaces them. They are
# not measurements and the script says so in its output and in project.conf.
# The form is peak_mb = MEM_BASE_MB + MEM_QUAD_MB_PER_KRES2 * (L/1000)^2, which
# is the pair representation term; the linear term is folded into the base.
# -----------------------------------------------------------------------------
MEM_BASE_MB=${MEM_BASE_MB:-2600}
MEM_QUAD_MB_PER_KRES2=${MEM_QUAD_MB_PER_KRES2:-34000}
SEC_BASE=${SEC_BASE:-120}
SEC_QUAD_PER_KRES2=${SEC_QUAD_PER_KRES2:-9000}
CURVE_SOURCE="provisional, not measured on this machine"

project_mem_mb() {  # project_mem_mb <length>
    awk -v L="$1" -v b="$MEM_BASE_MB" -v q="$MEM_QUAD_MB_PER_KRES2" \
        'BEGIN { printf "%d", b + q * (L/1000.0)^2 }'
}
project_sec() {     # project_sec <length>
    awk -v L="$1" -v b="$SEC_BASE" -v q="$SEC_QUAD_PER_KRES2" \
        'BEGIN { printf "%d", b + q * (L/1000.0)^2 }'
}
# The longest sequence whose projected peak fits the memory budget.
max_len_for_budget() {
    local budget_mb=$(( RAM_GB * 1024 ))
    awk -v b="$MEM_BASE_MB" -v q="$MEM_QUAD_MB_PER_KRES2" -v m="$budget_mb" \
        'BEGIN { if (m <= b || q <= 0) { print 0; exit } printf "%d", 1000 * sqrt((m - b) / q) }'
}

# -----------------------------------------------------------------------------
# 3. --calibrate: measure the curve rather than assume it.
# -----------------------------------------------------------------------------
if (( CALIBRATE == 1 )); then
    [[ -r "$CONF" ]] || { echo "[error] project.conf not found; run without --calibrate first" >&2; exit 3; }
    # The answers already in project.conf are the defaults for this pass, and
    # anything given on the command line wins. Sourcing overwrites every
    # variable, so the command-line values are put aside first and restored
    # afterwards; without that, a calibration pass would silently reset a
    # budget the operator had just changed.
    CLI_THREADS="$THREADS"; CLI_RAM="$RAM_GB"; CLI_DISK="$DISK_GB"; CLI_JOBS="$JOBS"
    CLI_HOURS="$HOURS"; CLI_ARMS="$ARMS"; CLI_MODELS="$MODELS"; CLI_MAX_TARGETS="$MAX_TARGETS"
    CLI_DATA="$DATA_DIR"; CLI_CACHE="$CACHE_DIR"
    # shellcheck source=/dev/null
    source "$CONF"
    HOURS="${TIME_BUDGET_HOURS:-72}"
    [[ -n "$CLI_THREADS"     ]] && THREADS="$CLI_THREADS"
    [[ -n "$CLI_RAM"         ]] && RAM_GB="$CLI_RAM"
    [[ -n "$CLI_DISK"        ]] && DISK_GB="$CLI_DISK"
    [[ -n "$CLI_JOBS"        ]] && JOBS="$CLI_JOBS"
    [[ -n "$CLI_HOURS"       ]] && HOURS="$CLI_HOURS"
    [[ -n "$CLI_MAX_TARGETS" ]] && MAX_TARGETS="$CLI_MAX_TARGETS"
    [[ -n "$CLI_DATA"        ]] && DATA_DIR="$CLI_DATA"
    [[ -n "$CLI_CACHE"       ]] && CACHE_DIR="$CLI_CACHE"
    # ARMS and MODELS have non-empty defaults at the top of this script, so
    # they are compared against those defaults rather than against empty.
    [[ "$CLI_ARMS"   != "$DEFAULT_ARMS"   ]] && ARMS="$CLI_ARMS"
    [[ "$CLI_MODELS" != "$DEFAULT_MODELS" ]] && MODELS="$CLI_MODELS"
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    # shellcheck source-path=SCRIPTDIR
    # shellcheck source=lib_common.sh
    source "${SCRIPT_DIR}/lib_common.sh"
    echo "[00_configure] calibrating on this machine. Three short targets, one"
    echo "               model, no MSA, so the measurement costs minutes rather"
    echo "               than hours and isolates the length term."
    csc_conda_init
    CAL_OUT="${RESULTS_DIR}/environment/memory_curve.tsv"
    mkdir -p "$(dirname "$CAL_OUT")"
    if ! bash "${SCRIPT_DIR}/05_predict.sh" --calibrate --out "$CAL_OUT"; then
        echo "[error] calibration run failed; project.conf keeps its provisional curve" >&2
        exit 5
    fi
    # Least squares on peak_rss_mb = a + b*(L/1000)^2, and the same for
    # seconds. The largest absolute residual is carried out with the
    # coefficients: with three points a quadratic fit always looks good, and a
    # residual of a few hundred megabytes on a measurement of a few thousand is
    # the difference between a ceiling worth trusting and one that is
    # arithmetic dressed up as evidence.
    read -r MEM_BASE_MB MEM_QUAD_MB_PER_KRES2 SEC_BASE SEC_QUAD_PER_KRES2 N_CAL MAX_RESID_MB < <(
        awk -F'\t' 'NR>1 && $3 ~ /^[0-9.]+$/ {
                x[n] = ($2/1000.0)^2; y[n] = $3; t[n] = $4
                sx += x[n]; sy += y[n]; sxx += x[n]*x[n]; sxy += x[n]*y[n]
                st += t[n]; sxt += x[n]*t[n]; n++
            }
            END {
                if (n < 2) { print "NA NA NA NA", n, "NA"; exit }
                d = n*sxx - sx*sx
                if (d == 0) { print "NA NA NA NA", n, "NA"; exit }
                b = (n*sxy - sx*sy) / d; a = (sy - b*sx) / n
                bt = (n*sxt - sx*st) / d; at = (st - bt*sx) / n
                worst = 0
                for (i = 0; i < n; i++) {
                    r = y[i] - (a + b*x[i]); if (r < 0) r = -r
                    if (r > worst) worst = r
                }
                printf "%d %d %d %d %d %d\n", a, b, at, bt, n, worst
            }' "$CAL_OUT"
    )
    [[ "$MEM_BASE_MB" == "NA" ]] && { echo "[error] calibration produced fewer than two usable points" >&2; exit 5; }
    CURVE_SOURCE="measured on this machine, ${N_CAL} points, largest residual ${MAX_RESID_MB} MB, ${CAL_OUT}"
    echo "[00_configure] fitted peak_rss_mb = ${MEM_BASE_MB} + ${MEM_QUAD_MB_PER_KRES2} * (L/1000)^2"
    echo "[00_configure] fitted seconds     = ${SEC_BASE} + ${SEC_QUAD_PER_KRES2} * (L/1000)^2"
    echo "[00_configure] largest memory residual over the ${N_CAL} measured points: ${MAX_RESID_MB} MB"
fi

MAX_SEQ_LEN="$(max_len_for_budget)"
if (( MAX_SEQ_LEN < 60 )); then
    echo "[error] the memory curve says nothing longer than ${MAX_SEQ_LEN} residues" >&2
    echo "        fits in ${RAM_GB} GB. Raise --ram or run this elsewhere." >&2
    exit 4
fi
# The set is built at 80 to 300 residues, so a ceiling above 300 costs nothing
# and a ceiling below it removes the long targets. Both cases are reported.
SET_MAX_LEN=300
EFFECTIVE_MAX=$(( MAX_SEQ_LEN < SET_MAX_LEN ? MAX_SEQ_LEN : SET_MAX_LEN ))

echo "  memory curve             : ${CURVE_SOURCE}"
echo "  peak RSS at 100 / 200 / 300 residues (projected): $(project_mem_mb 100) / $(project_mem_mb 200) / $(project_mem_mb 300) MB"
echo "  maximum sequence length  : ${MAX_SEQ_LEN} residues at ${RAM_GB} GB"
if (( MAX_SEQ_LEN < SET_MAX_LEN )); then
    echo "  [note] below the ${SET_MAX_LEN}-residue upper bound of the held-out set, so"
    echo "         targets between ${MAX_SEQ_LEN} and ${SET_MAX_LEN} residues will be refused"
    echo "         and listed in results/excluded.tsv with their projected peak."
fi
echo

# -----------------------------------------------------------------------------
# 4. Project disk and wall clock for the arms requested.
# -----------------------------------------------------------------------------
n_items() { tr ',' '\n' <<<"$1" | grep -c . ; }
N_ARMS="$(n_items "$ARMS")"
N_MODELS="$(n_items "$MODELS")"
N_TARGETS=${MAX_TARGETS:-0}
(( N_TARGETS == 0 )) && N_TARGETS=100      # planning figure until 02 has run

# Weights. One set is on disk at a time: 05_predict.sh fetches, runs, records
# the checksum and size, then deletes before the next model is fetched. So the
# figure that matters is the largest single set, not their sum.
#
# The sizes below are the bytes each publisher currently serves, rounded up to
# whole megabytes, taken from HTTP headers rather than from documentation.
# 01_install.sh records what actually landed, and if the two disagree the
# recorded figure is the one the README quotes.
#
#   alphafold2_ptm  one parameter tar of 3,722,752,000 bytes. ColabFold streams
#                   it rather than saving it, and extracts ten parameter files,
#                   of which the five this arm uses are about half. The peak is
#                   therefore the extracted set rather than tar plus contents.
#   boltz2          a structure checkpoint of 2,286,561,469 bytes, an affinity
#                   checkpoint of 2,062,139,170 that is fetched whether or not
#                   affinity is asked for, a molecule archive of 1,855,662,080,
#                   and that archive's extracted contents, which the code keeps
#                   alongside it. No affinity is predicted here, so the
#                   affinity checkpoint is pure cost; 01_install.sh tests
#                   whether it can be skipped and records the answer.
MB_WEIGHTS_AF2=3551
MB_WEIGHTS_BOLTZ=8100
MB_WEIGHTS_PEAK=0
case ",${MODELS}," in *,alphafold2_ptm,*) (( MB_WEIGHTS_AF2   > MB_WEIGHTS_PEAK )) && MB_WEIGHTS_PEAK=$MB_WEIGHTS_AF2 ;; esac
case ",${MODELS}," in *,boltz2,*)         (( MB_WEIGHTS_BOLTZ > MB_WEIGHTS_PEAK )) && MB_WEIGHTS_PEAK=$MB_WEIGHTS_BOLTZ ;; esac
# Per target, in the data directory rather than the cache: the alignment the
# server returns, one prediction per arm, the confidence arrays, the reference
# structure, and the scores. Alignments dominate and are capped and compressed.
MB_PER_TARGET=9
# The stage 1 clone, its prepared receptors for the docking subset, and the
# reference structures those need.
MB_STAGE1=400
PROJ_CACHE_MB=$(( MB_WEIGHTS_PEAK ))
PROJ_DATA_MB=$(( N_TARGETS * MB_PER_TARGET + MB_STAGE1 ))
PROJ_DISK_MB=$(( PROJ_CACHE_MB + PROJ_DATA_MB ))
PROJ_DISK_GB=$(awk -v m="$PROJ_DISK_MB" 'BEGIN { printf "%.1f", m/1024 }')

# Wall clock. One prediction per target per arm per model, at the median length
# of the set, plus the five-model subset. This is the number that decides
# whether the design fits, and it is the one most likely to be wrong, so the
# assumption is printed next to it.
MEDIAN_LEN=190
SEC_ONE="$(project_sec "$MEDIAN_LEN")"
PROJ_SEC=$(( N_TARGETS * SEC_ONE * N_ARMS ))
PROJ_HOURS=$(( PROJ_SEC / 3600 ))

echo "  arms requested           : ${N_ARMS} (${ARMS})"
echo "  models requested         : ${N_MODELS} (${MODELS})"
echo "  targets assumed          : ${N_TARGETS}"
echo "  projected disk           : ${PROJ_DISK_GB} GB = ${PROJ_CACHE_MB} MB of weights (one set at a time) plus ${PROJ_DATA_MB} MB of data"
echo "  projected wall clock     : ${PROJ_HOURS} h at ${SEC_ONE} s per prediction of a ${MEDIAN_LEN}-residue target"
echo

# The comparisons are in megabytes. Rounding each side to whole gigabytes first
# loses up to a gigabyte on a budget of ten, which on this machine is the
# difference between a run that fits and a refusal.
BUDGET_MB=$(( DISK_GB * 1024 ))
if (( PROJ_DISK_MB > BUDGET_MB )); then
    echo "[error] the requested arms project to ${PROJ_DISK_MB} MB and the budget is ${BUDGET_MB} MB." >&2
    printf '        Shortfall %.1f GB. Raise --disk, point --cache-dir at a larger\n' \
        "$(awk -v a="$(( PROJ_DISK_MB - BUDGET_MB ))" 'BEGIN{print a/1024}')" >&2
    echo "        filesystem, or drop a model. Refusing to write a project.conf that" >&2
    echo "        would die two thirds of the way through." >&2
    exit 4
fi
if [[ -n "$HOST_FREE_MB" ]] && (( PROJ_DISK_MB + HOST_RESERVE_MB > HOST_FREE_MB )); then
    echo "[error] projection is ${PROJ_DISK_MB} MB plus ${HOST_RESERVE_MB} MB held back for the WSL" >&2
    echo "        swap file, and the Windows drive has ${HOST_FREE_MB} MB free." >&2
    printf '        Shortfall %.1f GB. df inside WSL cannot see this; see data/README.md.\n' \
        "$(awk -v a="$(( PROJ_DISK_MB + HOST_RESERVE_MB - HOST_FREE_MB ))" 'BEGIN{print a/1024}')" >&2
    echo "        Dropping boltz2 from --models leaves the AlphaFold2 arms, whose" >&2
    echo "        weights are $(( MB_WEIGHTS_BOLTZ - MB_WEIGHTS_AF2 )) MB smaller." >&2
    exit 4
fi
if (( PROJ_HOURS > HOURS )); then
    if (( ASSUME_YES == 1 )); then
        echo "[warn] projected ${PROJ_HOURS} h exceeds the ${HOURS} h budget; proceeding because --yes was given."
        echo "       02_build_holdout_set.py will cut the set size to fit and record what it cut."
    else
        echo "[error] projected ${PROJ_HOURS} h exceeds the ${HOURS} h budget." >&2
        echo "        Pass --yes to proceed and let 02_build_holdout_set.py cut the set" >&2
        echo "        size, raise --hours, or drop an arm. Cutting the set is the" >&2
        echo "        preferred response: the arms are the experiment." >&2
        exit 4
    fi
fi
TARGET_BUDGET=$(( HOURS * 3600 / (SEC_ONE * N_ARMS) ))
(( TARGET_BUDGET < 1 )) && TARGET_BUDGET=1

# -----------------------------------------------------------------------------
# 5. Locate conda.
# -----------------------------------------------------------------------------
find_conda_sh() {
    local c root
    for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null || true)"; do
        [[ -n "$c" && -x "$c" ]] || continue
        root="$(dirname "$(dirname "$c")")"
        [[ -r "${root}/etc/profile.d/conda.sh" ]] && { echo "${root}/etc/profile.d/conda.sh"; return; }
    done
    for c in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/anaconda3" /opt/conda; do
        [[ -r "${c}/etc/profile.d/conda.sh" ]] && { echo "${c}/etc/profile.d/conda.sh"; return; }
    done
    echo ""
}
CONDA_SH="$(find_conda_sh)"
[[ -n "$CONDA_SH" ]] || echo "[warn] conda not found. 01_install.sh needs it; install miniforge first."

cat > "$CONF_TMP" <<CONF_EOF
# =============================================================================
#  project.conf - written by scripts/00_configure.sh on $(date -Iseconds)
#  Re-run scripts/00_configure.sh to change these. Do not edit by hand: every
#  later stage sources this file and a hand edit will not survive the next run.
# =============================================================================

# ---- hardware, as measured on this machine --------------------------------
THREADS=${THREADS}
RAM_GB=${RAM_GB}
DISK_GB=${DISK_GB}
JOBS=${JOBS}
# What the GPU probe found. Every hardware claim in the README is read from
# here or from logs/, never assumed. "none" means the predictions below are CPU
# numbers.
GPU_DETECTED="${DET_GPU}"
# Under WSL the Linux filesystem lives in a virtual disk on this Windows drive,
# and the disk gate in lib_common.sh measures the drive rather than the virtual
# filesystem. HOST_RESERVE_MB is held back for the WSL swap file, which lives
# on the same drive and grows under memory pressure.
HOST_DRIVE_MOUNT="${HOST_MOUNT}"
HOST_RESERVE_MB=${HOST_RESERVE_MB}
# What the host drive had free when this file was written. The gate compares
# against the drive as it is at run time, not against this; the figure is here
# so a later failure can be read against the conditions the plan assumed.
HOST_FREE_MB_AT_CONFIGURE=${HOST_FREE_MB:-NA}

# ---- the memory ceiling ---------------------------------------------------
# peak_rss_mb = MEM_BASE_MB + MEM_QUAD_MB_PER_KRES2 * (L/1000)^2
# seconds     = SEC_BASE    + SEC_QUAD_PER_KRES2    * (L/1000)^2
# Source of these coefficients: ${CURVE_SOURCE}
# The quadratic term is the pair representation, which is an L by L tensor.
# A fit over three short targets extrapolated to 300 residues is the weakest
# number this configuration produces, and 05_predict.sh records the actual peak
# for every target so the fit can be checked against the run it predicted.
MEM_BASE_MB=${MEM_BASE_MB}
MEM_QUAD_MB_PER_KRES2=${MEM_QUAD_MB_PER_KRES2}
SEC_BASE=${SEC_BASE}
SEC_QUAD_PER_KRES2=${SEC_QUAD_PER_KRES2}
# Longest sequence whose projected peak fits RAM_GB. 05_predict.sh refuses any
# target above this with the projected peak and the shortfall in gigabytes.
MAX_SEQ_LEN=${MAX_SEQ_LEN}
EFFECTIVE_MAX_SEQ_LEN=${EFFECTIVE_MAX}

# ---- what to run ----------------------------------------------------------
# af2_msa_notmpl   ColabFold AlphaFold2, remote MSA, templates off. The primary
#                  arm. Templates off because ColabFold's template search runs
#                  against a current PDB with no identity filter, so a held-out
#                  target can be handed its own structure as a template.
# af2_msa_tmpl     the same with templates on, to measure what that is worth
#                  and how often it is homology modelling rather than
#                  prediction. Every template hit is recorded with its identity
#                  and release date.
# af2_nomsa        the same weights and the same code with --msa-mode
#                  single_sequence. The no-MSA arm is AlphaFold2 rather than a
#                  separate single-sequence model because the alternatives do
#                  not fit in this machine's memory; 01_install.sh records the
#                  arithmetic and the README says the arm was chosen for that
#                  reason.
# boltz2_msa       Boltz-2 on the same alignment.
# null_template    the highest-identity PDB chain released before the model's
#                  cutoff, copied verbatim and scored as if it were a
#                  prediction. The honest comparator.
# null_unrelated   a chain of similar length from an unrelated entry. The
#                  trivial floor.
ARMS="${ARMS}"
MODELS="${MODELS}"
# 0 means take whatever the filters in 02_build_holdout_set.py produce.
MAX_TARGETS=${MAX_TARGETS}
# What the time budget allows at the projected per-prediction cost. Stage 02
# cuts the set to this rather than dropping an arm.
TARGET_BUDGET=${TARGET_BUDGET}
TIME_BUDGET_HOURS=${HOURS}

# ---- reproducibility ------------------------------------------------------
# One seed for everything stochastic. For AlphaFold2 through ColabFold that is
# --random-seed; for Boltz-2 it is --seed, which seeds the diffusion sampling.
# Neither framework promises bitwise determinism on CPU, and the seed
# replicates below are how much that matters here.
SEED=20260925
SEED_REPLICATES="20260925 20260926 20260927 20260928 20260929"
# One target predicted once per seed, in every arm, because run-to-run spread
# is part of the measurement and is almost never reported. Stage 1 found a
# 9.05 Angstrom spread in top-1 RMSD across five seeds on one complex.
SEED_VARIANCE_N_TARGETS=1

# ---- prediction protocol --------------------------------------------------
# Five models is what ColabFold runs by default and what the field reports, so
# it is the protocol here. On CPU it costs five times one model, which is why
# the full set runs the number below and a subset runs all five: the
# difference between them is a result rather than a quietly chosen default.
AF2_NUM_MODELS=1
AF2_NUM_MODELS_SUBSET=5
AF2_SUBSET_N_TARGETS=20
AF2_MODEL_TYPE=alphafold2_ptm
# Recycles left at the model's own default of 3 for alphafold2_ptm. Raising it
# improves accuracy and costs linear time; leaving it is what the published
# numbers used.
AF2_NUM_RECYCLE=3
AF2_RANK=plddt
# Boltz-2 sampling. The defaults its own documentation states, recorded here
# because a results table that does not say which sampling produced it is not
# reproducible.
BOLTZ_RECYCLING_STEPS=3
BOLTZ_SAMPLING_STEPS=200
BOLTZ_DIFFUSION_SAMPLES=1
BOLTZ_STEP_SCALE=1.5

# ---- held-out set construction --------------------------------------------
# The common cutoff is the latest of the models' training cutoffs, so one set
# is out of training for every model compared. 02_build_holdout_set.py reads
# the per-model cutoffs from config/model_cutoffs.tsv, which carries the source
# each was read from, and derives this date rather than trusting a constant.
SET_MIN_RESOLUTION=2.0
SET_MIN_LENGTH=80
SET_MAX_LENGTH=${SET_MAX_LEN}
SET_MAX_PROTEIN_ENTITIES=1
# Cluster identity for the leakage proxy. A target is dropped if any member of
# its cluster was released before the cutoff. This is a proxy: the training
# sets are not published, and the README says so plainly.
SET_CLUSTER_IDENTITY=30
# Ligand filters for the docking subset, following the published benchmark this
# repository hands its receptors to.
LIG_MIN_WEIGHT=100
LIG_MAX_WEIGHT=900
LIG_MIN_HEAVY_ATOMS=10
LIG_MIN_RSCC=0.95
LIG_MAX_RSR=0.2

# ---- scoring --------------------------------------------------------------
# lDDT is computed CA-only, because that is the quantity the confidence head
# was trained to predict. All-atom lDDT is computed as well and reported in its
# own column; the two are not interchangeable and the calibration plots use the
# CA-only one.
LDDT_INCLUSION_RADIUS=15.0
# Binding-site radius for the pocket-local measures. 5 Angstrom around the
# crystallographic ligand, which is the radius the closest published comparison
# used. It is arbitrary in the sense that no pocket has a 5 Angstrom edge, and
# 06_score_structures.py computes 6 and 8 Angstrom as well so the conclusion
# can be checked against the choice.
POCKET_RADIUS=5.0
POCKET_RADII="4.0 5.0 6.0 8.0"
# The pLDDT band whose reliability the README reports: the fraction of residues
# above this confidence whose actual lDDT-CA falls below the accuracy line.
PLDDT_HIGH=90
LDDT_TRUST=0.7
# Calibration binning, named because expected calibration error depends on it.
CALIBRATION_BINS=10

# ---- paths ----------------------------------------------------------------
REPO_DIR="${REPO_DIR}"
CONFIG_DIR="${REPO_DIR}/config"
SCRIPTS_DIR="${REPO_DIR}/scripts"
DATA_DIR="${DATA_DIR}"
CACHE_DIR="${CACHE_DIR}"
RESULTS_DIR="${REPO_DIR}/results"
FIGURES_DIR="${REPO_DIR}/figures"
LOG_DIR="${REPO_DIR}/logs"

# ---- tools ----------------------------------------------------------------
CONDA_SH="${CONDA_SH}"
# Four environments. ColabFold pins JAX and Boltz pins PyTorch, and neither
# tolerates the other's numpy; OpenStructure comes from bioconda and drags in a
# Qt and OpenMM stack that has no business in an inference environment. Stage 1
# needed four for the same kind of reason and called its tools by absolute path
# rather than activating, which is what happens here.
CONDA_ENV_ANALYSIS=csc_analysis
CONDA_ENV_COLABFOLD=csc_colabfold
CONDA_ENV_BOLTZ=csc_boltz
CONDA_ENV_OST=csc_ost

# ---- the stage 1 pipeline, for the docking handoff ------------------------
# 09_dock_into_predictions.sh clones this repository at this commit and calls
# it unmodified. config/sources.tsv records the same commit with the date it
# was pinned. Two fixes needed by this stage were made there rather than here,
# and this commit is the one that carries them.
STAGE1_REPO="https://github.com/Qasim-Hussain-Code/vina_gnina_pose_benchmark_pipeline.git"
STAGE1_COMMIT=d329fbe0c08e4b36ed89f8ce077febc606350254
STAGE1_DIR="${DATA_DIR}/stage1"
# Vina and Vinardo only. GNINA's release binary is 2.1 GB and needs a CUDA
# runtime environment of about the same size, and stage 1 measured its CNN
# ensemble as unable to finish one complex in eighteen minutes on this machine.
# Dropping it is a disk and time decision, recorded rather than silent.
STAGE1_METHODS=vina,vinardo
STAGE1_ARMS=A0_null,A2_genconf_refbox
CONF_EOF

mv "$CONF_TMP" "$CONF"
echo "  wrote ${CONF}"
echo "  threads=${THREADS} ram=${RAM_GB}GB disk=${DISK_GB}GB jobs=${JOBS} max_len=${MAX_SEQ_LEN}"
echo "  data_dir=${DATA_DIR}"
echo "  cache_dir=${CACHE_DIR}"
if (( CALIBRATE == 0 )); then
    echo
    echo "  The memory curve above is a projection. Run 01_install.sh, then"
    echo "  bash scripts/00_configure.sh --calibrate to replace it with a"
    echo "  measurement from this machine."
fi
