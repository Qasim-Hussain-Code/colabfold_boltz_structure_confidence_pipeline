#!/usr/bin/env bash
# =============================================================================
#  01_install.sh - build the four environments, record exactly what landed
# =============================================================================
#  Nothing here hardcodes a tool version except where a version is the point.
#  Each environment is solved on the day, and whatever the solver and pip pick
#  is written to files the README quotes:
#
#      results/environment/versions.tsv        one row per tool, version, source
#      results/environment/tool_licences.tsv   licence, as the packager records it
#      results/environment/install_cost.tsv    what each environment cost, measured
#      config/env_*.lock.yml                   each solved environment, pinned
#
#  Why four environments.
#
#  ColabFold pins JAX. Boltz pins PyTorch and numpy below 2. OpenStructure
#  comes from bioconda in a build that is not headless on this platform, so it
#  drags in Qt, mesa and OpenMM, and OpenMM drags in CUDA libraries that no
#  part of this benchmark can use. Putting any two of those together means one
#  of them silently gets a version it was not tested against. Stage 1 needed
#  four environments for the same kind of reason and found the problem only
#  because it recorded versions; the same table is written here.
#
#  Three version choices are deliberate and are recorded with their reasons.
#
#  ColabFold comes from PyPI rather than bioconda. The bioconda package is
#  several releases behind, still requires TensorFlow, and pins JAX below what
#  the current release needs. The PyPI version also carries a fix, in the
#  release used here, for a defect where single-sequence mode quietly added one
#  randomly chosen homologue per seed when given a cached alignment. The
#  no-MSA arm of this benchmark is exactly that code path, so an older version
#  would have measured something other than what it claimed.
#
#  Boltz comes from a pinned commit on the development branch, not from the
#  release. The current release runs the diffusion module in reduced precision
#  on CPU, which has been reported to produce distorted geometry and lower
#  confidence than the same input on a GPU. The fix is merged but unreleased,
#  and the version string did not change with it, so the commit is what
#  identifies this build. The fix is partial: the trunk still runs in mixed
#  precision on CPU. 07_validate_geometry.py measures whether the reported
#  failure mode appears here rather than assuming either way.
#
#  The CPU-only PyTorch wheel is installed before Boltz. Boltz requires torch
#  and does not pin it, so installing it first pulls the default build and
#  several gigabytes of CUDA libraries onto a machine with no GPU.
#
#  Usage:
#      bash scripts/01_install.sh                 # skip what is already present
#      bash scripts/01_install.sh --force         # rebuild every environment
#      bash scripts/01_install.sh --only ost      # one environment
#      bash scripts/01_install.sh --record-only   # re-record versions, install nothing
#
#  Options:
#      --force        remove and rebuild the environments
#      --only NAME    analysis, colabfold, boltz or ost
#      --record-only  skip installation and rewrite the version tables
#      -h, --help     this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
csc_load_conf

STAGE=01_install
FORCE=0; ONLY=""; RECORD_ONLY=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)       FORCE=1; shift ;;
        --only)        ONLY="$2"; shift 2 ;;
        --record-only) RECORD_ONLY=1; shift ;;
        -h|--help)     sed -n '2,60p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

(( RECORD_ONLY == 0 )) && { csc_skip_if_done "$STAGE" "$FORCE" && exit 0; }

ENVDIR="${RESULTS_DIR}/environment"
mkdir -p "$ENVDIR"

# Pinned versions. These are the only hardcoded versions in the repository, and
# each is here because the specific version is part of the method rather than
# an implementation detail. config/sources.tsv carries the same values with the
# date each was checked.
COLABFOLD_VERSION="1.6.3"
JAX_VERSION="0.10.2"
BOLTZ_COMMIT="b1ebfc46ecf57f5414e0d1a6f9027bbb122c53bc"
TORCH_CPU_INDEX="https://download.pytorch.org/whl/cpu"

csc_stage_start "$STAGE"
csc_conda_init
# conda-forge must win outright, for the reason stage 1 recorded: with flexible
# priority the solver is free to satisfy a dependency from another channel with
# an older build, and the result is a version nobody chose.
conda config --set channel_priority strict >/dev/null 2>&1 || true

# Appended, not rewritten. An environment is built once and reused, so a later
# run of this stage has nothing to measure; rewriting the header would throw
# away the only record of what the build cost.
COST="${ENVDIR}/install_cost.tsv"
[[ -s "$COST" ]] || printf 'environment\telapsed_s\thost_drive_mb\twsl_used_mb\tnote\trecorded\n' > "$COST"

host_free_mb() { csc_host_free_mb; }
wsl_used_mb()  { df -Pm / 2>/dev/null | awk 'NR==2 {print $3}'; }

# build_env <suffix>
# Creates the environment named in project.conf for that suffix, from the
# matching spec, and measures what it cost. An environment that already exists
# is reused and said so: re-solving it would change versions under a run that
# has already produced results with the old ones.
build_env() {
    local suffix="$1" env spec lock h0 h1 w0 w1 t0 t1
    case "$suffix" in
        analysis)  env="$CONDA_ENV_ANALYSIS" ;;
        colabfold) env="$CONDA_ENV_COLABFOLD" ;;
        boltz)     env="$CONDA_ENV_BOLTZ" ;;
        ost)       env="$CONDA_ENV_OST" ;;
        *) echo "[error] unknown environment ${suffix}" >&2; return 1 ;;
    esac
    [[ -n "$ONLY" && "$ONLY" != "$suffix" ]] && return 0
    spec="${CONFIG_DIR}/env_${suffix}.yml"
    lock="${CONFIG_DIR}/env_${suffix}.lock.yml"

    if (( FORCE == 1 )); then
        echo "[${STAGE}] removing ${env}"
        conda env remove -n "$env" -y >/dev/null 2>&1 || true
    fi
    if conda env list | awk '{print $1}' | grep -qx "$env"; then
        echo "[${STAGE}] ${env} already exists; not re-solving"
        return 0
    fi
    [[ -s "$spec" ]] || { echo "[error] ${spec} missing" >&2; return 1; }

    h0="$(host_free_mb)"; w0="$(wsl_used_mb)"; t0="$(date +%s)"
    if [[ -s "$lock" ]]; then
        echo "[${STAGE}] creating ${env} from the pinned lock file"
        csc_run "env_${suffix}_lock" conda env create -n "$env" -f "$lock" -y
    else
        echo "[${STAGE}] solving ${env} from ${spec}"
        csc_run "env_${suffix}_spec" conda env create -n "$env" -f "$spec" -y
    fi
    t1="$(date +%s)"; h1="$(host_free_mb)"; w1="$(wsl_used_mb)"
    printf '%s\t%d\t%s\t%d\t%s\t%s\n' "$env" "$(( t1 - t0 ))" \
        "$( [[ -n "$h0" && -n "$h1" ]] && echo $(( h0 - h1 )) || echo NA )" \
        "$(( w1 - w0 ))" "conda solve" "$(date -Iseconds)" >> "$COST"
}

if (( RECORD_ONLY == 0 )); then
    build_env analysis
    build_env colabfold
    build_env boltz
    build_env ost
fi

PY_ANALYSIS="$(csc_env_bin "$CONDA_ENV_ANALYSIS" python || true)"
PY_COLABFOLD="$(csc_env_bin "$CONDA_ENV_COLABFOLD" python || true)"
PY_BOLTZ="$(csc_env_bin "$CONDA_ENV_BOLTZ" python || true)"
PY_OST="$(csc_env_bin "$CONDA_ENV_OST" python || true)"
[[ -n "$PY_ANALYSIS" ]] || { echo "[error] ${CONDA_ENV_ANALYSIS} has no python after install" >&2; exit 1; }

# -----------------------------------------------------------------------------
# The pip half. conda supplies the interpreter and pip supplies the model code,
# because neither ColabFold nor Boltz has a usable conda package.
# -----------------------------------------------------------------------------
pip_install_colabfold() {
    [[ -n "$ONLY" && "$ONLY" != "colabfold" ]] && return 0
    [[ -n "$PY_COLABFOLD" ]] || { echo "[warn] ${CONDA_ENV_COLABFOLD} missing; skipping"; return 0; }
    if "$PY_COLABFOLD" -c 'import colabfold' >/dev/null 2>&1 && (( FORCE == 0 )); then
        echo "[${STAGE}] colabfold already installed in ${CONDA_ENV_COLABFOLD}"
        return 0
    fi
    # The prebuilt hhsearch and kalign binaries the templates arm needs are
    # compiled for AVX2. Without it the templates arm cannot run, and that is
    # worth knowing now rather than at stage 5.
    if ! grep -qm1 avx2 /proc/cpuinfo; then
        echo "[warn] this CPU does not report AVX2. The prebuilt template search"
        echo "       binaries will not run, so the templates-on arm will fail."
    fi
    local h0 t0 h1 t1; h0="$(host_free_mb)"; t0="$(date +%s)"
    echo "[${STAGE}] pip installing colabfold ${COLABFOLD_VERSION} with jax ${JAX_VERSION}"
    csc_run "pip_colabfold" "$PY_COLABFOLD" -m pip install --no-input -q \
        "colabfold[alphafold]==${COLABFOLD_VERSION}" \
        "jax==${JAX_VERSION}" "jaxlib==${JAX_VERSION}"
    t1="$(date +%s)"; h1="$(host_free_mb)"
    printf '%s\t%d\t%s\t%s\t%s\t%s\n' "${CONDA_ENV_COLABFOLD}/pip" "$(( t1 - t0 ))" \
        "$( [[ -n "$h0" && -n "$h1" ]] && echo $(( h0 - h1 )) || echo NA )" NA \
        "colabfold ${COLABFOLD_VERSION}" "$(date -Iseconds)" >> "$COST"
}

pip_install_boltz() {
    [[ -n "$ONLY" && "$ONLY" != "boltz" ]] && return 0
    [[ -n "$PY_BOLTZ" ]] || { echo "[warn] ${CONDA_ENV_BOLTZ} missing; skipping"; return 0; }
    if "$PY_BOLTZ" -c 'import boltz' >/dev/null 2>&1 && (( FORCE == 0 )); then
        echo "[${STAGE}] boltz already installed in ${CONDA_ENV_BOLTZ}"
        return 0
    fi
    local h0 t0 h1 t1; h0="$(host_free_mb)"; t0="$(date +%s)"
    echo "[${STAGE}] pip installing the CPU-only torch wheel first"
    csc_run "pip_torch_cpu" "$PY_BOLTZ" -m pip install --no-input -q \
        torch --index-url "$TORCH_CPU_INDEX"
    echo "[${STAGE}] pip installing boltz at ${BOLTZ_COMMIT:0:12}"
    csc_run "pip_boltz" "$PY_BOLTZ" -m pip install --no-input -q \
        "boltz @ git+https://github.com/jwohlwend/boltz@${BOLTZ_COMMIT}"
    t1="$(date +%s)"; h1="$(host_free_mb)"
    printf '%s\t%d\t%s\t%s\t%s\t%s\n' "${CONDA_ENV_BOLTZ}/pip" "$(( t1 - t0 ))" \
        "$( [[ -n "$h0" && -n "$h1" ]] && echo $(( h0 - h1 )) || echo NA )" NA \
        "torch cpu wheel plus boltz at ${BOLTZ_COMMIT:0:12}" "$(date -Iseconds)" >> "$COST"
}

if (( RECORD_ONLY == 0 )); then
    pip_install_colabfold
    pip_install_boltz
fi

# -----------------------------------------------------------------------------
# Prove each tool runs before anything depends on it. Stage 1 learned this the
# expensive way: a binary described as static turned out to need five CUDA
# libraries and exited before printing its version.
# -----------------------------------------------------------------------------
echo "[${STAGE}] checking that each tool starts"
CHECKS="${ENVDIR}/tool_checks.tsv"
printf 'tool\tenvironment\tstatus\tdetail\trecorded\n' > "$CHECKS"
# The environment name rather than the absolute path. The path is particular to
# one machine and one account, and every tracked file in this repository is
# meant to be readable by a stranger without telling them where it was built.
check() {  # check <name> <env> <detail> <command...>
    local name="$1" env="$2" detail="$3"; shift 3
    if "$@" >/dev/null 2>&1; then
        printf '%s\t%s\tok\t%s\t%s\n' "$name" "$env" "$detail" "$(date -Iseconds)" >> "$CHECKS"
        echo "  ${name}: ok"
    else
        printf '%s\t%s\tfailed\t%s\t%s\n' "$name" "$env" "$detail" "$(date -Iseconds)" >> "$CHECKS"
        echo "  ${name}: FAILED"
        return 1
    fi
}

CF_BIN="$(csc_env_bin "$CONDA_ENV_COLABFOLD" colabfold_batch || true)"
BOLTZ_BIN="$(csc_env_bin "$CONDA_ENV_BOLTZ" boltz || true)"
OST_BIN="$(csc_env_bin "$CONDA_ENV_OST" ost || true)"
USALIGN_BIN="$(csc_env_bin "$CONDA_ENV_OST" USalign || true)"
TMSCORE_BIN="$(csc_env_bin "$CONDA_ENV_OST" TMscore || true)"

RC=0
[[ -n "$CF_BIN"      ]] && { check colabfold_batch "$CONDA_ENV_COLABFOLD" "AlphaFold2 inference" "$CF_BIN" --help || RC=1; }
[[ -n "$BOLTZ_BIN"   ]] && { check boltz "$CONDA_ENV_BOLTZ" "Boltz-2 inference" "$BOLTZ_BIN" --help || RC=1; }
[[ -n "$OST_BIN"     ]] && { check ost "$CONDA_ENV_OST" "lDDT, TM-score, stereochemistry" "$OST_BIN" compare-structures --help || RC=1; }
[[ -n "$USALIGN_BIN" ]] && { check USalign "$CONDA_ENV_OST" "TM-score" "$USALIGN_BIN" -h || RC=1; }
[[ -n "$TMSCORE_BIN" ]] && { check TMscore "$CONDA_ENV_OST" "GDT-TS and GDT-HA" "$TMSCORE_BIN" -h || RC=1; }

# The accelerator each framework actually found. Every hardware claim in the
# README is a measurement, and this is where the CPU claim comes from.
ACCEL="${ENVDIR}/accelerator.tsv"
printf 'framework\tbackend\tdetail\trecorded\n' > "$ACCEL"
if [[ -n "$PY_COLABFOLD" ]]; then
    "$PY_COLABFOLD" - "$ACCEL" <<'PY' || true
import datetime, sys
try:
    import jax
    row = ("jax", jax.default_backend(), f"devices={jax.devices()}")
except Exception as exc:
    row = ("jax", "error", repr(exc)[:120])
with open(sys.argv[1], "a") as fh:
    fh.write("\t".join(row) + "\t" + datetime.datetime.now().astimezone().isoformat(timespec="seconds") + "\n")
PY
fi
if [[ -n "$PY_BOLTZ" ]]; then
    "$PY_BOLTZ" - "$ACCEL" <<'PY' || true
import datetime, sys
try:
    import torch
    backend = "cuda" if torch.cuda.is_available() else "cpu"
    row = ("torch", backend,
           f"version={torch.__version__} cuda_build={torch.version.cuda} threads={torch.get_num_threads()}")
except Exception as exc:
    row = ("torch", "error", repr(exc)[:120])
with open(sys.argv[1], "a") as fh:
    fh.write("\t".join(row) + "\t" + datetime.datetime.now().astimezone().isoformat(timespec="seconds") + "\n")
PY
fi

# -----------------------------------------------------------------------------
# Record what landed.
# -----------------------------------------------------------------------------
echo "[${STAGE}] recording versions"
VER="${ENVDIR}/versions.tsv"
printf 'tool\tversion\tsource\tdetail\trecorded\n' > "$VER"
rec() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$(date -Iseconds)" >> "$VER"; }

# Versions come from the package manager rather than from each tool's own
# --version, because several of these print a commit or nothing useful. Boltz
# is the clearest case: the development build reports the released version
# number, so the commit is recorded next to it.
record_conda_pkgs() {  # record_conda_pkgs <env> <pkg...>
    local env="$1"; shift
    local listing; listing="$(conda list -n "$env" 2>/dev/null)" || return 0
    local p v b
    for p in "$@"; do
        v="$(awk -v p="$p" '$1==p {print $2; exit}' <<<"$listing")"
        b="$(awk -v p="$p" '$1==p {print $3; exit}' <<<"$listing")"
        [[ -n "$v" ]] && rec "$p" "$v" "conda/${env}" "$b"
    done
    return 0
}
record_pip_pkgs() {  # record_pip_pkgs <python> <env label> <pkg...>
    local py="$1" label="$2"; shift 2
    [[ -n "$py" ]] || return 0
    local listing; listing="$("$py" -m pip list 2>/dev/null)" || return 0
    local p v
    for p in "$@"; do
        v="$(awk -v p="$p" 'tolower($1)==tolower(p) {print $2; exit}' <<<"$listing")"
        [[ -n "$v" ]] && rec "$p" "$v" "pip/${label}" ""
    done
    return 0
}

record_conda_pkgs "$CONDA_ENV_ANALYSIS" python gemmi biopython numpy scipy pandas \
                  matplotlib-base seaborn statsmodels requests quarto shellcheck
record_conda_pkgs "$CONDA_ENV_OST" python openstructure usalign numpy
record_pip_pkgs "$PY_COLABFOLD" "$CONDA_ENV_COLABFOLD" colabfold alphafold-colabfold \
                jax jaxlib numpy biopython dm-haiku colabfold-binaries \
                colabfold-legacy-kernels absl-py ml-collections
record_pip_pkgs "$PY_BOLTZ" "$CONDA_ENV_BOLTZ" boltz torch numpy pytorch-lightning \
                rdkit fairscale torchmetrics

# The Boltz commit, which the version string does not carry.
[[ -n "$PY_BOLTZ" ]] && rec boltz_commit "$BOLTZ_COMMIT" "git" \
    "development branch; the release of this version runs the diffusion module in reduced precision on CPU"
# Whether the precision fix is in the installed code, and whether the trunk
# precision it does not change is still mixed. Both are read from the installed
# files rather than assumed from the commit.
if [[ -n "$PY_BOLTZ" ]]; then
    BOLTZ_PKG="$("$PY_BOLTZ" -c 'import boltz, pathlib; print(pathlib.Path(boltz.__file__).parent)' 2>/dev/null || true)"
    if [[ -n "$BOLTZ_PKG" && -r "${BOLTZ_PKG}/model/models/boltz2.py" ]]; then
        if grep -q "autocast_device_type" "${BOLTZ_PKG}/model/models/boltz2.py"; then
            rec boltz_cpu_precision_fix present "source" "autocast disabled for the active device in boltz2.py"
        else
            rec boltz_cpu_precision_fix absent "source" "boltz2.py still disables autocast for cuda only"
        fi
    fi
    if [[ -r "${BOLTZ_PKG}/main.py" ]]; then
        TRUNK="$(grep -m1 'precision=' "${BOLTZ_PKG}/main.py" | sed 's/^ *//' | cut -c1-90)"
        rec boltz_trunk_precision "${TRUNK:-unknown}" "source" \
            "the merged fix covers the diffusion and affinity modules, not the trunk"
    fi
fi

# The chemical component dictionary OpenStructure will score against. Its
# snapshot date decides which ligands and modified residues are recognised, and
# it is a property of the package rather than of this repository.
if [[ -n "$PY_OST" ]]; then
    "$PY_OST" - "$VER" <<'PY' || true
import datetime, sys
stamp = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
rows = []
try:
    from ost import conop
    lib = conop.GetDefaultLib()
    if lib is None:
        rows.append(("ost_compound_lib", "absent", "openstructure",
                     "compare-structures will refuse to run without one"))
    else:
        rows.append(("ost_compound_lib", str(lib.GetCreationDate()), "openstructure",
                     f"built with OpenStructure {lib.GetOSTVersionUsed()}"))
except Exception as exc:
    rows.append(("ost_compound_lib", "error", "openstructure", repr(exc)[:110]))
with open(sys.argv[1], "a") as fh:
    for r in rows:
        fh.write("\t".join(r) + "\t" + stamp + "\n")
PY
fi

# The platform itself.
rec conda  "$(conda --version 2>&1 | awk '{print $2}')" "base" "$(command -v conda)"
rec kernel "$(uname -r)" "host" "$(uname -s) $(uname -m)"
rec glibc  "$(ldd --version 2>/dev/null | head -1 | awk '{print $NF}')" "host" "wheel compatibility floor"
rec cpu_flags "$(grep -o -m1 -E 'avx512f|avx2' /proc/cpuinfo | sort -u | tr '\n' ',' | sed 's/,$//')" \
    "host" "the prebuilt template search binaries require avx2"
rec gpu "${GPU_DETECTED:-none}" "00_configure.sh" "every timing in this repository is a CPU timing"

# Licences, read from each environment's own package metadata. This is what the
# run used, which is the claim the LICENSE file has to support.
LIC="${ENVDIR}/tool_licences.tsv"
printf 'package\tversion\tlicence\tenvironment\tchecked\n' > "$LIC"
"$PY_ANALYSIS" - "$LIC" "$(csc_conda_root)" \
    "$CONDA_ENV_ANALYSIS" "$CONDA_ENV_OST" <<'PY' || true
# conda records a licence string in each environment's conda-meta JSON. `conda
# list` does not print it, and `conda search --info` hits the network and can
# describe a different build from the installed one.
import datetime, json, pathlib, sys

out, root = sys.argv[1], sys.argv[2]
envs = sys.argv[3:]
wanted = {"openstructure", "usalign", "python", "gemmi", "biopython", "numpy",
          "scipy", "pandas", "matplotlib-base", "seaborn", "statsmodels",
          "quarto", "shellcheck", "requests"}
stamp = datetime.date.today().isoformat()
rows = 0
with open(out, "a") as fh:
    for env in envs:
        meta = pathlib.Path(root) / "envs" / env / "conda-meta"
        if not meta.is_dir():
            continue
        for p in sorted(meta.glob("*.json")):
            try:
                d = json.loads(p.read_text())
            except Exception:
                continue
            if d.get("name") in wanted:
                fh.write(f"{d['name']}\t{d.get('version','NA')}\t"
                         f"{d.get('license') or 'NA'}\t{env}\t{stamp}\n")
                rows += 1
print(f"[01_install] wrote {rows} conda licence rows")
PY
# pip records licences in wheel metadata, which is a different format.
for pair in "${PY_COLABFOLD}:${CONDA_ENV_COLABFOLD}" "${PY_BOLTZ}:${CONDA_ENV_BOLTZ}"; do
    py="${pair%%:*}"; env="${pair##*:}"
    [[ -n "$py" && -x "$py" ]] || continue
    "$py" - "$LIC" "$env" <<'PY' || true
import datetime, sys
from importlib import metadata

out, env = sys.argv[1], sys.argv[2]
wanted = {"colabfold", "alphafold-colabfold", "jax", "jaxlib", "boltz", "torch",
          "pytorch-lightning", "rdkit", "numpy", "dm-haiku", "colabfold-binaries"}
stamp = datetime.date.today().isoformat()


def field(meta, key):
    # importlib.metadata warns when a missing key is read with __getitem__ and
    # will raise for it in a later version, so every read goes through get().
    return meta.get(key) or ""


with open(out, "a") as fh:
    for dist in metadata.distributions():
        meta = dist.metadata
        name = field(meta, "Name").lower()
        if name not in wanted:
            continue
        lic = field(meta, "License-Expression") or field(meta, "License")
        # A few projects put the whole licence text in that field. Where it is
        # long, the classifier is the short name worth recording.
        if not lic or len(lic) > 60:
            lic = next((c.split("::")[-1].strip()
                        for c in meta.get_all("Classifier") or []
                        if c.startswith("License ::")), lic[:60] or "NA")
        fh.write(f"{name}\t{dist.version}\t{lic}\t{env}\t{stamp}\n")
PY
done

# Solved environments, pinned. --no-builds keeps the version pins and drops the
# build strings, which is what travels between machines; the explicit list next
# to it keeps the build strings for exact reproduction on this platform.
for pair in "${CONDA_ENV_ANALYSIS}:analysis" "${CONDA_ENV_OST}:ost" \
            "${CONDA_ENV_COLABFOLD}:colabfold" "${CONDA_ENV_BOLTZ}:boltz"; do
    env="${pair%%:*}"; base="${pair##*:}"
    conda env list | awk '{print $1}' | grep -qx "$env" || continue
    conda env export -n "$env" --no-builds 2>/dev/null | grep -v '^prefix:' \
        > "${CONFIG_DIR}/env_${base}.lock.yml" || true
    conda list -n "$env" --explicit > "${ENVDIR}/explicit_${base}_$(uname -m).txt" 2>/dev/null || true
done
# pip's half of the two mixed environments, which conda env export does record
# but only as a nested pip block; a plain list is easier to act on.
#
# --format=freeze rather than pip freeze. A plain freeze writes any package
# that conda installed as a name followed by a local file URL pointing at the
# build directory on whichever machine built the conda-forge package. That is
# useless for reproducing anything, and it puts a build machine's directory
# layout into a tracked file. The name and version are what a reader needs,
# so the suffix is cut.
record_pip_freeze() {  # record_pip_freeze <python> <outfile>
    local py="$1" out="$2"
    [[ -n "$py" && -x "$py" ]] || return 0
    "$py" -m pip list --format=freeze 2>/dev/null | sed 's/ @ .*//' > "$out" || true
}
record_pip_freeze "$PY_COLABFOLD" "${ENVDIR}/pip_freeze_colabfold.txt"
record_pip_freeze "$PY_BOLTZ" "${ENVDIR}/pip_freeze_boltz.txt"

echo "[${STAGE}] versions recorded in ${VER}"
column -t -s "$(printf '\t')" "$VER" | sed 's/^/    /'

csc_stage_end "rc=${RC}"
(( RC == 0 )) || { echo "[error] a tool did not start; see ${CHECKS}" >&2; exit 5; }
csc_mark_done "$STAGE"
