#!/usr/bin/env bash
# =============================================================================
#  03_fetch_references.sh - the experimental structures every score is against
# =============================================================================
#  One mmCIF per target, and for the docking subset one ligand per target as
#  SDF. Both come from the archive rather than from a mirror, both are
#  checksummed on arrival, and the size and checksum of every file is recorded
#  so a rerun can tell whether it got the same bytes.
#
#  mmCIF rather than the legacy format. The legacy format cannot represent a
#  chain identifier longer than one character or an entry with more than
#  99,999 atoms, and it carries the author numbering only. The scoring stage
#  needs the label numbering to line a prediction numbered from one up against
#  a deposited model numbered from anything, so the format that carries both is
#  the one fetched. Where a later stage needs the legacy format, it converts
#  and the conversion is visible.
#
#  The ligands come from the model server rather than from the entry, because
#  an SDF cut from an entry by hand loses the bond orders, and a ligand without
#  bond orders is a different molecule to every tool that reads it. Stage 1
#  learned this the expensive way: a pose read back without its bond orders
#  could not be matched against the reference and silently left the denominator.
#
#  Usage:
#      bash scripts/03_fetch_references.sh
#      bash scripts/03_fetch_references.sh --limit 5
#
#  Options:
#      --limit N   stop after N targets
#      --force     refetch what is already on disk
#      -h, --help  this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
csc_load_conf

STAGE=03_fetch_references
LIMIT=0; FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --limit)   LIMIT="$2"; shift 2 ;;
        --force)   FORCE=1; shift ;;
        -h|--help) sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

csc_skip_if_done "$STAGE" "$FORCE" && exit 0

MANIFEST="${CONFIG_DIR}/targets.tsv"
[[ -s "$MANIFEST" ]] || { echo "[error] ${MANIFEST} missing; run 02_build_holdout_set.py" >&2; exit 3; }
REF_DIR="${DATA_DIR}/reference"
LIG_DIR="${DATA_DIR}/reference_ligands"
mkdir -p "$REF_DIR" "$LIG_DIR"

PROV="${RESULTS_DIR}/reference_provenance.tsv"
[[ -s "$PROV" ]] || printf 'target_id\tkind\tfile\tbytes\tsha256\tsource\tstatus\treason\trecorded\n' > "$PROV"

csc_stage_start "$STAGE"

fetch_to() {  # fetch_to <url> <dest>
    curl -sS --fail --max-time 180 --retry 3 --retry-all-errors -o "$2" "$1"
}

record() {  # record <target> <kind> <file> <source> <status> <reason>
    local f="$3" bytes="" sha=""
    if [[ -s "$f" ]]; then
        bytes="$(stat -c %s "$f")"
        sha="$(sha256sum "$f" | awk '{print $1}')"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$(basename "$f")" "$bytes" "$sha" "$4" "$5" "$6" "$(date -Iseconds)" >> "$PROV"
}

mapfile -t ROWS < <(awk -F'\t' '
    NR==1 { for (i = 1; i <= NF; i++) { if ($i == "target_id") id = i; if ($i == "pdb_id") p = i } next }
    id && p { print $id "\t" $p }' "$MANIFEST")
(( LIMIT > 0 )) && ROWS=("${ROWS[@]:0:LIMIT}")

N_OK=0; N_SKIP=0; N_FAIL=0
for row in "${ROWS[@]}"; do
    target="${row%%$'\t'*}"; pdb="${row##*$'\t'}"
    dest="${REF_DIR}/${target}.cif.gz"
    if [[ -s "$dest" && $FORCE -eq 0 ]]; then
        N_SKIP=$(( N_SKIP + 1 )); continue
    fi
    if fetch_to "https://files.rcsb.org/download/${pdb}.cif.gz" "$dest"; then
        record "$target" structure "$dest" "files.rcsb.org" ok ""
        N_OK=$(( N_OK + 1 ))
    else
        rm -f "$dest"
        record "$target" structure "$dest" "files.rcsb.org" failed "download failed"
        csc_record_exclusion "$target" "$STAGE" "" "the reference structure could not be downloaded" ""
        N_FAIL=$(( N_FAIL + 1 ))
    fi
done
echo "[${STAGE}] structures: ${N_OK} fetched, ${N_SKIP} already present, ${N_FAIL} failed"

# --- ligands for the docking subset -----------------------------------------
SUBSET="${CONFIG_DIR}/docking_subset.tsv"
L_OK=0; L_SKIP=0; L_FAIL=0
if [[ -s "$SUBSET" ]]; then
    mapfile -t LROWS < <(awk -F'\t' '
        NR==1 { for (i = 1; i <= NF; i++) { h[$i] = i } next }
        { print $h["target_id"] "\t" $h["pdb_id"] "\t" $h["label_asym_id"] "\t" $h["comp_id"] }' "$SUBSET")
    (( LIMIT > 0 )) && LROWS=("${LROWS[@]:0:LIMIT}")
    for row in "${LROWS[@]}"; do
        IFS=$'\t' read -r target pdb asym comp <<< "$row"
        [[ -n "$asym" ]] || continue
        dest="${LIG_DIR}/${target}__${comp}.sdf"
        if [[ -s "$dest" && $FORCE -eq 0 ]]; then
            L_SKIP=$(( L_SKIP + 1 )); continue
        fi
        # The instance is named explicitly. The endpoint returns the first
        # group that matches whatever it was given, so asking by component
        # alone would return an arbitrary copy when an entry holds several.
        url="https://models.rcsb.org/v1/${pdb}/ligand?label_asym_id=${asym}&encoding=sdf"
        if fetch_to "$url" "$dest" && [[ -s "$dest" ]]; then
            record "$target" ligand "$dest" "models.rcsb.org" ok "instance ${asym} of ${comp}"
            L_OK=$(( L_OK + 1 ))
        else
            rm -f "$dest"
            record "$target" ligand "$dest" "models.rcsb.org" failed "download failed"
            L_FAIL=$(( L_FAIL + 1 ))
        fi
    done
    echo "[${STAGE}] ligands: ${L_OK} fetched, ${L_SKIP} already present, ${L_FAIL} failed"
else
    echo "[${STAGE}] no docking subset yet; ligands not fetched"
fi

csc_stage_end "structures=${N_OK} ligands=${L_OK} failed=$(( N_FAIL + L_FAIL ))"
csc_mark_done "$STAGE"
