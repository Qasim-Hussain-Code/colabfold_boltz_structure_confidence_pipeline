# What belongs in data/

Nothing here is tracked by git except this file. Everything else is either
downloaded from a source recorded in `config/sources.tsv` or written by a script
in `scripts/`, so it is all reproducible from a clean clone.

Two directories, not one, and the difference matters on the machine this was
developed on. `DATA_DIR` holds alignments, predictions, references and scores,
which are small and are needed for the whole run. `CACHE_DIR` holds model
weights, which are gigabytes and are fetched, used and deleted one set at a
time. `scripts/00_configure.sh --data-dir` and `--cache-dir` place them, and
`project.conf` records where they went.

## If you are running under a Windows virtual machine, read this first

The Linux filesystem lives inside a virtual disk file on the Windows drive.
That file grows when anything is written inside it and does not shrink when the
file is deleted. So `df` inside the virtual machine is not measuring what you
think it is measuring.

On this machine, `df /` reported 896 GB free while the Windows drive holding
the virtual disk had 6 GB. A previous stage of this work recorded the same trap
and described it; this stage measures the Windows drive directly, in
`csc_host_free_mb` in `scripts/lib_common.sh`, and every disk gate refuses
against that figure rather than against the one from inside.

Two consequences shape the layout here:

A file written under `/mnt/c` is on the Windows drive itself and its space
comes back when it is deleted. A file written anywhere else costs Windows disk
permanently until the virtual disk is compacted, which needs the virtual
machine shut down and an administrator. So `--cache-dir` points at `/mnt/c`:
the weights are the only thing here large enough for the difference to matter,
and they are deleted between models.

Per-run scratch goes to shared memory rather than to either, for the same
reason. A previous stage measured about 3.5 MB of permanent virtual-disk growth
per run from scratch files that no longer existed.

What actually happened on this machine, recorded in `logs/wsl_disk_reclaim.tsv`:
the Windows drive reached 30 MB free during setup. Recovering it took deleting
regenerated trace logs, clearing package caches, removing a 2 GB binary from
the previous stage that its own documentation calls re-fetchable, and two
compactions of the virtual disk. Sparse virtual disks, which would reclaim
space automatically, are refused by this version of the platform without an
unsafe flag, so they are not used.

One consumer is outside all of this and worth watching: the Windows page file
grew by more than 2 GB during the first prediction runs and is managed by the
operating system rather than by anything here.

## Layout

```
data/
  cache/rcsb/            cached archive responses, so a second run of the set
                         builder makes no requests at all. Safe to delete; it
                         will be refetched.
  reference/             one mmCIF per target, as deposited, gzipped. The
                         structure every score is measured against.
  reference_ligands/     one SDF per target in the docking subset, cut from the
                         deposited entry by the archive's own model server so
                         the bond orders are correct.
  msa/                   alignments as returned, uncompressed, for the tools to
                         read. The compressed copies in results/msa/ are the
                         tracked ones.
  predictions/<arm>/     one structure per target per arm, gzipped. Small; the
                         confidence arrays next to them in results/ are what
                         the analysis actually reads.
  floors/                deposited entries fetched to build the two null
                         floors. One per source entry, shared between targets.
  pedv/                  the application arm: deposited entries, the construct
                         derived from them, and its predictions.
  stage1/                a clone of the docking pipeline at the commit pinned
                         in project.conf, with its own data and results
                         directories beside it so that its tables and this
                         stage's tables never mix.
  work/                  per-run scratch, when shared memory is unavailable.
                         Deleted as each run finishes; a killed run leaves a
                         directory that the next run removes.

<cache dir>/
  colabfold/params/      one parameter set, about 1.8 GB once the half this
                         benchmark never loads has been deleted. The peak
                         during extraction is twice that.
  boltz/                 the structure checkpoint and the component data.
```

## Size

Measured on this machine rather than estimated. The current figures are in
`logs/*.resources.tsv` and `results/environment/install_cost.tsv`, which is
what the README quotes; the table below is for planning.

| What | Size | Re-fetchable |
|---|---|---|
| the four environments | 4.9 GB, of which 2.5 GB is the scoring stack | yes, from `config/env_*.yml` |
| AlphaFold2 parameters, during extraction | 3.5 GB | yes |
| AlphaFold2 parameters, once trimmed | 1.8 GB | yes |
| Boltz-2 weights | about 4.2 GB after the unused checkpoint is skipped | yes |
| 150 reference structures and 15 ligands | 59 MB | yes, from the archive |
| alignments | depends entirely on the target; see `results/msa_depth.tsv` | from the committed copies in `results/msa/` |
| predictions and confidence arrays | tens of kilobytes per target per arm | yes, from stage 5 |

The gate in `scripts/00_configure.sh` refuses to write a configuration whose
projection does not fit, and the gates in the prediction stage refuse each
weight download the same way. The projection takes the larger of two moments
rather than their sum: the weights are at their largest before any result has
been written, and the results are at their largest after the weights have been
trimmed. Adding the two together describes a moment that never happens, and an
earlier version of that arithmetic refused a run that fits.

## Deleting it

Safe to remove entirely between analyses:

```bash
source project.conf
rm -rf "$DATA_DIR" "$CACHE_DIR"
```

Everything in `results/`, `figures/` and `logs/` survives, and those are what
the README and the report are built from. The alignments are the one thing
worth keeping: they are committed compressed under `results/msa/`, because the
server that produced them versions its databases independently of this code, so
the same query next year may return a different alignment and a different
structure. Re-running with those in place reproduces the inputs this run used
rather than whatever the server holds on the day.

Re-running `bash run_all.sh` rebuilds everything else. The stage stamps in
`logs/` need removing first if you want a finished stage to run again rather
than skip.
