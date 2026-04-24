# AntiFold Vendor Record

This directory contains a vendored copy of the AntiFold source tree, imported
into RFantibody so the package resolves through the standard installer and does
not depend on runtime `sys.path` injection.

## Upstream

- **Repository**: https://github.com/oxpig/AntiFold
- **Branch**: `master`
- **HEAD commit**: `789d46786624c01eb44f177ef4c0deeeb6e77469`
  (`Update README.md`, Magnus Haraldson Høie, 2025-12-09)
- **License**: BSD 3-Clause (see `LICENSE`)
- **Vendor date**: 2026-04-24
- **Vendored by**: RFantibody integration (commit `658ecdb` on `antifold-dev`)

## Scope of Vendored Content

Imported paths (relative to this directory):
- `antifold/` — the Python package (model, ESM backbone, ESM-IF1 utilities,
  IMGT helpers)
- `LICENSE`
- `README.md`

Excluded on purpose (not required by the RFantibody runner):
- `test/`, `examples/`, `notebooks/`, `output/`, `data/`
- model weights (`*.pt`) — downloaded separately into `weights/`
- CI configuration, lint configs

## Local Modifications on Top of Upstream HEAD

The vendored tree is **not** a clean snapshot of commit `789d4678`.  Two files
carry local patches for forward-compatibility with newer `biotite` releases,
where `filter_backbone` was renamed to `filter_peptide_backbone`:

```
antifold/esm/inverse_folding/util.py      | 2 call sites
antifold/esm_util_custom.py               | 3 call sites
```

The patch is mechanical:

```diff
-from biotite.structure import filter_backbone, get_chains
+from biotite.structure import filter_peptide_backbone, get_chains
 ...
-bbmask = filter_backbone(structure)
+bbmask = filter_peptide_backbone(structure)
```

No other functional changes were applied.  Verified via `diff -r` and per-file
SHA256 against the upstream working tree at vendor time.

## Numerical Equivalence Verification

Performed at vendor time (A100, CUDA 12.x, torch 2.x):

| Check | Result |
|---|---|
| Source files (excl. local patches above) | `diff -r --brief` exit 0 |
| Loaded `state_dict` SHA256 (upstream load vs vendored load) | identical |
| Forward `logits` on same seed, same input | `max |Δ| = 7.63e-06` |
| Same-path re-run (CUDA noise baseline) | `max |Δ| = 6.44e-06` |

The vendored path's deviation is within CUDA kernel non-determinism and has no
effect on downstream argmax/multinomial sampling.  Full evidence trail: see
conversation log for `antifold-dev`, 2026-04-24.

## Re-syncing with Upstream

To refresh this copy against a newer upstream commit:

```bash
# 1. Clone or pull upstream into a scratch location
git clone https://github.com/oxpig/AntiFold /tmp/antifold_upstream
cd /tmp/antifold_upstream
git checkout <new-commit>

# 2. Re-apply the biotite compatibility patch if upstream has not merged it
#    (search for `filter_backbone` — if still present, patch as above)

# 3. Mirror the antifold/ package into this directory
rsync -av --delete \
  --exclude='__pycache__' --exclude='*.pyc' \
  /tmp/antifold_upstream/antifold/ \
  include/antifold/antifold/

# 4. Refresh LICENSE / README if changed
cp /tmp/antifold_upstream/LICENSE include/antifold/
cp /tmp/antifold_upstream/README.md include/antifold/

# 5. Update this file: HEAD commit, vendor date, patch list
# 6. Re-run numerical equivalence check (see test/antifold/)
```

## Why Vendor Instead of Submodule or PyPI Dependency

- AntiFold is not published on PyPI.
- A git submodule would require network access during `pip install -e .` and
  pin users to a specific upstream revision via `.gitmodules`, which has failed
  reproducibility expectations in past containerized builds on this project.
- The upstream package is small (~750 KiB source) and rarely changes; the
  maintenance cost of vendoring is outweighed by the deterministic builds and
  the ability to carry targeted compatibility patches (e.g. the biotite API
  rename above) without waiting for upstream releases.
