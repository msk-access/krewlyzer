# Release Guide

This guide documents the process for releasing new versions of Krewlyzer following Git Flow.

---

## Prerequisites

- Git LFS installed and configured
- Access to push to `origin`
- All tests passing on develop branch

!!! important
    **Version Format**: Use `0.5.2` (no `v` prefix) everywhere - code, filenames, and git tags.

---

## Git Flow Overview

```mermaid
gitGraph
    commit id: "develop"
    branch release/X.Y.Z
    commit id: "bump version"
    commit id: "update CHANGELOG"
    checkout main
    merge release/X.Y.Z tag: "X.Y.Z"
    checkout develop
    merge release/X.Y.Z
```

---

## Phase 1: Create Release Branch

```bash
# Ensure you're on latest develop
git checkout develop
git pull origin develop

# Create release branch
git checkout -b release/X.Y.Z

# Verify
git branch --show-current
```

---

## Phase 2: Update Version Files

### Version Locations

Line numbers are deliberately omitted: the previous table pointed at
`wrapper.py:674` and `feature_serializer.py:54,291`, none of which held a
version by the time anyone read it. A stale pointer in a release checklist is
worse than no pointer, because it invites editing the wrong line.

| Category | File | Notes |
|----------|------|-------|
| **Python** | `src/krewlyzer/__init__.py` | `__version__` — the single source of truth |
| **Python** | `pyproject.toml` | packaging metadata |
| **Rust** | `rust/Cargo.toml` | crate version |
| **Rust** | `rust/Cargo.lock` | auto-updated by `cargo check` |
| **Nextflow** | `nextflow/nextflow.config` | container tag |
| **Nextflow** | `nextflow/main.nf` | container tag |
| **Modules** | `nextflow/modules/local/krewlyzer/*/main.nf` | 2 per module (container + `versions.yml`) |

Everything on the Python side other than `__init__.py` imports `__version__`.
`wrapper.py` and `core/feature_serializer.py` used to keep their own copies;
`tests/unit/test_version_stamp.py` fails if one reappears.

### Quick Update Script

```bash
VERSION="X.Y.Z"
OLD_VERSION="A.B.C"  # Current version before bump

# Python
# src/krewlyzer/__init__.py is the single source of truth for the Python side.
# wrapper.py and core/feature_serializer.py used to carry their own copies and
# needed their own sed lines; they now import __version__, so there is nothing
# to substitute. tests/unit/test_version_stamp.py fails if a copy reappears.
sed -i '' "s/__version__ = \".*\"/__version__ = \"${VERSION}\"/g" src/krewlyzer/__init__.py
sed -i '' "s/version = \"${OLD_VERSION}\"/version = \"${VERSION}\"/g" pyproject.toml

# Rust
sed -i '' "s/version = \"${OLD_VERSION}\"/version = \"${VERSION}\"/g" rust/Cargo.toml
cd rust && cargo check && cd ..

# Nextflow
sed -i '' "s/${OLD_VERSION}/${VERSION}/g" nextflow/nextflow.config
sed -i '' "s/${OLD_VERSION}/${VERSION}/g" nextflow/main.nf
find nextflow/modules -name "main.nf" -exec sed -i '' "s/${OLD_VERSION}/${VERSION}/g" {} \;
```

!!! warning "Verify the bump landed"
    Every command above depends on `OLD_VERSION` being exactly right, and
    nothing in `sed` reports a pattern that matched nothing. Run this
    afterwards — it fails if `__init__.py`, `pyproject.toml` and
    `rust/Cargo.toml` disagree, or if any module has reintroduced a version
    literal of its own:

    ```bash
    pytest tests/unit/test_version_stamp.py -q
    ```

    Then confirm nothing is left behind:

    ```bash
    git grep -n "${OLD_VERSION}" -- . ':!CHANGELOG.md' ':!docs'
    ```

    Remaining hits in `CHANGELOG.md` and `docs/` are expected: those are
    historical references to how an older release behaved, and rewriting them
    would falsify the record.

---

## Phase 2.5: Update Documentation Versions

Docker image versions are referenced in documentation files:

| File | Version Location |
|------|------------------|
| `docs/getting-started/installation.md` | Docker/Singularity pull commands |
| `docs/getting-started/quickstart.md` | Docker pull example |
| `docs/nextflow/examples.md` | Container image references |

### Update Script

```bash
OLD_VERSION="0.5.1"
VERSION="X.Y.Z"

# Update installation docs
sed -i '' "s/${OLD_VERSION}/${VERSION}/g" docs/getting-started/installation.md
sed -i '' "s/${OLD_VERSION}/${VERSION}/g" docs/getting-started/quickstart.md

# Verify no :latest tags remain (we don't publish :latest)
grep -r ":latest" docs/ && echo "WARNING: :latest tags found!" || echo "✓ No :latest tags"

# Verify changes
grep -n "ghcr.io/msk-access/krewlyzer" docs/getting-started/*.md
```

!!! warning "No :latest Tag"
    We do NOT publish a `:latest` tag. Always use explicit version tags (e.g., `:0.5.3`).
    Replace `X.Y.Z` with the version from [releases](https://github.com/msk-access/krewlyzer/releases).

---

## Phase 2.6: Regenerate the aggregated documentation

`krewlyzer_all_docs.md` is a single-file concatenation of `docs/`, generated
rather than hand-maintained. Any doc edit in the release — including the
version bumps in Phase 2.5 — leaves it stale.

```bash
python scripts/build_all_docs.py
```

CI runs `--check` and fails if it is out of date, so this cannot be skipped
silently.

---

## Phase 2.7: Stamp the bundled PONs

The version-update script in Phase 2 uses `sed`, and a PON is a Parquet file —
so the models are the one place a version literal does **not** get updated by
it. They record the version of whatever built them, which is a `develop`
checkout still reporting the previous release.

```bash
krewlyzer stamp-pon src/krewlyzer/data/pon/GRCh37/*/*.parquet --version X.Y.Z
krewlyzer validate-pon src/krewlyzer/data/pon/GRCh37/*/*.parquet
```

`stamp-pon` refuses to stamp a model that fails `validate-pon`, so a broken
model cannot be blessed by this step. Re-run `validate-pon` afterwards anyway:
the file that ships should be the file that was checked.

Commit the restamped models with `git lfs push --all` **before** pushing the
branch, or the pointers land without the objects.

## Phase 2.8: Does this release change what a PON means?

Almost always **no**, and then there is nothing to do here.

`MIN_PON_VERSION` in `src/krewlyzer/pon/provenance.py` is a compatibility
floor, not the package version. It is a tuple — `(0, 9, 0)` — precisely so the
`sed` in Phase 2 cannot move it, and so a PON stamped 0.9.0 keeps working at
krewlyzer 1.0 and beyond. Ordinary releases leave it alone.

Raise it **only** when this release changes what an existing feature *means*,
such that a PON built before it would score samples against a different
quantity. The 0.9.0 examples:

- `wps_background` held a hardcoded `167.0 / 5.0` for every group
- six σ floors turned "no spread measured" into a divisor
- region-MDS was fitted over 65–400 bp while samples are measured over
  65–1000 bp — a median bias of +1.15 σ in every gene

Adding a *new* feature or block is not that: an older PON simply lacks it, and
`validate-pon`'s packing-list check reports the absence.

If you do raise it, every bundled PON must be rebuilt — not merely re-stamped.
Stamping a model that predates the change would launder exactly the
incompatibility the floor exists to catch.

## Phase 3: Update CHANGELOG

Add new entry at the top of `CHANGELOG.md`:

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added
- Feature descriptions

### Changed
- Breaking changes and modifications

### Fixed
- Bug fixes

### Documentation
- Doc updates
```

---

## Phase 4: Verify and Commit

```bash
# Run tests
pytest tests/ -v --tb=short

# Verify Rust compiles
cd rust && cargo check && cd ..

# Check version
python -c "from krewlyzer import __version__; print(__version__)"

# Commit
git add -A
git commit -m "chore: bump version to X.Y.Z"

# Push release branch for review
git push -u origin release/X.Y.Z
```

---

## Phase 5: Finalize Release

After review and approval:

```bash
# Merge to main
git checkout main
git pull origin main
git merge --no-ff release/X.Y.Z -m "Release X.Y.Z"

# Create annotated tag
git tag -a X.Y.Z -m "Release X.Y.Z"

# Push main and tag (triggers CI release)
git push origin main
git push origin X.Y.Z

# Merge back to develop
git checkout develop
git merge --no-ff release/X.Y.Z -m "Merge release/X.Y.Z back to develop"
git push origin develop

# Delete release branch
git branch -d release/X.Y.Z
git push origin --delete release/X.Y.Z
```

---

## Git LFS

Large files are tracked via Git LFS (see `.gitattributes`):

```
src/krewlyzer/data/**/*.gz filter=lfs diff=lfs merge=lfs -text
src/krewlyzer/data/**/*.parquet filter=lfs diff=lfs merge=lfs -text
src/krewlyzer/data/**/*.bed filter=lfs diff=lfs merge=lfs -text
```

Ensure Git LFS is installed before cloning:

```bash
git lfs install
git clone https://github.com/msk-access/krewlyzer.git
```

---

## Nextflow Module Locations

Each module has 2 version references:

| Module | Container Line | versions.yml Line |
|--------|----------------|-------------------|
| build_pon | 13 | 63 |
| extract | 13 | 68 |
| fsc | 14 | 67 |
| fsd | 13 | 60 |
| fsr | 13 | 63 |
| mfsd | 13 | 62 |
| motif | 13 | 55 |
| ocf | 13 | 63 |
| region_entropy | 14 | 73 |
| region_mds | 14 | 79 |
| runall | 18 | 183 |
| uxm | 13 | 55 |
| wps | 14 | 76 |
