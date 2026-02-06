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
    merge release/X.Y.Z tag: "vX.Y.Z"
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

### Version Locations (37 total)

| Category | File | Line(s) |
|----------|------|---------|
| **Python** | `src/krewlyzer/__init__.py` | 3 |
| **Python** | `pyproject.toml` | 3 |
| **Python** | `src/krewlyzer/wrapper.py` | 674 |
| **Python** | `src/krewlyzer/core/feature_serializer.py` | 54, 291 |
| **Rust** | `rust/Cargo.toml` | 3 |
| **Rust** | `rust/Cargo.lock` | Auto-updated |
| **Nextflow** | `nextflow/nextflow.config` | 171 |
| **Nextflow** | `nextflow/main.nf` | 36 |
| **Modules** | `nextflow/modules/local/krewlyzer/*/main.nf` | 2 per module |

### Quick Update Script

```bash
VERSION="X.Y.Z"
OLD_VERSION="0.3.2"  # Current version

# Python
sed -i '' "s/__version__ = \".*\"/__version__ = \"${VERSION}\"/g" src/krewlyzer/__init__.py
sed -i '' "s/version = \"${OLD_VERSION}\"/version = \"${VERSION}\"/g" pyproject.toml
sed -i '' "s/version=\"${OLD_VERSION}\"/version=\"${VERSION}\"/g" src/krewlyzer/wrapper.py
sed -i '' "s/\"${OLD_VERSION}\"/\"${VERSION}\"/g" src/krewlyzer/core/feature_serializer.py

# Rust
sed -i '' "s/version = \"${OLD_VERSION}\"/version = \"${VERSION}\"/g" rust/Cargo.toml
cd rust && cargo check && cd ..

# Nextflow
sed -i '' "s/${OLD_VERSION}/${VERSION}/g" nextflow/nextflow.config
sed -i '' "s/${OLD_VERSION}/${VERSION}/g" nextflow/main.nf
find nextflow/modules -name "main.nf" -exec sed -i '' "s/${OLD_VERSION}/${VERSION}/g" {} \;
```

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
    We do NOT publish a `:latest` tag. Always use explicit version tags like `:0.5.2`.

---

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
git merge --no-ff release/X.Y.Z -m "Release vX.Y.Z"

# Create annotated tag
git tag -a vX.Y.Z -m "Release vX.Y.Z"

# Push main and tag (triggers CI release)
git push origin main
git push origin vX.Y.Z

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
