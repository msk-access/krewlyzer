---
name: qa-check
description: Pre-merge sweep for krewlyzer
---

# Pre-merge check

```bash
black --check src/krewlyzer/ tests/
ruff check src/krewlyzer/ tests/
mypy src/krewlyzer/ --ignore-missing-imports --no-error-summary
python -m pytest tests/ -q
cargo test --manifest-path rust/Cargo.toml --lib
cargo clippy --manifest-path rust/Cargo.toml -- -D warnings \
    -A clippy::too_many_arguments -A clippy::type_complexity
python scripts/check_output_format.py
bash scripts/check_phi_guard.sh
git status && git diff --stat
```

Assert "all pass", never a hardcoded count — exact test numbers drift.

## If you changed `rust/src/`

`maturin develop --release --manifest-path rust/Cargo.toml` before trusting the
Python suite. It imports the last-built extension, so without a rebuild it is
testing the previous binary while `cargo test` reads the current source. The
two can disagree, and the Python result is the weaker signal. A cold build is
~10 minutes.

## If a gate is inconvenient

Fix the cause. Every check here exists because something shipped without it:
`cargo test` was compiled but never run, `check_output_format.py` was
documented as CI-enforced and wired into no workflow, and two test modules
skipped everywhere except one developer's machine. A gate that is easy to
bypass is a gate that is not there.
