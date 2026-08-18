#!/usr/bin/env bash
#
# Run the linters at exactly the versions CI runs.
#
# CI installs `.[dev]`, so its black, ruff and mypy are the `==` pins in
# pyproject.toml. A developer's shell usually has whatever their environment
# happened to install, and the two drift apart silently: you run mypy, see
# green, push, and CI fails on an error your version no longer reports.
#
# That is not hypothetical. PR #71 failed CI on
#
#     error: Library stubs not installed for "markdown"  [import-untyped]
#
# which mypy 1.19.1 reports and 2.3.1 does not. Local mypy was 2.3.1 and local
# black was 26.5.1 against a 26.1.0 pin -- two of three tools adrift, with no
# indication anywhere.
#
# Versions are read from pyproject.toml rather than repeated here, so this
# cannot fall out of step with what CI installs.
# `tests/unit/test_lint_script_matches_ci.py` checks the rest: that every
# pinned tool is actually run below, with the same flags as the workflow.
#
# Usage:  bash scripts/lint.sh          # check only, like CI
#         bash scripts/lint.sh --fix    # apply black and ruff fixes
set -euo pipefail

cd "$(dirname "$0")/.."

FIX=0
[ "${1:-}" = "--fix" ] && FIX=1

pin() {
    python3 - "$1" <<'PY'
import re, sys
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
tool = sys.argv[1]
with open("pyproject.toml", "rb") as fh:
    dev = tomllib.load(fh)["project"]["optional-dependencies"]["dev"]
for spec in dev:
    m = re.fullmatch(rf"{tool}==([\w.]+)", spec.strip())
    if m:
        print(m.group(1))
        break
else:
    sys.exit(f"no == pin for {tool} in [dev]")
PY
}

# uvx fetches and caches a pinned tool without touching the active
# environment. Without it we can still run, but we say what we ran with --
# a silently different version is the whole problem.
if command -v uvx >/dev/null 2>&1; then
    run_pinned() { local tool=$1; shift; uvx "${tool}@$(pin "$tool")" "$@"; }
else
    echo "warning: uvx not found, so these run at whatever is on PATH, which" >&2
    echo "         may not be what CI uses. Install uv to pin them:" >&2
    echo "         curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
    run_pinned() {
        local tool=$1; shift
        printf '  (unpinned %s: %s)\n' "$tool" "$("$tool" --version 2>&1 | head -1)" >&2
        "$tool" "$@"
    }
fi

status=0
step() { printf '\n=== %s ===\n' "$1"; }

step "black $(pin black)"
if [ "$FIX" = 1 ]; then
    run_pinned black src/krewlyzer/ tests/ || status=1
else
    run_pinned black --check src/krewlyzer/ tests/ || status=1
fi

step "ruff $(pin ruff)"
if [ "$FIX" = 1 ]; then
    run_pinned ruff check --fix src/krewlyzer/ tests/ || status=1
else
    run_pinned ruff check src/krewlyzer/ tests/ || status=1
fi

# mypy imports krewlyzer._core, so the extension has to be built. uvx runs in
# its own environment and cannot see it; --python points mypy at the
# interpreter that can.
step "mypy $(pin mypy)"
if command -v uvx >/dev/null 2>&1; then
    uvx "mypy@$(pin mypy)" --python-executable "$(command -v python3)" \
        src/krewlyzer/ --ignore-missing-imports --no-error-summary || status=1
else
    run_pinned mypy src/krewlyzer/ --ignore-missing-imports --no-error-summary || status=1
fi

# Not pinned in pyproject -- clippy ships with the Rust toolchain. Same flags
# as the workflow: without the two -A the repo reports 11 errors CI does not.
step "cargo clippy"
cargo clippy --manifest-path rust/Cargo.toml -- \
    -D warnings \
    -A clippy::too_many_arguments \
    -A clippy::type_complexity || status=1

step "output-format call sites"
python3 scripts/check_output_format.py || status=1

if [ "$status" = 0 ]; then
    printf '\nall lint checks passed at CI versions\n'
else
    printf '\nlint failed -- see above\n' >&2
fi
exit "$status"
