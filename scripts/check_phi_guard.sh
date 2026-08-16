#!/usr/bin/env bash
# Verify that the PHI/secret guard still works -- in BOTH directions.
#
# This script does not scan for PHI; gitleaks and .gitleaks.toml do that. This
# checks that the scanner would still fire if it had to. A scanner that never
# fires is worse than no scanner: it produces a false sense of safety, and
# nothing about a silent pass distinguishes "clean" from "broken".
#
# Assertions:
#   1. NEGATIVE  -- the TRACKED tree is clean. Guards against a rule that
#      false-positives and blocks every commit. Tracked-only on purpose: a
#      filesystem scan drowns in nextflow/.nextflow/ (hundreds of cached sample
#      ids, gitignored, never pushed).
#   2. POSITIVE  -- synthetic PHI is detected and EVERY rule fires. Guards
#      against a rule that silently stops matching after a typo, a bad regex or
#      an over-broad allowlist. The expected rule list is read from the config
#      rather than hardcoded, so adding a rule cannot leave it unprobed.
#   3. PLACEHOLDER -- the documented P-0000000 / C-000000-L000 placeholders are
#      NOT flagged, so docs and tests can show the shape of an identifier.
#   4. CONFIG SELF-CHECK -- .gitleaks.toml itself carries no real identifier.
#      gitleaks skips its own config, so anything pasted into a comment there is
#      invisible to the scanner and would be published unnoticed. Only this
#      assertion catches it. (It has already caught one.)
#   5. HOME PATHS -- no tracked file leaks an absolute /Users or /home path.
#      gitleaks has no rule for these and they identify a person as surely as a
#      sample id does.
#
# Exit: 0 = guard healthy, 1 = guard broken (the message says which assertion).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$REPO_ROOT/.gitleaks.toml"

fail() { echo "FAIL: $*" >&2; exit 1; }

command -v gitleaks >/dev/null 2>&1 \
    || fail "gitleaks is not installed (brew install gitleaks)."
[ -f "$CONFIG" ] || fail "missing $CONFIG"

TMP="$(mktemp -d)"
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT

scan() {  # scan <path> -> exit 1 if findings
    gitleaks detect --source "$1" --config "$CONFIG" --no-git \
        --no-banner --redact --exit-code 1 >/dev/null 2>&1
}

# ---------------------------------------------------------------------------
echo "  [1/6] tracked tree is clean ..."
# git archive gives exactly the tracked files, so gitignored local artifacts
# cannot cause a false alarm.
#
# GIT_LFS_SKIP_SMUDGE exports the 132-byte LFS pointers instead of the hydrated
# blobs. Without it `git archive` runs the smudge filter and this scan covers
# 660 MB of parquet, gzipped BED and FASTA -- 75 seconds of the script's ~100,
# to search compressed binary for a regex that could only match plaintext.
#
# What this gives up, stated plainly: an identifier *inside* an LFS payload is
# no longer scanned. That is an accepted trade. Those payloads are machine-built
# reference data (PON models, GC references, panel BEDs), none of it carrying
# free text, and every one is compressed -- gitleaks was never going to match a
# patient id in a zstd page anyway. The pointers ARE what git stores, so this
# scan now covers exactly the tracked content.
mkdir -p "$TMP/tracked"
(cd "$REPO_ROOT" && GIT_LFS_SKIP_SMUDGE=1 git archive HEAD) | tar -x -C "$TMP/tracked"
if ! scan "$TMP/tracked"; then
    echo "        findings:" >&2
    gitleaks detect --source "$TMP/tracked" --config "$CONFIG" --no-git \
        --no-banner --redact --exit-code 0 2>&1 | sed 's/^/        /' >&2 || true
    fail "assertion 1: the tracked tree contains PHI or a secret, or a rule is
      over-broad. Either remove the identifier or narrow the rule -- do not
      widen the allowlist to make this pass."
fi
echo "        OK -- no findings in tracked files"

# ---------------------------------------------------------------------------
echo "  [2/6] every rule still fires ..."
# Read the rule ids out of the config so a newly added rule cannot go unprobed.
# Not `mapfile`: macOS ships bash 3.2, so that would pass CI on Ubuntu and fail
# for every developer running this locally -- the worst way for a guard to break.
RULE_IDS=()
while IFS= read -r rule_id; do
    [ -n "$rule_id" ] && RULE_IDS+=("$rule_id")
done < <(grep -E '^[[:space:]]*id[[:space:]]*=' "$CONFIG" | sed -E 's/.*"([^"]+)".*/\1/')
[ "${#RULE_IDS[@]}" -gt 0 ] || fail "assertion 2: no rule ids found in $CONFIG"

# Probe values are assembled at runtime. Written as literals they would make
# this file itself a finding, and assertion 1 would fail on the guard.
{
  printf 'patient  P-%s-T01-XS2\n' 1234567
  printf 'cmo      C-%s-L001-d01\n' ABC123
  printf 'ssn      %s-%s-%s\n'      123 45 6789
  printf 'mrn      MRN: %s\n'       9876543
  printf 'dob      DOB: %s/%s/%s\n' 03 14 1972
} > "$TMP/probe.md"

gitleaks detect --source "$TMP" --config "$CONFIG" --no-git --no-banner --redact \
    --exit-code 0 --report-format json --report-path "$TMP/probe.json" >/dev/null 2>&1

fired="$(python3 -c "
import json,sys
try:
    data = json.load(open('$TMP/probe.json'))
except Exception:
    sys.exit('could not read the probe report')
print(' '.join(sorted({f['RuleID'] for f in data})))
")"

for rule in "${RULE_IDS[@]}"; do
    case " $fired " in
        *" $rule "*) ;;
        *) fail "assertion 2: rule '$rule' did not fire on synthetic PHI. It has
      stopped matching -- a typo, a broken regex, or an allowlist that now
      swallows it. Real identifiers of that shape would pass unnoticed.
      Rules that fired: ${fired:-none}" ;;
    esac
done
echo "        OK -- all ${#RULE_IDS[@]} rule(s) fired: $fired"

# ---------------------------------------------------------------------------
echo "  [3/6] documented placeholders are allowed ..."
mkdir -p "$TMP/placeholder"
{
  echo 'Example output: P-0000000-T01-XS1.FSC.gene.parquet'
  echo 'Example sample: C-000000-L000-d01.bam'
} > "$TMP/placeholder/doc.md"
scan "$TMP/placeholder" \
    || fail "assertion 3: the documented placeholders are being flagged, so docs
      and tests cannot show the shape of an identifier. Check [allowlist] in
      $CONFIG."
echo "        OK -- P-0000000 and C-000000-L000 pass"

# ---------------------------------------------------------------------------
echo "  [4/6] the config itself carries no identifier ..."
# gitleaks refuses to scan its own config, so renaming it is the only way in.
mkdir -p "$TMP/configcopy"
cp "$CONFIG" "$TMP/configcopy/config-copy.txt"
scan "$TMP/configcopy" \
    || fail "assertion 4: $CONFIG contains a real identifier, most likely in a
      comment. gitleaks skips its own config, so nothing else would ever catch
      it and it would be published. Replace it with a placeholder."
echo "        OK -- no identifiers in the config"

# ---------------------------------------------------------------------------
echo "  [5/6] no absolute home paths in tracked files ..."
# LFS payloads excluded for the same reason as assertion 1: grepping 660 MB of
# compressed reference data for "/Users/" costs 13 seconds and cannot match.
(cd "$REPO_ROOT" && git lfs ls-files -n 2>/dev/null || true) | sort > "$TMP/lfs-paths"
LEAKED="$(cd "$REPO_ROOT" && git ls-files \
    | { [ -s "$TMP/lfs-paths" ] && grep -vxF -f "$TMP/lfs-paths" || cat; } \
    | tr '\n' '\0' \
    | xargs -0 grep -lE '(/Users/|/home/)[A-Za-z0-9._-]+' 2>/dev/null || true)"
if [ -n "$LEAKED" ]; then
    echo "        offending files:" >&2
    printf '        %s\n' $LEAKED >&2
    fail "assertion 5: a tracked file contains an absolute home path. That names
      a person and hardcodes one machine. Use a relative path, an environment
      variable, or a CLI argument."
fi
echo "        OK -- no home paths in tracked files"

# ---------------------------------------------------------------------------
echo "  [6/6] the hooks on disk are the hooks in the repo ..."
DRIFT=""
for hook in $(cd "$REPO_ROOT" && git ls-files .githooks/); do
    disk="$REPO_ROOT/$hook"
    [ -f "$disk" ] || { DRIFT="$DRIFT
        $hook -- MISSING from the working tree"; continue; }
    # Worktree vs INDEX, not vs HEAD. A hook you are deliberately editing is
    # staged, and should pass; an overwrite by `git lfs install` is not staged,
    # and should fail. Comparing against HEAD would flag every legitimate edit
    # and train people to ignore this assertion.
    if ! (cd "$REPO_ROOT" && git diff --quiet -- "$hook") 2>/dev/null; then
        DRIFT="$DRIFT
        $hook -- differs from the staged version"
    fi
done

# The guard is only armed if pre-push still delegates to git-lfs. Without this,
# `git lfs install` has a reason to run, and installing is what replaces the
# hook -- so a missing delegation is a latent version of assertion 6 failing.
grep -q 'git lfs pre-push' "$REPO_ROOT/.githooks/pre-push" 2>/dev/null \
    || DRIFT="$DRIFT
        .githooks/pre-push -- no longer delegates to git-lfs"

if [ -n "$DRIFT" ]; then
    echo "        drift:$DRIFT" >&2
    fail "assertion 6: a git hook on disk does not match the one in the repo.
      This is how the guard was lost during 0.9.0: checking out a commit that
      predates .githooks/ removes the hooks, and the next \`git lfs install\`
      writes its own three-line pre-push in place of the PHI scan. Nothing
      reports it -- the push simply stops being checked. Restore with
      \`git checkout -- .githooks/\` and confirm \`git config core.hooksPath\`
      is still .githooks."
fi
echo "        OK -- hooks match the repo and pre-push still chains git-lfs"

echo "PHI guard healthy."
