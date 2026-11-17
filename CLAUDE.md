# krewlyzer — Claude Code entry

> Canonical standing context is `AGENTS.md` (cross-tool). This file imports it so
> Claude Code loads the same rules every other tool reads. Keep Claude-only notes
> below the import; put shared rules in `AGENTS.md`.

@AGENTS.md

## Claude-only notes

- **Rules** live in `.agents/rules/*.md` and are read on demand, not every turn.
  `architecture.md` has the Rust/Python boundary table; `development.md` has the
  build and QC commands; `output-contract.md` has what downstream relies on;
  `validation-gates.md` has the order to run the gates in and what each one is
  blind to.
- **Shared memory** is in `.agents/memory/` — conventions, the modelling
  philosophy, and dated decision records. It holds only what no rule, test or
  doc already states; `.agents/memory/README.md` says what belongs there and
  what belongs in your own per-user memory instead.
- **Hooks** are git-level and tool-agnostic, in `.githooks/`, so they apply
  whether a human, Claude, Cursor or Codex drives the commit. They are opt-in
  per clone — `git config core.hooksPath .githooks` — because git offers no way
  to enable them from a checkout. CI repeats every check as the backstop.
- **The extension can be stale.** `src/krewlyzer/_core*.so` is gitignored and
  built by `maturin develop --release`. Python tests import whatever was built
  last, so after any `rust/src/` change the suite is testing the *old* binary
  until you rebuild. `cargo test` reads the source directly and is unaffected —
  which means the two can disagree, and Python passing is the weaker signal.

## Behavioural guidelines (Karpathy)

From [multica-ai/andrej-karpathy-skills](https://github.com/multica-ai/andrej-karpathy-skills)
(`skills/karpathy-guidelines/SKILL.md`, MIT), derived from
[Karpathy's notes on LLM coding pitfalls](https://x.com/karpathy/status/2015883857489522876).
Condensed rather than copied whole, because two of the four already exist in
`AGENTS.md` in a stronger form — noted inline so neither reads as the whole rule.

They bias toward caution over speed; for trivial tasks, use judgement.

1. **Think before coding.** State assumptions explicitly. Where several
   readings exist, present them rather than silently picking one. If a simpler
   approach exists, say so. If something is unclear, stop and name it.

2. **Simplicity first.** The minimum code that solves the problem: no features
   beyond the ask, no abstractions for single-use code, no configurability
   nobody requested, no error handling for impossible states. If 200 lines
   could be 50, rewrite it.

3. **Surgical changes.** Touch only what the request requires; don't "improve"
   adjacent code, comments or formatting, and match existing style even where
   you'd choose differently. Remove the imports and helpers *your* change
   orphaned — mention pre-existing dead code rather than deleting it.
   Every changed line should trace to the request.

4. **Goal-driven execution.** Turn a task into a verifiable goal before
   starting — "add validation" becomes "write tests for invalid inputs, then
   make them pass". For multi-step work, state the plan with the check for
   each step.

Where this repo already goes further, and the stricter rule wins:

- On verification (4), `AGENTS.md` names the specific gates and
  `rules/validation-gates.md` says what each one is blind to; invariant #1
  additionally requires asserting that *two different inputs produce two
  different outputs*, which "write a test and make it pass" does not imply.
- On surgical changes (3), the Git Flow convention requires docs and CHANGELOG
  in the **same commit** as the change they describe, so the minimal diff is
  not always the smallest one.
