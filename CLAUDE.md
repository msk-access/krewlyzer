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
- **Hooks** are git-level and tool-agnostic, in `.githooks/`, so they apply
  whether a human, Claude, Cursor or Codex drives the commit. They are opt-in
  per clone — `git config core.hooksPath .githooks` — because git offers no way
  to enable them from a checkout. CI repeats every check as the backstop.
- **The extension can be stale.** `src/krewlyzer/_core*.so` is gitignored and
  built by `maturin develop --release`. Python tests import whatever was built
  last, so after any `rust/src/` change the suite is testing the *old* binary
  until you rebuild. `cargo test` reads the source directly and is unaffected —
  which means the two can disagree, and Python passing is the weaker signal.
