# Shared memory

Durable context that belongs to the **project**, not to one person's machine.
Checked in, so a new contributor or a fresh agent session starts with it.

## What goes here

Things that are true about how this project is built and reasoned about, and
that no other file already says.

- **Conventions** — how work is expected to be done here.
- **Decisions** — a choice that had real alternatives, with the numbers that
  drove it, so it is not relitigated from scratch or quietly reversed.

## What does *not* go here

**Anything `.agents/rules/` already covers.** Rules are normative and read on
demand; memory is context. Restating a rule here creates a second copy that
drifts, which is the failure this repository keeps finding. Two candidates were
dropped for exactly this reason: the kreview parquet-only contract lives in
`rules/output-contract.md`, and the E1 sparsity numbers live in
`docs/features/core/fsc.md`.

**Anything the code, tests or CHANGELOG already record.** A fact with a home in
a test is better off there, where it fails when it stops being true.

**Anything local to one machine or one operator** — tooling quirks, permission
constraints, paths. Those belong in the per-user memory directory
(`~/.claude/projects/<project>/memory/`), which is not shared and not reviewed.

## Why the distinction matters

A shared file is reviewed and can be corrected by anyone; a local one cannot.
Putting a machine-specific quirk here would present one person's environment as
a project fact. Putting a project decision in local memory hides it from
everyone else.
