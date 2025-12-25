# Contributing

See the full [Contributing Guide](https://github.com/msk-access/krewlyzer/blob/main/CONTRIBUTING.md) for detailed instructions.

## Quick Start

```bash
# Clone and setup
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev,test]"

# Make changes, test, and format
pytest
ruff format src/ tests/

# Commit and PR
git checkout -b feature/my-feature
git commit -m "feat: add new feature"
git push origin feature/my-feature
```

## Key Commands

| Task | Command |
|------|---------|
| Run tests | `pytest` |
| Format Python | `ruff format src/ tests/` |
| Lint Python | `ruff check src/ tests/` |
| Build Rust | `cd rust && maturin develop --release` |
| Preview docs | `mkdocs serve` |
| Build docs | `mkdocs build --strict` |
