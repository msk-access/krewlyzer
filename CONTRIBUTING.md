# Contributing to Krewlyzer

Thank you for your interest in contributing to Krewlyzer! This guide will help you get started.

## Development Setup

### Prerequisites

- Python 3.10+
- Rust toolchain ([rustup](https://rustup.rs/))
- C compiler (clang recommended)
- [uv](https://github.com/astral-sh/uv) (recommended)

### Clone and Install

```bash
git clone https://github.com/msk-access/krewlyzer.git
cd krewlyzer

# Create virtual environment
uv venv .venv
source .venv/bin/activate

# Install with dev dependencies
uv pip install -e ".[dev,test]"

# Verify installation
python -c "from krewlyzer import _core; print('Rust core:', _core.version())"
```

---

## Code Style

### Python

We use [Ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
# Format
ruff format src/ tests/

# Lint
ruff check src/ tests/

# Auto-fix
ruff check --fix src/ tests/
```

### Rust

We use `cargo fmt` and `cargo clippy`:

```bash
cd rust
cargo fmt
cargo clippy
```

---

## Running Tests

### All Tests

```bash
pytest
```

### With Coverage

```bash
pytest --cov=krewlyzer --cov-report=html
open htmlcov/index.html
```

### Specific Test File

```bash
pytest tests/integration/test_fsc_cli.py -v
```

### Real Data Tests

Tests in `tests/integration/test_real_data.py` use fixtures in `tests/data/fixtures/`:

```bash
pytest tests/integration/test_real_data.py -v
```

---

## Making Changes

### Workflow

1. **Fork** the repository
2. **Create a branch**: `git checkout -b feature/my-feature`
3. **Make changes** and add tests
4. **Run tests**: `pytest`
5. **Format code**: `ruff format src/ tests/`
6. **Commit**: Use conventional commits (e.g., `feat:`, `fix:`, `docs:`)
7. **Push** and open a Pull Request

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add new FSR output column
fix: resolve GC correction edge case
docs: update WPS usage example
test: add integration test for mFSD
refactor: consolidate PON loading logic
```

---

## Rust Development

### Building the Extension

```bash
cd rust
maturin develop --release
```

### Adding a New Function

1. Create or edit a module in `rust/src/`
2. Export via PyO3 in `rust/src/lib.rs`
3. Test from Python:

```python
from krewlyzer import _core
result = _core.my_module.my_function(...)
```

### Debugging Rust

```bash
RUST_LOG=debug krewlyzer extract sample.bam -r hg19.fa -o out/
```

---

## Documentation

### Local Preview

```bash
mkdocs serve
# Open http://127.0.0.1:8000
```

### Adding a Page

1. Create `.md` file in `docs/`
2. Add to `mkdocs.yml` nav section
3. Verify build: `mkdocs build --strict`

---

## Pull Request Guidelines

- **One feature per PR** - Keep changes focused
- **Include tests** - All new features need tests
- **Update docs** - Add/update documentation as needed
- **Pass CI** - All checks must pass
- **Describe changes** - Clear PR description with context

### CI Checks

- Python linting (Ruff)
- Unit and integration tests (pytest)
- Rust build (maturin)
- Docker build
- MkDocs build

---

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/msk-access/krewlyzer/issues)
- **Discussions**: [GitHub Discussions](https://github.com/msk-access/krewlyzer/discussions)

---

## License

By contributing, you agree that your contributions will be licensed under the AGPL-3.0 license.
