# Local cohort gate

Opt-in checks against a real scored cohort. **Not run in CI, and no data from a
cohort is ever committed here.**

```bash
KREWLYZER_TEST_CORPUS=/path/to/scored/cohort python -m pytest tests/real_data -q
```

Unset the variable and everything skips.

## Why this exists

The unit suite proves a scorer *can* run. Only a cohort proves it *did* — and
only a cohort can distinguish a metric that varies from one that is constant.
Four PON blocks were built, shipped and consumed by nothing for a year;
`test_pon_coverage.py` marks those `xfail(strict=True)`, so wiring them up
turns the test red and forces the marker off.

## PHI

Sample directories are named for the patient. Nothing here may put an
identifier in a test name, a parameter id, an assertion message, or a written
artifact. Address samples positionally, or use `sample_label()` for a stable
non-reversible handle.
