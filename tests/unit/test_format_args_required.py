"""A writer or reader of sample outputs takes its format from the caller, never from a default.

Invariant #6: single-tool and run-all output must match. A ``output_format="tsv"`` default lets a
caller that forgets the argument write -- or read -- TSV beside a Parquet run, silently. With no
default the signature enforces it: mypy reports the omission in every annotated caller, and any
path a test runs raises ``TypeError``. `scripts/check_output_format.py` covers the older writers
that still default.
"""

import importlib
import inspect

import pytest

REQUIRED = {
    "krewlyzer.core.sample_processor": (
        "write_motif_tables",
        "write_motif_outputs",
        "process_sample",
    ),
    "krewlyzer.core.output_provenance": (
        "record_stage_outcomes",
        "record_pon_used",
        "refuse_a_second_pon",
        "record_path",
        "recorded_status",
        "record_names_pon",
        "stamp_metadata",
        "refuse_a_different_configuration",
        "recorded_settings",
    ),
    "krewlyzer.uxm": ("uxm_stage",),
    "krewlyzer.mfsd": ("mfsd_stage",),
    "krewlyzer.region_mds": ("region_mds_stage",),
    "krewlyzer.core.pon_integration": ("load_requested_pon",),
    "krewlyzer.core.unified_processor": ("run_features",),
    "krewlyzer.core.output_utils": ("convert_rust_tsv",),
}

CASES = [(module, name) for module, names in REQUIRED.items() for name in names]


@pytest.mark.parametrize("module,name", CASES, ids=[n for _, n in CASES])
@pytest.mark.parametrize("arg", ["output_format", "compress"])
def test_format_argument_is_keyword_only_with_no_default(module, name, arg):
    params = inspect.signature(
        getattr(importlib.import_module(module), name)
    ).parameters
    assert arg in params, f"{name} takes no {arg}"
    p = params[arg]
    assert p.kind is inspect.Parameter.KEYWORD_ONLY, f"{name}({arg}) is {p.kind.name}"
    assert (
        p.default is inspect.Parameter.empty
    ), f"{name}({arg}) defaults to {p.default!r}"
