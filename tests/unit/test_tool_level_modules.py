"""The tool-level Nextflow modules give each tool what run-all gives that stage, and publish.

`--use_runall false` runs each tool as its own process. Measured on 22 Nov, it could not produce a
sample: eight modules required a `.tsv` output under the Parquet default, `motif` was handed the
BED where it takes a BAM, nothing was published, and five tools were given no assay, the
pipeline-wide targets or PON instead of the sample's, or no `--skip-pon`. So a panel sample's
tables, had any been written, would not have been run-all's (invariant #6).

Text checks (Nextflow is not installed; `tests/nextflow_modules.py`). What the tools actually write
against these declarations is `tests/invariants/test_tool_level_outputs.py`.
"""

from __future__ import annotations

import re

import pytest
import typer.main

import nextflow_modules as nf
from krewlyzer.cli import app

pytestmark = pytest.mark.unit

_CLI = typer.main.get_command(app)


def _options(command: str) -> set:
    params = _CLI.commands[command].params
    return {
        o
        for p in params
        for o in (*getattr(p, "opts", []), *getattr(p, "secondary_opts", []))
        if o.startswith("--")
    }


@pytest.mark.parametrize("module, command", nf.TOOL_MODULES, ids=lambda v: str(v))
def test_every_option_run_all_gives_a_stage_the_module_gives_its_tool(module, command):
    """Each option run-all's module passes and the tool takes, the tool's module passes too."""
    expected = nf.flags(nf.script("runall")) & _options(command)
    missing = sorted(expected - nf.flags(nf.script(module)))
    assert not missing, (
        f"{module} does not pass {missing}, which run-all's module passes and "
        f"`krewlyzer {command}` takes: the tool-level table would differ from run-all's"
    )


#: Settings run-all takes per sample (`INPUT_CHECK` resolves them from the samplesheet and the
#: assay), never pipeline-wide.
_PIPELINE_WIDE = (
    "params.targets",
    "params.pon_model",
    "params.wps_anchors",
    "params.assay",
)


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_a_sample_setting_comes_from_the_sample(module):
    used = [
        p for p in _PIPELINE_WIDE if re.search(rf"{re.escape(p)}\b", nf.script(module))
    ]
    assert not used, (
        f"{module} reads {used}: run-all takes these per sample, so a sample whose samplesheet "
        "row or assay names another one would be computed differently"
    )


def _assay_rule(module: str) -> str:
    found = re.findall(r"^\s*def assay_arg = (.*)$", nf.script(module), re.M)
    return found[0].strip() if found else ""


@pytest.mark.parametrize(
    "module",
    [m for m, c in nf.TOOL_MODULES if "--assay" in _options(c)],
)
def test_the_assay_is_given_as_run_all_gives_it(module):
    """`meta.assay`, lower-cased, and none for UNKNOWN or WGS."""
    assert _assay_rule("runall"), "run-all's module no longer builds assay_arg"
    assert _assay_rule(module) == _assay_rule("runall")


def test_motif_is_given_the_bam():
    """`krewlyzer motif` reads a BAM. It was handed the BED, in a tuple one item short."""
    assert re.search(
        r"tuple val\(meta\), path\(bam\), path\(bai\)", nf.text("motif")
    ), "motif's input is no longer the BAM"
    subworkflow = nf.TOOL_LEVEL.read_text()
    call = re.search(r"KREWLYZER_MOTIF\(\s*(\w+)", subworkflow)
    assert call, "tool_level no longer calls KREWLYZER_MOTIF"
    source = re.search(rf"^\s*{call.group(1)}\s*=\s*(.*)$", subworkflow, re.M)
    assert source and "ch_extract" in source.group(
        1
    ), f"KREWLYZER_MOTIF is given {call.group(1)}, which is not built from the BAM channel"


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_no_required_output_names_one_format(module):
    """Under `--output-format parquet` (the default) a required `*.X.tsv` fails the task."""
    single = [
        o.glob
        for o in nf.outputs(module)
        if not o.optional and re.search(r"\.tsv(\.gz)?$", o.glob)
    ]
    assert not single, f"{module} requires {single}, which only a TSV run writes"


def test_every_emit_the_subworkflow_reads_is_declared():
    """A renamed emit is found only when Nextflow runs; Groovy resolves `.out.x` at runtime."""
    declared = {
        m: {o.emit for o in nf.outputs(m)} | {"versions"} for m, _ in nf.TOOL_MODULES
    }
    reads = re.findall(r"KREWLYZER_(\w+)\.out\.(\w+)", nf.TOOL_LEVEL.read_text())
    assert reads, "tool_level reads no module output"
    for process, emit in reads:
        assert (
            emit in declared[process.lower()]
        ), f"tool_level reads {process}.out.{emit}, which the module does not declare"


_RUN_ALL_PUBLISH = nf.publish_dir("runall")


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES if m != "extract"])
def test_every_tool_publishes_into_the_sample_folder_as_run_all_does(module):
    """No tool-level module published: a green run left `--outdir` empty."""
    assert _RUN_ALL_PUBLISH, "run-all's module no longer publishes"
    assert nf.publish_dir(module) == _RUN_ALL_PUBLISH


def test_extract_publishes_everything_but_its_record():
    """`{sample}.metadata.parquet` is downstream's completion marker. Extraction's record has
    no `krewlyzer_version` or `pon_applied`, and would appear before any feature table.
    """
    line = nf.publish_dir("extract")
    assert line.startswith(_RUN_ALL_PUBLISH.split("saveAs")[0]), line
    assert "filename.contains('.metadata.') ? null" in line, line
    assert any(".metadata." in o.glob for o in nf.outputs("extract"))


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES if m != "extract"])
def test_no_tool_declares_the_record(module):
    """Only extraction's record exists in tool-level mode, and it is not the marker."""
    assert not [o.glob for o in nf.outputs(module) if ".metadata." in o.glob]


def _stub(module: str) -> str:
    source = nf.text(module)
    return source[source.index("\n    stub:") :]


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_the_stub_writes_every_required_output(module):
    """`nextflow -stub-run` fails a task whose stub leaves a required output unwritten, and the
    stub is the only end-to-end check of the wiring that runs without the data."""
    files = [
        "S." + name for name in re.findall(r"\$\{prefix\}\.([\w.]+)", _stub(module))
    ]
    missing = [
        o.glob
        for o in nf.outputs(module)
        if not o.optional and not any(o.matches(f) for f in files)
    ]
    assert not missing, f"{module}'s stub writes nothing matching {missing}"


def test_region_mds_runs_only_for_a_sample_with_an_assay():
    """run-all runs Region-MDS only with an assay (its gene BED). Without `--assay` the tool falls
    back to the WGS gene BED, so tool-level wrote `MDS.exon/gene` that run-all never writes.
    """
    subworkflow = nf.TOOL_LEVEL.read_text()
    call = subworkflow[subworkflow.index("KREWLYZER_REGION_MDS(") :]
    call = call[: call.index("\n    )")]
    assert ".filter" in call and "'UNKNOWN'" in call and "'WGS'" in call, call
