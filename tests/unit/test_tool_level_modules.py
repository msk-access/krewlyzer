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

import fnmatch
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
    # A word match: `ch_extracted`, the BED channel, contains `ch_extract`.
    assert source and re.search(
        r"\bch_extract\b", source.group(1)
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


#: What run-all's module declares that no tool-level module does, deliberately, with the reason.
_RUN_ALL_ONLY = {
    # Written by run-all's own stamp and validation; tool-level has neither (the record and
    # the marker are the finalize step's).
    "*.features.json": "run-all's summary",
    "*.validation.json": "run-all's own validation",
    "*.fingerprint.json": "run-all's cohort fingerprint",
    # The single tool writes no `.sync` tables in the run's format (msk-access/krewlyzer#100).
    "*.TFBS.sync.*": "#100",
    "*.TFBS.ontarget.sync.*": "#100",
    "*.ATAC.sync.*": "#100",
    "*.ATAC.ontarget.sync.*": "#100",
}


def test_every_table_run_all_declares_a_tool_module_declares():
    """Read from the declarations, not from a run: a table that only some data produces
    (`OCF.ontarget`, which needs on-target fragments in an OCR) is never written by the
    synthetic sample, so the drift test cannot see it undeclared.
    """
    tools = [o for m, _ in nf.TOOL_MODULES for o in nf.outputs(m)]
    uncovered = [
        name
        for name in (o.glob.replace("*", "S") for o in nf.outputs("runall"))
        if not any(t.matches(name) for t in tools)
        and not any(fnmatch.fnmatchcase(name, g) for g in _RUN_ALL_ONLY)
    ]
    assert not uncovered, (
        f"run-all's module declares {uncovered}, which no tool-level module does: "
        "tool-level mode would write them and drop them from --outdir"
    )


def _items(listing: str) -> int:
    """How many items a Groovy list literal's body has, counting only top-level commas."""
    depth, count = 0, 1
    for ch in listing:
        depth += ch in "[({"
        depth -= ch in "])}"
        count += ch == "," and depth == 0
    return count


def _last_list(expression: str) -> str:
    """The body of the last `-> [ ... ]` a closure in `expression` returns."""
    start = expression.rindex("-> [") + len("-> [")
    depth = 1
    for i, ch in enumerate(expression[start:], start):
        depth += ch == "["
        depth -= ch == "]"
        if depth == 0:
            return expression[start:i]
    raise ValueError(f"unclosed list in {expression!r}")


def _channel_size(subworkflow: str, argument: str) -> int:
    """The tuple size of the channel `argument` names, or builds inline."""
    if "-> [" in argument:
        return _items(_last_list(argument))
    defined = re.search(
        rf"^\s*{argument}\s*=(.*?)(?=^\s*$|^\s*//|^\s*KREWLYZER_)",
        subworkflow,
        re.M | re.S,
    )
    if defined:
        return _items(_last_list(defined.group(1)))
    take = re.search(rf"^\s*{argument}\s*//\s*\[([^\]]*)\]", subworkflow, re.M)
    assert take, f"cannot tell what {argument} carries"
    return _items(take.group(1))


def _call_argument(subworkflow: str, process: str) -> str:
    """The first argument of `KREWLYZER_<process>(...)`, up to its top-level comma."""
    start = subworkflow.index(f"KREWLYZER_{process}(") + len(f"KREWLYZER_{process}(")
    depth = 0
    for i, ch in enumerate(subworkflow[start:], start):
        if ch in "[({":
            depth += 1
        elif ch in "])}":
            if depth == 0:
                return subworkflow[start:i].strip()
            depth -= 1
        elif ch == "," and depth == 0:
            return subworkflow[start:i].strip()
    raise ValueError(process)


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_each_process_is_given_a_tuple_of_its_input_size(module):
    """motif was given a 2-tuple where it takes three. Groovy finds that only when it runs."""
    expected = len(
        re.findall(
            r"\b(?:val|path)\(",
            re.search(r"^\s*tuple (.*)$", nf.text(module), re.M).group(1),
        )
    )
    subworkflow = nf.TOOL_LEVEL.read_text()
    given = _channel_size(subworkflow, _call_argument(subworkflow, module.upper()))
    assert given == expected, f"{module} takes a {expected}-tuple and is given {given}"


_INPUT_CHECK = nf.TOOL_LEVEL.parents[1] / "input_check" / "main.nf"


def _maf_channel(name: str) -> str:
    source = _INPUT_CHECK.read_text()
    start = source.index(f"{name} = ch_parsed")
    return source[start : source.index("\n\n", start)]


@pytest.mark.parametrize("channel", ["ch_maf_multi", "ch_maf_single"])
def test_mfsd_reads_the_dedicated_mfsd_bam_as_run_all_does(channel):
    """Run-all reads mFSD from the samplesheet's `mfsd_bam` (the duplex BAM) when given, with
    duplex weighting on. Tool-level read the main BAM, with it off."""
    body = _maf_channel(channel)
    # Not `mfsd_bai ?: bai`: an unindexed duplex BAM would be given the main BAM's index.
    assert "mfsd_bam ?: bam" in body and "mfsd_bam ? mfsd_bai : bai" in body, body
    assert "mfsd_duplex" in body, body
    assert "meta.mfsd_duplex" in nf.script("mfsd")


def test_mfsd_writes_no_table_run_all_does_not():
    """Run-all calls mFSD with `output_distributions=False`."""
    assert "--output-distributions" not in nf.script("mfsd")


def _conditions(module: str) -> dict:
    """``{flag: the condition under which the module passes it}``, for `def x_arg = c ? "--f"`."""
    found = re.findall(
        r'^\s*def \w+_arg = (.+?) \? "(--[a-z][a-z0-9-]*)', nf.script(module), re.M
    )
    return {flag: condition.strip() for condition, flag in found}


#: A tool module that passes a flag under another condition than run-all's, with the reason.
_OWN_CONDITION = {
    # Run-all switches duplex on inside the CLI when given `--mfsd-bam`; `krewlyzer mfsd` is
    # given that BAM as its only input, so the module says so.
    ("mfsd", "--duplex"): "run-all decides it in the CLI from --mfsd-bam",
}


@pytest.mark.parametrize("module", [m for m, _ in nf.TOOL_MODULES])
def test_a_flag_is_passed_under_run_alls_condition(module):
    """Passing the flag is not enough: `wps` passed `--bait-padding` only when it was truthy, so
    `--bait_padding 0` reached run-all and not the tool, which then used 50."""
    run_all = _conditions("runall")
    differ = {
        flag: (condition, run_all[flag])
        for flag, condition in _conditions(module).items()
        if flag in run_all
        and condition != run_all[flag]
        and (module, flag) not in _OWN_CONDITION
    }
    assert (
        not differ
    ), f"{module} passes {differ} (tool, run-all) under other conditions"
