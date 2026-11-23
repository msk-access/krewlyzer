"""The samplesheet's ``vcf`` column must reach mFSD in both pipeline modes.

``INPUT_CHECK`` parsed the column and built ``variants = maf ?: (vcf ?: [])``
-- in an emit, ``runall``, that the main workflow never read. Run-all mode
built its tuple from ``maf_for_filter``, which carried the MAF only, so a
VCF-only sample went to FILTER_MAF's header-only placeholder. Tool-level mode
selected mFSD samples by their MAF, so it dropped them too. Nothing failed: the
sample ran and simply had no mFSD.

Pure text analysis, like the other Nextflow tests here: Nextflow is not
installed locally or in CI. That makes these structural checks, so each one
names the behaviour it stands for.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List

import pytest

_ROOT = Path(__file__).resolve().parents[2]
_NEXTFLOW = _ROOT / "nextflow"
_INPUT_CHECK = _NEXTFLOW / "subworkflows" / "local" / "input_check" / "main.nf"
_WORKFLOW = _NEXTFLOW / "workflows" / "krewlyzer" / "main.nf"

pytestmark = [
    pytest.mark.unit,
    pytest.mark.skipif(
        not _NEXTFLOW.is_dir(), reason="nextflow/ not present in this checkout"
    ),
]


def _block(text: str, pattern: str) -> str:
    """The brace-balanced body following the first match of ``pattern``."""
    m = re.search(pattern, text)
    assert m, f"pattern not found: {pattern}"
    start = text.index("{", m.start())
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1 : i]
    raise AssertionError(f"unbalanced braces after: {pattern}")


def _statement(text: str, lhs: str) -> str:
    """``lhs = ...`` up to the next blank line: one channel definition."""
    m = re.search(rf"^\s*{re.escape(lhs)}\s*=(.*?)(?:\n\s*\n)", text, re.M | re.S)
    assert m, f"no assignment to {lhs}"
    return m.group(1)


def _strip_comments(text: str) -> str:
    return re.sub(r"//[^\n]*", "", text)


def _input_check() -> str:
    return _strip_comments(_INPUT_CHECK.read_text())


def _workflow() -> str:
    return _strip_comments(_WORKFLOW.read_text())


def _runall_mode() -> str:
    return _block(_workflow(), r"if\s*\(\s*params\.use_runall\s*\)\s*\{")


def _tool_level_mode() -> str:
    wf = _workflow()
    runall = _runall_mode()
    after = wf.index(runall) + len(runall)
    return _block(wf[after:], r"else\s*\{")


def _parse_map() -> str:
    return _block(_input_check(), r"\.splitCsv\(header:\s*true\)\s*\.map\s*\{")


def _parsed_slots() -> List[str]:
    """Names of the positions in the tuple ``ch_parsed`` carries.

    Read from the parse map's own output list, so the indexes the filters use
    are checked against the tuple rather than against a number copied here.
    """
    tuples = re.findall(r"^\s*\[(meta,[^\]]*)\]\s*$", _parse_map(), re.M)
    assert tuples, "parse map has no output tuple"
    return [s.strip() for s in tuples[-1].split(",")]


def _emit_section() -> str:
    return _input_check().split("emit:", 1)[1]


def _emits() -> List[str]:
    return re.findall(r"^\s*(\w+)\s*=", _emit_section(), re.M)


def _emitted_as(channel: str) -> str:
    m = re.search(rf"^\s*(\w+)\s*=\s*{re.escape(channel)}\s*$", _emit_section(), re.M)
    assert m, f"{channel} is not emitted by INPUT_CHECK"
    return m.group(1)


def _map_output(statement: str) -> str:
    """The tuple a ``ch_parsed.filter{}.map{}`` channel emits, whitespace-free."""
    body = _block(statement, r"\.map\s*\{")
    out = re.findall(r"^\s*(\[[^\n]*\])\s*$", body, re.M)
    assert out, "map has no output tuple"
    return re.sub(r"\s+", "", out[-1])


def test_every_input_check_emit_is_read_by_the_workflow():
    """An emit nothing reads is logic that looks wired and is not.

    ``runall`` carried the only code that routed a VCF to mFSD, so the column
    was parsed, documented and ignored. The general rule catches the next one.
    """
    wf = _workflow()
    unread = [name for name in _emits() if f"INPUT_CHECK.out.{name}" not in wf]
    assert not unread, (
        f"INPUT_CHECK emits {unread} but the main workflow never reads them. "
        "Wire them in or delete them."
    )


def test_runall_mode_gives_a_vcf_only_sample_its_vcf():
    """Run-all: a sample with a VCF and no MAF must reach RUNALL with the VCF.

    Three things have to hold. The VCF rides through ``maf_for_filter`` in the
    meta, as every other path there does (the workflow is join-free). The split
    sends it down its own branch, AFTER the MAF branch: ``branch`` takes the
    first match, which is what makes a MAF win over a VCF. And that branch goes
    straight into the RUNALL tuple, not into FILTER_MAF, which only knows MAF
    and would replace the VCF with its header-only placeholder.
    """
    maf_for_filter = re.search(r"maf_for_filter:(.*?)\[new_meta", _input_check(), re.S)
    assert maf_for_filter, "no maf_for_filter branch in INPUT_CHECK"
    assert re.search(
        r"_vcf\s*:\s*vcf\b", maf_for_filter.group(1)
    ), "maf_for_filter's meta does not carry the samplesheet vcf"

    runall = _runall_mode()
    split = _block(runall, r"INPUT_CHECK\.out\.maf_for_filter\s*\.branch\s*\{")
    labels = re.findall(r"^\s*(\w+)\s*:\s*(.+?)\s*$", split, re.M)
    names = [n for n, _ in labels]
    vcf_labels = [n for n, cond in labels if "_vcf" in cond]
    assert vcf_labels, "the MAF split has no branch for samples with a VCF"
    vcf_label = vcf_labels[0]
    maf_label = next(n for n, cond in labels if re.fullmatch(r"maf", cond))
    assert names.index(maf_label) < names.index(
        vcf_label
    ), "the VCF branch precedes the MAF branch, so a VCF would beat a MAF"

    filter_calls = re.findall(r"FILTER_MAF\(([^)]*)\)", runall)
    assert not any(
        vcf_label in c for c in filter_calls
    ), "VCF samples are sent to FILTER_MAF, which only reads MAF"

    runall_tuple = _statement(runall, "ch_runall")
    assert (
        f".{vcf_label}" in runall_tuple and "meta._vcf" in runall_tuple
    ), "the VCF branch does not feed the RUNALL tuple with meta._vcf"


def test_tool_level_mode_sends_a_vcf_only_sample_to_mfsd():
    """Tool-level: a sample with a VCF and no MAF must reach KREWLYZER_MFSD.

    Its BAM is chosen exactly as the MAF channels choose theirs, so mFSD does
    not read a different BAM depending on which variant format a row used.
    That comparison is deliberate: if the MAF channels change how they pick
    the BAM, this fails until the VCF channel follows.
    """
    ic = _input_check()
    slots = _parsed_slots()
    vcf_i, maf_i = slots.index("vcf"), slots.index("maf")

    vcf_channels = []
    for m in re.finditer(
        r"^\s*(ch_\w+)\s*=\s*ch_parsed\s*\.filter\s*\{([^}]*)\}", ic, re.M
    ):
        cond = m.group(2)
        if f"it[{vcf_i}]" in cond and re.search(
            rf"(it\[{maf_i}\]\s*==\s*null|!\s*it\[{maf_i}\])", cond
        ):
            vcf_channels.append(m.group(1))
    assert vcf_channels, (
        f"no ch_parsed channel selects samples with a vcf (it[{vcf_i}]) and no "
        f"maf (it[{maf_i}])"
    )
    vcf_channel = vcf_channels[0]
    emitted = _emitted_as(vcf_channel)

    mfsd_all = _statement(_tool_level_mode(), "ch_mfsd_all")
    assert (
        f"INPUT_CHECK.out.{emitted}" in mfsd_all
    ), f"INPUT_CHECK.out.{emitted} is not mixed into the mFSD input"

    vcf_out = _map_output(_statement(ic, vcf_channel))
    maf_out = _map_output(_statement(ic, "ch_maf_single"))
    assert vcf_out == re.sub(r"\bmaf\]$", "vcf]", maf_out), (
        f"VCF channel emits {vcf_out}, MAF channels emit {maf_out}: the BAM "
        "mFSD reads must not depend on the variant format"
    )


def test_a_row_with_both_maf_and_vcf_says_which_one_was_used():
    """MAF wins over VCF; the VCF being ignored is said once, in the log."""
    assert re.search(
        r"if\s*\(\s*(maf\s*&&\s*vcf|vcf\s*&&\s*maf)\s*\)\s*\{?\s*log\.warn",
        _parse_map(),
    ), "a row with both maf and vcf uses the MAF without saying so"
