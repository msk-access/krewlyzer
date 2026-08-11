"""Numeric constants that documentation, code and this registry must agree on.

Most of the audit that produced the output-contract gate was not finding broken
code. It was finding places where a document stated a number the implementation
had since changed -- a tolerance of 20 documented as 50, a search band of
140-200 documented as 150-250, size channels that had been retired a year
earlier. Each was individually harmless and collectively meant the docs could
not be trusted to describe the tool.

Prose cannot enforce that. A test can: every entry below names a constant, the
value it is supposed to have, and where the implementation defines it. Change
the code without changing the registry and the test fails; change the registry
without changing the code and it fails too.

Two registries, because reality has warts:

* :data:`CLAIMS` -- constants that must hold. The ordinary case.
* :data:`KNOWN_DIVERGENCES` -- places where two parts of the codebase
  disagree *today*, recorded with why it matters. These assert that the
  disagreement is still present, so resolving one fails the test and forces
  the entry to be removed rather than left to rot.
"""

from __future__ import annotations

import importlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Union

_RUST_SRC = Path(__file__).resolve().parents[3] / "rust" / "src"


class ConstantNotFound(LookupError):
    """A referenced constant is absent.

    Raised rather than skipped: a renamed constant is exactly the drift this
    registry exists to catch, and silently passing would defeat it.
    """


@dataclass(frozen=True)
class RustConst:
    """A `const NAME: type = value;` in a Rust source file."""

    file: str
    name: str

    def resolve(self) -> str:
        path = _RUST_SRC / self.file
        if not path.exists():
            raise ConstantNotFound(f"{self.file} does not exist")
        pattern = re.compile(
            rf"const\s+{re.escape(self.name)}\s*:\s*[A-Za-z0-9_:<>]+\s*=\s*([^;]+);"
        )
        match = pattern.search(path.read_text())
        if match is None:
            raise ConstantNotFound(f"const {self.name} not found in {self.file}")
        return match.group(1).strip()


@dataclass(frozen=True)
class RustPattern:
    """A literal that must appear in a Rust file, for things that are not consts.

    The FSC channel boundaries are inline comparisons rather than named
    constants, so there is nothing to look up -- only a shape to pin.
    """

    file: str
    pattern: str
    description: str

    def resolve(self) -> str:
        path = _RUST_SRC / self.file
        if not path.exists():
            raise ConstantNotFound(f"{self.file} does not exist")
        found = re.findall(self.pattern, path.read_text())
        if not found:
            raise ConstantNotFound(
                f"pattern for {self.description} not found in {self.file}"
            )
        return ",".join(found)


@dataclass(frozen=True)
class PyConst:
    """A module-level constant in the Python package."""

    module: str
    name: str

    def resolve(self) -> str:
        try:
            module = importlib.import_module(self.module)
        except ImportError as exc:
            # A moved or renamed module is drift too; surface it as such rather
            # than letting an ImportError read as an environment problem.
            raise ConstantNotFound(f"cannot import {self.module}: {exc}") from exc
        if not hasattr(module, self.name):
            raise ConstantNotFound(f"{self.module}.{self.name} not found")
        return str(getattr(module, self.name))


ImplRef = Union[RustConst, RustPattern, PyConst]


@dataclass(frozen=True)
class Claim:
    id: str
    value: str
    impl: ImplRef
    why: str


@dataclass(frozen=True)
class Divergence:
    """Two places that should agree and do not."""

    id: str
    left: ImplRef
    left_value: str
    right: ImplRef
    right_value: str
    why: str


# ---------------------------------------------------------------------------
# Constants that must hold
# ---------------------------------------------------------------------------

CLAIMS: Tuple[Claim, ...] = (
    # -- WPS background periodicity ----------------------------------------
    # The whole NRL family was constant before these values changed. Pinning
    # them means a future edit to the window cannot silently re-degenerate the
    # metric: the DFT grid is n_fft * BIN_SIZE / i, so the search band and the
    # padding together decide how many periods are even reachable.
    # Pinned as the expression rather than 2000: the derivation is the point,
    # and FLANK_BP and ALU_BODY_BP are pinned separately below, so the
    # effective 850*2 + 300 = 2000 is fully constrained either way.
    Claim(
        "wps.background.profile_length",
        "Self::FLANK_BP * 2 + Self::ALU_BODY_BP",
        RustConst("wps.rs", "PROFILE_LENGTH"),
        "300bp of Alu body alone admitted exactly one period inside the "
        "nucleosomal band, pinning nrl_bp at 150.0 for every sample",
    ),
    Claim(
        "wps.background.flank_bp",
        "850",
        RustConst("wps.rs", "FLANK_BP"),
        "spans ~10 nucleosome repeats either side of the Alu body",
    ),
    Claim(
        "wps.background.alu_body_bp",
        "300",
        RustConst("wps.rs", "ALU_BODY_BP"),
        "length of the Alu consensus the profile is stacked on",
    ),
    Claim(
        "wps.background.pad_factor",
        "8",
        RustConst("wps.rs", "PAD_FACTOR"),
        "zero-padding sets the DFT grid spacing; without it the reachable "
        "periods are too coarse to resolve 190bp",
    ),
    Claim(
        "wps.expected_nrl_bp",
        "190.0",
        RustConst("wps.rs", "EXPECTED_NRL_BP"),
        "healthy chromatin nucleosome repeat length (Snyder et al. 2016)",
    ),
    Claim(
        "wps.nrl_tolerance_bp",
        "50.0",
        RustConst("wps.rs", "TOLERANCE_BP"),
        "documented tolerance; the code used 20 while docs said 50, which "
        "drove adjusted_score to 0 for every sample",
    ),
    # -- GC correction ------------------------------------------------------
    Claim(
        "gc.length_bins",
        "188",
        RustConst("gc_correction.rs", "NUM_BINS"),
        "5bp bins spanning 60-999bp; see the divergence below",
    ),
    # -- FSD ------------------------------------------------------------------
    Claim(
        "fsd.bin_range",
        "65..400",
        RustPattern(
            "fsd.rs",
            r"for s in \((\d+\.\.\d+)\)\.step_by\(5\)",
            "the FSD histogram's bin range",
        ),
        "the histogram bins [65, 400) in 67 steps of 5, while the length "
        "*filter* admits 65-1000. meaning.py quoted the filter range as the "
        "histogram range, which makes `total` look like it covers fragments "
        "it excludes entirely",
    ),
    # -- mFSD ---------------------------------------------------------------
    Claim(
        "mfsd.min_for_ks",
        "2",
        RustConst("mfsd.rs", "MIN_FOR_KS"),
        "mfsd.md claimed 5 in two places while the code used 2, so a worked "
        "example printed a KS_Valid the tool cannot produce",
    ),
    # -- MDS ----------------------------------------------------------------
    Claim(
        "mds.entropy_normaliser",
        "8.0",
        RustPattern(
            "motif_utils.rs",
            r"entropy\s*/\s*([\d.]+)\s*\n\}",
            "MDS entropy normaliser",
        ),
        "log2(256): divides raw 4-mer entropy so MDS lands in [0, 1]. "
        "region-mds.md documented the raw 6.0-8.0 bit scale in its formula "
        "section while its own clinical table quoted 0.95-1.0 -- the page "
        "disagreed with itself, and a threshold built from the former is out "
        "by a factor of 8",
    ),
    # -- UXM ----------------------------------------------------------------
    # These were 0.5/0.5 at the call site for a year, which made the X class
    # unreachable because the backend tests >= methy first.
    Claim(
        "uxm.methylated_threshold",
        "0.75",
        PyConst("krewlyzer.uxm", "METHY_THRESHOLD"),
        "above this a fragment is M; equal thresholds make X unreachable",
    ),
    Claim(
        "uxm.unmethylated_threshold",
        "0.25",
        PyConst("krewlyzer.uxm", "UNMETHY_THRESHOLD"),
        "below this a fragment is U; the gap between the two is the X class",
    ),
    # -- FSC channels -------------------------------------------------------
    # Non-overlapping by design: overlapping channels would make the ML
    # features collinear.
    #
    # There is one entry here because there is now one definition. The genome
    # and gene paths used to carry separate inline comparison chains and had
    # drifted apart; both now call `SizeChannel::of`. If a second set of bounds
    # reappears anywhere, this claim keeps matching the shared one and says
    # nothing -- which is why `fsc_gene_ratios_sum_to_one` checks the emitted
    # values too.
    Claim(
        "fsc.channel_upper_bounds",
        "100,149,220,260,400,1000",
        RustPattern(
            "fsc.rs",
            r"\d+\.\.=(\d+)\s*=>\s*Some\(Self::",
            "SizeChannel::of upper bounds",
        ),
        "the six documented FSC channels, shared by genome bins and genes",
    ),
    # -- PON --------------------------------------------------------------
    # `docs/cli/index.md` quotes this in the validate-pon table, and the
    # builder and the gate enforce it in different files -- so the three can
    # drift unless something ties them together.
    Claim(
        "pon.min_samples_per_key",
        "3",
        PyConst("krewlyzer.pon.build", "MIN_SAMPLES_PER_KEY"),
        "an entry backed by fewer samples is an anecdote, and its sigma is "
        "the difference between two numbers",
    ),
    Claim(
        "pon.gate.min_samples",
        "3",
        PyConst("krewlyzer.validate.pon_gate", "MIN_SAMPLES"),
        "the gate must enforce the same floor the builder applies, or one of "
        "them is decoration",
    ),
)


# ---------------------------------------------------------------------------
# Disagreements that exist today
# ---------------------------------------------------------------------------

KNOWN_DIVERGENCES: Tuple[Divergence, ...] = (
    Divergence(
        "gc.length_bin_ceiling",
        left=RustConst("gc_correction.rs", "NUM_BINS"),
        left_value="188",
        right=RustPattern(
            "extract_motif.rs", r"/\s*5\)\.min\((\d+)\)", "GC length-bin cap"
        ),
        right_value="67",
        why=(
            "extract_motif.rs caps the GC-observation length bin at 67 (400bp) "
            "with a comment claiming it matches gc_correction.rs, which "
            "declares 188 bins spanning 60-999bp. Every fragment of 400bp or "
            "more therefore collapses into a single bin when GC observations "
            "are accumulated, while the correction model expects the full "
            "range. The ultra_long channel (401-1000bp) was added deliberately "
            "for necrosis detection, so this silently under-resolves exactly "
            "the fragments it was added to measure."
        ),
    ),
)

# Resolved, kept as a note rather than an entry:
#
# `fsc.channel_bounds` recorded that the gene/region aggregator split at
# <100/<150/<260/<400 while the genome-bin counter split at
# <=100/<=149/<=220/<=260/<=400 plus an ultra_long channel, so `mono_nucl`
# meant 150-220bp in FSC.parquet and 150-259bp in FSC.gene.parquet. Both paths
# now call `SizeChannel::of`, pinned by `fsc.channel_upper_bounds` above, and
# the divergence cannot recur without deleting that call.
