"""The single-sample HTML report.

For internal use. It contains one patient's actual measurements, so it is
generated on demand and not committed or published; ``describe-output`` covers
the shareable, structural case.

Structure
---------
**Summary** carries everything sample-level — the cross-axis verdict, whether a
PON was supplied, what is present and what is not. **Everything after it is
organised by table**, because the question a reader actually arrives with is
"what is this file and what is in it", not "what is the verdict". A chart lives
inside the section for the table it was computed from; a chart three screens
away from its data is a chart nobody connects to the data.

Each table section answers **why** the measurement exists, **how** the number is
arrived at, and **what** to look at — from ``meaning.py``, so the report and the
contract cannot drift apart.

Theme
-----
Material for MkDocs' palette, so the report and the docs site look like one
product: ``#ef5552`` — "krew" is blood in Polish — on Material's neutral scale,
Roboto and Roboto Mono. Light and dark are token-level, with an explicit toggle
that persists, and the Plotly figures re-theme with it.
"""

from __future__ import annotations

from datetime import datetime
from html import escape
from typing import Dict, List, Optional

import pandas as pd

from krewlyzer.validate.describe import SampleReport, TableFacts, human_bytes
from krewlyzer.validate.plots import Chart, charts_by_suffix
from krewlyzer.validate.verdict import Verdict

_CSS = """
:root {
  --accent:      #ef5552;
  --accent-dark: #c62828;
  --accent-soft: #ef55521a;
  --oppose:      #3b7ea1;
  --warn:        #b9770e;

  --fg:          rgba(0,0,0,.87);
  --fg-soft:     rgba(0,0,0,.62);
  --fg-faint:    rgba(0,0,0,.40);
  --fg-ghost:    rgba(0,0,0,.20);
  --bg:          #fff;
  --bg-sunk:     #f5f5f5;
  --rule:        rgba(0,0,0,.10);
  --shadow:      0 1px 2px rgba(0,0,0,.06), 0 4px 16px rgba(0,0,0,.05);
  --plot-paper:  #fff;
  --plot-grid:   rgba(0,0,0,.09);
  --plot-ink:    rgba(0,0,0,.62);
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --accent: #ff7875; --accent-dark: #ef5552; --accent-soft: #ff78751f;
    --oppose: #6cb2d6; --warn: #e0a458;
    --fg: rgba(255,255,255,.88); --fg-soft: rgba(255,255,255,.62);
    --fg-faint: rgba(255,255,255,.40); --fg-ghost: rgba(255,255,255,.18);
    --bg: #1f2129; --bg-sunk: #191b22; --rule: rgba(255,255,255,.12);
    --shadow: none;
    --plot-paper: #1f2129; --plot-grid: rgba(255,255,255,.10);
    --plot-ink: rgba(255,255,255,.62);
  }
}
:root[data-theme="dark"] {
  --accent: #ff7875; --accent-dark: #ef5552; --accent-soft: #ff78751f;
  --oppose: #6cb2d6; --warn: #e0a458;
  --fg: rgba(255,255,255,.88); --fg-soft: rgba(255,255,255,.62);
  --fg-faint: rgba(255,255,255,.40); --fg-ghost: rgba(255,255,255,.18);
  --bg: #1f2129; --bg-sunk: #191b22; --rule: rgba(255,255,255,.12);
  --shadow: none;
  --plot-paper: #1f2129; --plot-grid: rgba(255,255,255,.10);
  --plot-ink: rgba(255,255,255,.62);
}

*, *::before, *::after { box-sizing: border-box; }
html { scroll-behavior: smooth; }
@media (prefers-reduced-motion: reduce) {
  html { scroll-behavior: auto; }
  * { animation: none !important; transition: none !important; }
}

body {
  margin: 0; background: var(--bg); color: var(--fg);
  font: 400 15px/1.65 Roboto, -apple-system, "Segoe UI", Helvetica, sans-serif;
  -webkit-font-smoothing: antialiased;
}
code, .mono, td.num, th { font-family: "Roboto Mono", ui-monospace, SFMono-Regular, Menlo, monospace; }

.shell { display: grid; grid-template-columns: 17rem minmax(0,1fr); }
@media (max-width: 64rem) { .shell { grid-template-columns: 1fr; } }

/* ---------- navigation ---------- */
.rail {
  position: sticky; top: 0; height: 100vh; overflow-y: auto;
  border-right: 1px solid var(--rule); background: var(--bg-sunk);
  padding: 1.5rem 1rem 3rem 1.5rem;
}
@media (max-width: 64rem) {
  .rail { position: static; height: auto; border-right: 0; border-bottom: 1px solid var(--rule); }
}
.rail .brand { display: flex; align-items: baseline; gap: .5rem; margin-bottom: .1rem; }
.rail .brand b { color: var(--accent); font: 700 15px/1 Roboto, sans-serif; letter-spacing: -.01em; }
.rail .brand span { font: 400 10px/1 "Roboto Mono", monospace; color: var(--fg-faint); letter-spacing: .08em; }
.rail .sample { font: 500 12.5px/1.45 "Roboto Mono", monospace; color: var(--fg-soft); word-break: break-all; margin-bottom: 1.25rem; }
.rail h2 {
  font: 500 10.5px/1.4 "Roboto Mono", monospace; letter-spacing: .12em; text-transform: uppercase;
  color: var(--fg-faint); margin: 1.4rem 0 .4rem;
}
.rail ul { list-style: none; margin: 0; padding: 0; }
.rail li { margin: 0; }
.rail a {
  display: block; padding: .26rem .6rem; border-radius: 4px;
  color: var(--fg-soft); text-decoration: none; font-size: 13.5px;
  border-left: 2px solid transparent; margin-left: -.6rem;
}
.rail a:hover { background: var(--accent-soft); color: var(--fg); }
.rail a:focus-visible { outline: 2px solid var(--accent); outline-offset: -2px; }
.rail a.on { color: var(--accent); border-left-color: var(--accent); background: var(--accent-soft); }
.rail a.dim { color: var(--fg-ghost); }
.rail a code { font-size: 12.5px; }

.themer {
  margin-top: 2rem; display: flex; gap: .35rem; align-items: center;
  font: 500 10.5px/1 "Roboto Mono", monospace; letter-spacing: .1em;
  text-transform: uppercase; color: var(--fg-faint);
}
.themer button {
  font: inherit; color: var(--fg-soft); background: none; cursor: pointer;
  border: 1px solid var(--rule); border-radius: 4px; padding: .3rem .5rem;
}
.themer button[aria-pressed="true"] { color: var(--accent); border-color: var(--accent); }
.themer button:focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }

/* ---------- main ---------- */
main { padding: 2.5rem 3rem 7rem; max-width: 64rem; }
@media (max-width: 48rem) { main { padding: 1.75rem 1.15rem 4rem; } }

h1 { font: 400 2rem/1.2 Roboto, sans-serif; letter-spacing: -.01em; margin: 0 0 .35rem; text-wrap: balance; }
.sub { color: var(--fg-soft); margin: 0 0 2rem; font-size: .92rem; }
.sub code { font-size: 12.5px; }

.act {
  font: 500 10.5px/1 "Roboto Mono", monospace; letter-spacing: .14em;
  text-transform: uppercase; color: var(--accent); margin: 3rem 0 .6rem;
}
section { scroll-margin-top: 1rem; margin-bottom: 2.5rem; }
section > h2 { font: 400 1.4rem/1.25 Roboto, sans-serif; margin: 0 0 .15rem; letter-spacing: -.01em; }
section > h2 code { font-size: 1.15rem; font-weight: 500; }
.meta { color: var(--fg-faint); font: 400 12.5px/1.5 "Roboto Mono", monospace; margin: 0 0 1rem; }
.lede { color: var(--fg-soft); margin: 0 0 1rem; max-width: 66ch; }

/* ---------- why / how / what ---------- */
.whw { display: grid; gap: .1rem; margin: 0 0 1.25rem; max-width: 74ch; }
.whw div { display: grid; grid-template-columns: 3.6rem 1fr; gap: .9rem; padding: .45rem 0; border-top: 1px solid var(--rule); }
.whw dt {
  font: 500 10px/1.7 "Roboto Mono", monospace; letter-spacing: .1em;
  text-transform: uppercase; color: var(--accent); margin: 0;
}
.whw dd { margin: 0; color: var(--fg-soft); font-size: .93rem; }

/* ---------- callouts ---------- */
.note { border-left: 3px solid var(--fg-ghost); padding: .6rem 0 .6rem 1rem; margin: 1rem 0; color: var(--fg-soft); font-size: .9rem; max-width: 70ch; }
.note.warn { border-left-color: var(--warn); }
.note.stop { border-left-color: var(--accent); }
.note strong { color: var(--fg); font-weight: 500; }

/* ---------- run facts ---------- */
.facts {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(9rem, 1fr));
  gap: 0; margin: 0 0 1.25rem; border: 1px solid var(--rule); border-radius: 6px;
  overflow: hidden;
}
.facts > div { padding: .7rem .9rem; border-right: 1px solid var(--rule); background: var(--bg-sunk); }
.facts > div:last-child { border-right: 0; }
.facts dt {
  font: 500 9.5px/1.6 "Roboto Mono", monospace; letter-spacing: .1em;
  text-transform: uppercase; color: var(--fg-faint); margin: 0;
}
.facts dd { margin: 0; font: 500 13.5px/1.45 "Roboto Mono", monospace; color: var(--fg); }

/* ---------- verdict ---------- */
.card { background: var(--bg); border: 1px solid var(--rule); border-radius: 6px; box-shadow: var(--shadow); padding: 1.5rem 1.6rem; }
.card .head { font: 400 1.5rem/1.2 Roboto, sans-serif; margin: 0 0 .4rem; letter-spacing: -.01em; }
.card .say { color: var(--fg-soft); margin: 0; max-width: 64ch; }
.axes { display: grid; gap: .4rem; margin-top: 1.4rem; }
.ax { display: grid; grid-template-columns: 3px 1fr auto; gap: .9rem; align-items: center; background: var(--bg-sunk); border-radius: 4px; padding: .65rem .9rem .65rem 0; }
.ax i { display: block; width: 3px; height: 100%; min-height: 2.4rem; border-radius: 2px; background: var(--fg-ghost); }
.ax.agrees i { background: var(--accent); }
.ax.opposes i { background: var(--oppose); }
.ax.quiet i { background: var(--fg-ghost); }
.ax .nm { font-weight: 500; }
.ax .dt { display: block; color: var(--fg-faint); font-size: .8rem; margin-top: .12rem; }
.ax .dt code { font-size: 11px; }
.ax .vl { text-align: right; font: 500 1.05rem/1 "Roboto Mono", monospace; font-variant-numeric: tabular-nums; }
.ax.agrees .vl { color: var(--accent); }
.ax.opposes .vl { color: var(--oppose); }
.ax.notassessable .vl { font-size: .78rem; font-weight: 400; color: var(--fg-faint); max-width: 14rem; }
.ax .st { display: block; font: 500 9px/1.6 "Roboto Mono", monospace; letter-spacing: .1em; text-transform: uppercase; color: var(--fg-faint); }
.ax.agrees .st { color: var(--accent); }
.ax.opposes .st { color: var(--oppose); }

/* ---------- figures ---------- */
figure { margin: 1.25rem 0; }
figure .figtitle {
  font: 500 1rem/1.35 Roboto, sans-serif; margin: 0 0 .35rem; color: var(--fg);
}
figure .cap { color: var(--fg-soft); font-size: .89rem; margin: .5rem 0 0; max-width: 72ch; }
.plot { border: 1px solid var(--rule); border-radius: 6px; overflow: hidden; }
.nodraw { border: 1px dashed var(--rule); border-radius: 6px; padding: .85rem 1rem; color: var(--fg-faint); font-size: .89rem; background: var(--bg-sunk); }
.nodraw strong { color: var(--fg-soft); font-weight: 500; }

/* ---------- tables ---------- */
.tw { overflow-x: auto; border: 1px solid var(--rule); border-radius: 6px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { text-align: left; padding: .4rem .75rem; white-space: nowrap; border-bottom: 1px solid var(--rule); }
tr:last-child td { border-bottom: 0; }
th { font: 500 9.5px/2 "Roboto Mono", monospace; letter-spacing: .1em; text-transform: uppercase; color: var(--fg-faint); background: var(--bg-sunk); position: sticky; top: 0; }
td.num { text-align: right; font-variant-numeric: tabular-nums; font-size: 12px; }
tbody tr:hover { background: var(--accent-soft); }
td .pon { color: var(--warn); font: 500 9px/1 "Roboto Mono", monospace; letter-spacing: .06em; }

.tag { display: inline-block; font: 500 9px/1.7 "Roboto Mono", monospace; letter-spacing: .07em; text-transform: uppercase; padding: 0 .4rem; border-radius: 3px; border: 1px solid currentColor; }
.tag.read { color: var(--accent); }
.tag.unread { color: var(--fg-faint); }

.grid2 { columns: 2; font: 12.5px/1.9 "Roboto Mono", monospace; color: var(--fg-faint); }
@media (max-width: 48rem) { .grid2 { columns: 1; } }

footer { border-top: 1px solid var(--rule); margin-top: 4rem; padding-top: 1.25rem; color: var(--fg-faint); font-size: .82rem; max-width: 70ch; }
"""

_JS = """
(function () {
  var root = document.documentElement;
  var KEY = 'krewlyzer-report-theme';

  function apply(mode) {
    if (mode) { root.setAttribute('data-theme', mode); }
    else { root.removeAttribute('data-theme'); }
    document.querySelectorAll('.themer button').forEach(function (b) {
      b.setAttribute('aria-pressed', String(b.dataset.mode === (mode || 'auto')));
    });
    restyle();
  }

  function dark() {
    var set = root.getAttribute('data-theme');
    if (set) { return set === 'dark'; }
    return window.matchMedia('(prefers-color-scheme: dark)').matches;
  }

  // Plotly draws to a canvas of its own, so it cannot inherit CSS variables.
  // Re-layout every figure whenever the theme changes.
  function restyle() {
    if (!window.Plotly) { return; }
    var isDark = dark();
    var paper = isDark ? '#1f2129' : '#ffffff';
    var grid  = isDark ? 'rgba(255,255,255,.10)' : 'rgba(0,0,0,.09)';
    var ink   = isDark ? 'rgba(255,255,255,.62)' : 'rgba(0,0,0,.62)';
    document.querySelectorAll('.js-plotly-plot').forEach(function (el) {
      window.Plotly.relayout(el, {
        'paper_bgcolor': paper, 'plot_bgcolor': paper,
        'font.color': ink,
        'xaxis.gridcolor': grid, 'xaxis.linecolor': grid, 'xaxis.tickfont.color': ink,
        'yaxis.gridcolor': grid, 'yaxis.linecolor': grid, 'yaxis.tickfont.color': ink
      });
    });
  }

  document.querySelectorAll('.themer button').forEach(function (b) {
    b.addEventListener('click', function () {
      var mode = b.dataset.mode === 'auto' ? null : b.dataset.mode;
      if (mode) { localStorage.setItem(KEY, mode); } else { localStorage.removeItem(KEY); }
      apply(mode);
    });
  });
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', restyle);
  apply(localStorage.getItem(KEY));

  // Highlight the section currently in view.
  var links = {};
  document.querySelectorAll('.rail a[href^="#"]').forEach(function (a) {
    links[a.getAttribute('href').slice(1)] = a;
  });
  var seen = new Set();
  var obs = new IntersectionObserver(function (entries) {
    entries.forEach(function (e) {
      if (e.isIntersecting) { seen.add(e.target.id); } else { seen.delete(e.target.id); }
    });
    Object.keys(links).forEach(function (id) { links[id].classList.remove('on'); });
    var first = Object.keys(links).find(function (id) { return seen.has(id); });
    if (first) { links[first].classList.add('on'); }
  }, { rootMargin: '-8% 0px -80% 0px' });
  document.querySelectorAll('section[id]').forEach(function (s) { obs.observe(s); });
})();
"""


def _slug(text: str) -> str:
    return "t-" + "".join(c if c.isalnum() else "-" for c in text.lower()).strip("-")


def _axes_html(verdict: Verdict) -> str:
    out = []
    for a in verdict.axes:
        cls = a.support.value.replace(" ", "")
        value = f"{a.value:+.3f}" if a.value is not None else escape(a.reason or "—")
        bits = [escape(a.direction)]
        if a.detail:
            bits.append(escape(a.detail))
        out.append(
            f'<div class="ax {cls}"><i></i><span><span class="nm">{escape(a.name)}</span>'
            f'<span class="dt">{" · ".join(bits)} · <code>{escape(a.source)}</code></span></span>'
            f'<span><span class="vl">{value}</span><span class="st">{escape(a.support.value)}</span></span></div>'
        )
    return "".join(out)


def _figure_html(chart: Chart, index: int, with_runtime: bool = False) -> str:
    # A drawn chart needs its title too. An earlier version rendered the title
    # only on the not-drawn branch, so every figure that worked arrived without
    # a heading -- the caption alone left the reader to infer what they were
    # looking at.
    head = f'<h3 class="figtitle">{escape(chart.title)}</h3>'
    cap = f'<figcaption class="cap">{escape(chart.caption)}</figcaption>'
    if not chart.drawn:
        return (
            f'<figure>{head}<div class="nodraw"><strong>Not drawn.</strong> '
            f"{escape(chart.reason or '')}</div>{cap}</figure>"
        )
    import plotly.io as pio

    div = pio.to_html(
        chart.figure,
        include_plotlyjs="inline" if with_runtime else False,
        full_html=False,
        div_id=f"plot-{index}",
        config={"displaylogo": False, "responsive": True, "scrollZoom": True},
    )
    return f'<figure>{head}<div class="plot">{div}</div>{cap}</figure>'


def _whw(meaning) -> str:
    rows = []
    for label, text in (
        ("Why", meaning.why),
        ("How", meaning.how),
        ("What", meaning.what),
    ):
        if text:
            rows.append(f"<div><dt>{label}</dt><dd>{escape(text)}</dd></div>")
    return f'<dl class="whw">{"".join(rows)}</dl>' if rows else ""


def _table_section(
    t: TableFacts,
    charts: List[Chart],
    start_index: int,
    has_pon: bool,
    runtime_pending: List[bool],
) -> str:
    pon_cols = set(t.meaning.pon_columns) if t.meaning else set()
    meta = [
        f"{t.n_rows:,} rows × {t.n_cols} columns",
        human_bytes(t.size_bytes),
    ]
    # Deliberately not "read downstream". Whether kreview reads a table is a
    # fact about another repository that nothing here keeps in sync, and a
    # label that quietly goes stale is worse than none. This says what this
    # repo actually does: a gated table failing the contract check fails the
    # run; an inventory one is reported and never gates.
    tag = (
        '<span class="tag read">gated</span>'
        if t.consumed
        else '<span class="tag unread">inventory</span>'
    )

    parts = [
        f'<section id="{_slug(t.family)}">',
        f"<h2><code>{escape(t.family)}</code></h2>",
        f'<p class="meta">{" · ".join(meta)} &nbsp; {tag}</p>',
    ]
    if t.meaning:
        parts.append(f'<p class="lede">{escape(t.meaning.measures)}</p>')
        parts.append(_whw(t.meaning))
        if t.meaning.cancer_direction:
            parts.append(
                f'<p class="note"><strong>Cancer direction:</strong> '
                f"{escape(t.meaning.cancer_direction)}</p>"
            )
        if t.meaning.caveat:
            parts.append(f'<p class="note warn">{escape(t.meaning.caveat)}</p>')
    if pon_cols and not has_pon:
        parts.append(
            '<p class="note stop"><strong>No panel of normals was used.</strong> '
            f"{escape(', '.join(sorted(pon_cols)))} "
            + ("is" if len(pon_cols) == 1 else "are")
            + " absent or uninterpretable without one — the raw values below "
            "cannot be compared against a healthy baseline.</p>"
        )

    for offset, chart in enumerate(charts):
        first = runtime_pending[0] and chart.drawn
        if first:
            runtime_pending[0] = False
        parts.append(_figure_html(chart, start_index + offset, with_runtime=first))

    head = (
        "<tr><th>Column</th><th>Type</th><th>Range / example</th>"
        "<th>Distinct</th><th>Null</th></tr>"
    )
    body = []
    for c in t.columns:
        if c.minimum is not None and c.maximum is not None:
            value = f"{c.minimum:.6g} … {c.maximum:.6g}"
        else:
            value = f"<code>{escape(c.example)}</code>"
        flag = ' <span class="pon">PON</span>' if c.name in pon_cols else ""
        distinct = "—" if c.n_distinct < 0 else f"{c.n_distinct:,}"
        body.append(
            f"<tr><td><code>{escape(c.name)}</code>{flag}</td><td>{escape(c.dtype)}</td>"
            f'<td class="num">{value}</td><td class="num">{distinct}</td>'
            f'<td class="num">{c.n_null:,}</td></tr>'
        )
    parts.append(
        f'<div class="tw"><table><thead>{head}</thead><tbody>{"".join(body)}</tbody></table></div>'
    )
    parts.append("</section>")
    return "\n".join(parts)


def run_facts(report: SampleReport, tables: Dict) -> List[tuple]:
    """The run's own configuration, read back out of its output.

    Every number in this report is conditioned on these, and none of them is
    otherwise visible: whether a PON was used decides if z-scores exist at all,
    and panel versus WGS decides which tables can exist.
    """
    facts: List[tuple] = []
    meta = tables.get(".metadata.parquet")

    def cell(name: str):
        if meta is None or meta.empty or name not in meta.columns:
            return None
        value = meta[name].iloc[0]
        return None if pd.isna(value) else value

    total = cell("total_fragments")
    if total is not None:
        facts.append(("Fragments", f"{int(total):,}"))

    for label, column in (("Genome", "genome"), ("Assay", "assay")):
        value = cell(column)
        if value is not None:
            facts.append((label, str(value)))

    panel = cell("panel_mode")
    if panel is not None:
        facts.append(("Mode", "panel" if bool(panel) else "WGS"))
    elif any(t.family.startswith("FSC.gene") for t in report.present):
        # Gene-level output only exists in panel mode, so its presence answers
        # the question when the metadata does not.
        facts.append(("Mode", "panel (inferred from gene-level output)"))

    facts.append(("Panel of normals", "used" if detect_pon(report) else "none"))

    optional = {
        "mFSD": ".mFSD.parquet",
        "UXM": ".UXM.parquet",
        "OCF": ".OCF.offtarget.parquet",
    }
    have = [k for k, v in optional.items() if tables.get(v) is not None]
    facts.append(
        (
            "Optional inputs",
            ", ".join(have) if have else "none — no variants, no methylation BAM",
        )
    )
    return facts


def detect_pon(report: SampleReport) -> bool:
    """Was a panel of normals used?

    Inferred from the output rather than asked for: a PON-derived column is
    present only when one was supplied. Getting this wrong in the safe
    direction — reporting "no PON" when there was one — would add a warning
    nobody needs; the other way round would let a reader take a raw value for a
    normalised one.
    """
    for table in report.present:
        declared = set(table.meaning.pon_columns) if table.meaning else set()
        if declared & {c.name for c in table.columns}:
            return True
    return False


def render_html(
    report: SampleReport,
    verdict: Verdict,
    charts: List[Chart],
    tables: Optional[Dict] = None,
) -> str:
    present = report.present
    tables = tables or {}
    grouped: Dict[str, List[Chart]] = charts_by_suffix(charts)
    has_pon = detect_pon(report)
    drawn = sum(1 for c in charts if c.drawn)
    facts_html = "".join(
        f"<div><dt>{escape(k)}</dt><dd>{escape(str(v))}</dd></div>"
        for k, v in run_facts(report, tables)
    )

    # Charts whose table is absent still need somewhere to say so.
    orphan = [c for c in charts if c.suffix not in {t.suffix for t in present}]

    nav_tables = []
    index = 0
    sections = []
    # The plotly runtime rides along with whichever figure renders first.
    runtime_pending = [True]
    for t in present:
        mine = grouped.get(t.suffix, [])
        sections.append(_table_section(t, mine, index, has_pon, runtime_pending))
        index += len(mine)
        marker = " ●" if any(c.drawn for c in mine) else ""
        nav_tables.append(
            f'<li><a href="#{_slug(t.family)}"><code>{escape(t.family)}</code>{marker}</a></li>'
        )

    pon_banner = ""
    if not has_pon:
        pon_banner = (
            '<div class="note stop" style="margin:0 0 1.5rem"><strong>No panel of '
            "normals in this run.</strong> Every z-score and PON log-ratio is "
            "absent, which is most of the interpretable surface — three of the "
            "four verdict axes below cannot be assessed. Raw values are still "
            "correct, but they carry library prep and batch effects that a PON "
            "is there to cancel, so they are not comparable between samples. "
            "Re-run with <code>--pon-model</code> for anything comparative.</div>"
        )

    orphan_html = ""
    if orphan:
        orphan_html = (
            '<div class="act">Not produced</div><section id="t-absent">'
            "<h2>Charts without data</h2>"
            '<p class="lede">These need an output this run did not produce. '
            "Absence is information: mFSD needs a variant list, UXM needs a "
            "methylation-aware BAM, OCF is hg19-only, gene-level outputs need "
            "panel mode.</p>"
            + "".join(_figure_html(c, 900 + i) for i, c in enumerate(orphan))
            + "</section>"
        )

    missing_html = ""
    if report.missing:
        items = "".join(f"<div>{escape(t.family)}</div>" for t in report.missing)
        missing_html = (
            '<section id="t-missing"><h2>Absent tables</h2>'
            f'<p class="lede">{len(report.missing)} of {len(report.tables)} '
            "contracted tables were not produced.</p>"
            f'<div class="grid2">{items}</div></section>'
        )

    marker_note = ""
    if not report.has_completion_marker:
        marker_note = (
            '<div class="note stop"><strong>No <code>.metadata.parquet</code>.</strong> '
            "The downstream consumer treats it as a completion marker and drops "
            "the sample from the cohort silently — no warning, no error.</div>"
        )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(report.sample_id)} — krewlyzer</title>
<style>{_CSS}</style>
</head>
<body>
<div class="shell">
<div class="rail">
  <div class="brand"><b>krewlyzer</b><span>REPORT</span></div>
  <div class="sample">{escape(report.sample_id)}</div>

  <h2>Summary</h2>
  <ul>
    <li><a href="#t-summary">Cross-axis verdict</a></li>
  </ul>

  <h2>Tables <span style="color:var(--fg-ghost)">● has chart</span></h2>
  <ul>{"".join(nav_tables)}</ul>
  {'<h2>Other</h2><ul>' + ('<li><a href="#t-absent">Charts without data</a></li>' if orphan else '') + ('<li><a href="#t-missing">Absent tables</a></li>' if report.missing else '') + '</ul>' if (orphan or report.missing) else ''}

  <div class="themer">
    <span>Theme</span>
    <button data-mode="auto" aria-pressed="true">Auto</button>
    <button data-mode="light" aria-pressed="false">Light</button>
    <button data-mode="dark" aria-pressed="false">Dark</button>
  </div>
</div>

<main>
  <h1>Fragmentomics output</h1>
  <p class="sub"><code>{escape(str(report.directory))}</code> · {len(present)} of
     {len(report.tables)} tables · {human_bytes(report.total_bytes)} ·
     {drawn} of {len(charts)} charts drawn</p>

  {pon_banner}

  <section id="t-summary">
    <dl class="facts">{facts_html}</dl>
    <div class="card">
      <p class="head">{escape(verdict.headline)}</p>
      <p class="say">{escape(verdict.summary)}</p>
      <div class="axes">{_axes_html(verdict)}</div>
    </div>
    <p class="note">Axes are flagged at |z| &ge; {verdict.z_threshold:g} — a
      conventional value and a tunable one, not a clinical cut-off. An axis
      with no PON z-score reads <em>not assessable</em> rather than counting as
      disagreement, so a thinner run does not look like a healthier one. This
      summarises convergence between measurements; what that means for a
      patient is not something this tool decides.</p>
    {marker_note}
  </section>

  <div class="act">Tables</div>
  {"".join(sections)}

  {orphan_html}
  {missing_html}

  <footer>
    Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} by
    <code>krewlyzer report</code>. Internal use — this document contains one
    sample's measurements. Identifier columns are redacted; nothing else is.
  </footer>
</main>
</div>
<script>{_JS}</script>
</body>
</html>"""
