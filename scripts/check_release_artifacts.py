#!/usr/bin/env python3
"""Verify that the newest tag was actually published everywhere it should be.

Every other gate in this repository checks *intent* -- that the code is right,
that a workflow contains a deploy step, that the docs are current in git. None
of them checks *reality*: whether the thing a release is supposed to produce is
actually out there.

That gap is not hypothetical. 0.9.0 was tagged, PyPI and the container
published, and:

* the documentation site kept serving 0.8.3 -- the release whose defects 0.9.0
  exists to fix -- because the deploy had been deleted from CI three weeks
  earlier and nothing noticed for a year;
* the GitHub Release was never created at all, because the workflow had never
  created one and the Releases page was maintained by hand.

Both were green everywhere. A config-level test cannot catch either: the first
was a deleted step, the second a step that never existed, and neither would be
caught by a check that a *different* step still runs. The only thing that
distinguishes "published" from "looks published" is asking the registries.

Exit codes
----------
0   every artifact matches the newest tag
1   at least one artifact is missing or stale -- a real finding
2   at least one source could not be reached, and nothing was found stale

Two is separate on purpose. A network failure reported as "in sync" is the
silent pass this whole file exists to prevent, and a network failure reported
as "out of sync" sends someone hunting for a release problem that is not there.
A missing scanner is not a clean scan.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from typing import Any, Callable, Optional

REPO = "msk-access/krewlyzer"
PACKAGE = "krewlyzer"
TIMEOUT = 30


class Unreachable(RuntimeError):
    """The source could not be consulted. Distinct from 'it disagreed'."""


def _get_json(url: str, headers: Optional[dict] = None) -> Any:
    req = urllib.request.Request(url, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            return json.load(resp)
    except urllib.error.HTTPError as exc:
        raise Unreachable(f"{url} -> HTTP {exc.code}") from exc
    except Exception as exc:  # noqa: BLE001 - any failure to reach is the same class
        raise Unreachable(f"{url} -> {exc}") from exc


# ---------------------------------------------------------------------------
# What "the newest release" is
# ---------------------------------------------------------------------------


def newest_tag() -> str:
    """The highest release tag, from git if available, else the GitHub API.

    Sorted with `version:refname` rather than by commit date: a tag can be
    pushed out of order, and 0.10.0 must beat 0.9.0 rather than losing a string
    comparison.
    """
    try:
        out = subprocess.run(
            ["git", "tag", "--list", "--sort=-version:refname"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.split()
        tags = [t for t in out if t and not t.startswith("v")]
        if tags:
            return tags[0]
    except Exception:  # noqa: BLE001 - fall through to the API
        pass

    data = _get_json(f"https://api.github.com/repos/{REPO}/tags?per_page=100")
    names = [t["name"] for t in data if not t["name"].startswith("v")]
    if not names:
        raise Unreachable("no release tags found in git or on GitHub")
    return sorted(names, key=lambda s: [int(p) for p in s.split(".")], reverse=True)[0]


# ---------------------------------------------------------------------------
# The four places a release is supposed to land
# ---------------------------------------------------------------------------


def pypi_latest() -> str:
    data = _get_json(f"https://pypi.org/pypi/{PACKAGE}/json")
    return data["info"]["version"]


def ghcr_tags() -> list[str]:
    # Anonymous pull scope is enough to list tags on a public package, so this
    # needs no secret and runs identically on a laptop and in CI.
    auth = _get_json(
        f"https://ghcr.io/token?scope=repository:{REPO}:pull&service=ghcr.io"
    )
    token = auth.get("token")
    if not token:
        raise Unreachable("ghcr.io issued no anonymous pull token")
    data = _get_json(
        f"https://ghcr.io/v2/{REPO}/tags/list",
        headers={"Authorization": f"Bearer {token}"},
    )
    return list(data.get("tags") or [])


def docs_stable() -> str:
    """The version the `stable` alias resolves to on the published site.

    Read from the live site rather than the gh-pages branch: the branch having
    the right content and the CDN serving it are different claims, and it is
    the second one readers experience.
    """
    org, repo = REPO.split("/")
    data = _get_json(f"https://{org}.github.io/{repo}/versions.json")
    for entry in data:
        if "stable" in (entry.get("aliases") or []):
            return entry["version"]
    raise Unreachable("no version carries the 'stable' alias")


def github_release() -> str:
    data = _get_json(f"https://api.github.com/repos/{REPO}/releases/latest")
    name = data.get("tag_name")
    if not name:
        raise Unreachable("the latest release carries no tag_name")
    return name


# ---------------------------------------------------------------------------


@dataclass
class Check:
    name: str
    resolve: Callable[[], Any]
    matches: Callable[[Any, str], bool]
    render: Callable[[Any], str]
    remedy: str


CHECKS = (
    Check(
        "PyPI",
        pypi_latest,
        lambda got, tag: got == tag,
        lambda got: str(got),
        "the Release workflow's Publish step failed, or the tag was never pushed",
    ),
    Check(
        "GHCR container",
        ghcr_tags,
        lambda got, tag: tag in got,
        lambda got: ", ".join(sorted(got)[-4:]),
        "the Release workflow's Docker step failed",
    ),
    Check(
        "docs `stable`",
        docs_stable,
        lambda got, tag: got == tag,
        lambda got: str(got),
        "the Docs workflow did not deploy for this tag; publish it with "
        "Actions -> Docs -> Run workflow -> publish_stable",
    ),
    Check(
        "GitHub Release",
        github_release,
        lambda got, tag: got == tag,
        lambda got: str(got),
        "no Release object for this tag; release.yml creates one from the "
        "CHANGELOG section, so check that the section exists",
    ),
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tag",
        help="version to check for (default: the newest release tag)",
    )
    args = ap.parse_args()

    try:
        tag = args.tag or newest_tag()
    except Unreachable as exc:
        print(f"could not determine the newest tag: {exc}", file=sys.stderr)
        return 2

    print(f"Checking published artifacts against {tag}\n")

    stale: list[str] = []
    unreachable: list[str] = []

    for check in CHECKS:
        try:
            got = check.resolve()
        except Unreachable as exc:
            print(f"  ?  {check.name:<16} could not check -- {exc}")
            unreachable.append(check.name)
            continue

        if check.matches(got, tag):
            print(f"  ok {check.name:<16} {check.render(got)}")
        else:
            print(f"  XX {check.name:<16} {check.render(got)}  (expected {tag})")
            stale.append(f"{check.name}: {check.remedy}")

    print()
    if stale:
        print(f"{len(stale)} artifact(s) do not match {tag}:", file=sys.stderr)
        for line in stale:
            print(f"  - {line}", file=sys.stderr)
        if unreachable:
            print(
                f"\nAlso could not check: {', '.join(unreachable)}. "
                "There may be more.",
                file=sys.stderr,
            )
        return 1

    if unreachable:
        print(
            f"Could not check {', '.join(unreachable)}. Nothing was found stale, "
            "but this is NOT a clean result -- rerun before trusting it.",
            file=sys.stderr,
        )
        return 2

    print(f"All published artifacts match {tag}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
