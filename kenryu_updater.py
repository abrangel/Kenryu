"""KENRYU self-update checker.

Answers the question every local-first tool gets asked: "if the data lives
on my machine, when does it update?"

This module records which version of each reference database KENRYU
currently uses, checks the authoritative source for a newer release, and
reports what is outdated. It is the watcher that keeps the local data
honest: it fetches the authoritative release pages, compares versions, and
tells the user what to download. Downloading and installing a new release
is a separate, explicit step (never silent, never automatic), because
swapping a reference database changes results.

Rules:
- Every check fails soft. No network means status "unknown", never an
  exception and never a false "up to date".
- Versions are compared as dotted numbers ("22" < "23.1").
- A source whose local version is not pinned (miRTarBase, fetched live
  from Harmonizome) reports "unpinned" instead of pretending to be
  current.
- Live APIs (Enrichr) report "live": they are current whenever reachable.
"""

import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import httpx

BASE_DIR = Path(__file__).resolve().parent
VERSIONS_FILE = BASE_DIR / "local_db" / "data_versions.json"
_SYSTEM_CA = Path("/etc/ssl/certs/ca-certificates.crt")

# What KENRYU uses today, and where the truth lives. "latest_known" is the
# newest release verified by a human (see data_versions.json for the notes);
# the live check re-verifies it against the authoritative page.
DATA_SOURCES = [
    {
        "key": "mirbase",
        "name": "miRBase",
        "role": "miRNA names and annotation (sanitization)",
        "local": "22",
        "latest_known": "23.1",
        "check_url": "https://www.mirbase.org/download/CURRENT/README/",
        "parser": "mirbase_readme",
        "download_url": "https://www.mirbase.org/download/",
        "note": "KENRYU validates miRNA names against miRBase v22 (thesis). "
                "v23 (Aug 2026) added 8,847 hairpins; v23.1 (Sep 2026) "
                "renamed/deleted entries.",
    },
    {
        "key": "targetscan",
        "name": "TargetScanHuman",
        "role": "predicted miRNA targets (local database)",
        "local": "8.0",
        "latest_known": "8.0",
        "check_url": "https://www.targetscan.org/vert_80/",
        "parser": "targetscan_page",
        "download_url": "https://www.targetscan.org/vert_80/",
        "note": "Local database (targetscan_full.json.zip, 8.3 MB). "
                "No v9.0 release found as of 2026-11-18.",
    },
    {
        "key": "mirtarbase",
        "name": "miRTarBase",
        "role": "experimentally validated miRNA targets",
        "local": "live via Harmonizome (unpinned)",
        "latest_known": "10.0",
        "check_url": "https://mirtarbase.cuhk.edu.cn/php/download.php",
        "parser": "mirtarbase_page",
        "download_url": "https://mirtarbase.cuhk.edu.cn/php/download.php",
        "note": "Validated targets are fetched live from Harmonizome; the "
                "underlying miRTarBase version is not pinned. miRTarBase "
                "10.0 (Sep 2024) holds 3,817,550 validated interactions "
                "from 13,690 articles.",
    },
    {
        "key": "enrichr",
        "name": "Enrichr",
        "role": "functional enrichment (KEGG/Reactome/WikiPathways/GO)",
        "local": "live API",
        "latest_known": "live",
        "check_url": "https://maayanlab.cloud/Enrichr/datasetStatistics",
        "parser": "enrichr_stats",
        "download_url": "https://maayanlab.cloud/Enrichr/",
        "note": "Enrichment is queried live; libraries are always current "
                "when online. No local copy yet (see local-enrichment work).",
    },
]


def load_local_versions(path=VERSIONS_FILE):
    """Pinned local versions from data_versions.json (falls back to the
    registry defaults above when the file is missing). Always returns a
    flat {source_key: local_version} dict, whatever the file's shape."""
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return {
            k: (v.get("local") if isinstance(v, dict) else v)
            for k, v in data.items()
        }
    except (OSError, ValueError):
        return {s["key"]: s["local"] for s in DATA_SOURCES}


def version_tuple(version):
    """'23.1' -> (23, 1); '22' -> (22,). Returns None when unparseable."""
    try:
        return tuple(int(p) for p in str(version).strip().split("."))
    except ValueError:
        return None


def compare_versions(local, latest):
    """-1 if local is older, 0 if equal, 1 if local is newer, None if the
    comparison cannot be made (unpinned / live / unparseable)."""
    lv, rv = version_tuple(local), version_tuple(latest)
    if lv is None or rv is None:
        return None
    return (lv > rv) - (lv < rv)


def parse_mirbase_readme(text):
    """First line of the miRBase CURRENT README:
    'The miRBase Sequence Database -- Release 23.1 (Sep 2026)'."""
    m = re.search(r"Release\s+(\d+(?:\.\d+)*)", text or "")
    return m.group(1) if m else None


def parse_targetscan_page(text):
    """TargetScanHuman 8.0 on https://www.targetscan.org/vert_80/."""
    m = re.search(r"TargetScanHuman\s+(\d+(?:\.\d+)*)", text or "")
    return m.group(1) if m else None


def parse_mirtarbase_page(text):
    """'miRTarBase 10.0' on the official download page."""
    m = re.search(r"miRTarBase\s+(\d+(?:\.\d+)*)", text or "")
    return m.group(1) if m else None


def parse_enrichr_stats(text):
    """datasetStatistics returns JSON with a 'statistics' list; any valid
    payload means the live API is reachable and current."""
    try:
        data = json.loads(text or "")
    except ValueError:
        return None
    return "live" if isinstance(data, dict) and "statistics" in data else None


_PARSERS = {
    "mirbase_readme": parse_mirbase_readme,
    "targetscan_page": parse_targetscan_page,
    "mirtarbase_page": parse_mirtarbase_page,
    "enrichr_stats": parse_enrichr_stats,
}


def check_source(entry, client):
    """Check one source. Never raises: any failure becomes 'unknown'."""
    result = {
        "key": entry["key"],
        "name": entry["name"],
        "role": entry["role"],
        "local": entry["local"],
        "latest_known": entry["latest_known"],
        "latest_checked": None,
        "status": "unknown",
        "check_url": entry["check_url"],
        "download_url": entry["download_url"],
        "note": entry["note"],
    }
    # No claim is made before the fetch succeeds: an unreachable source
    # stays "unknown", it is never reported as current.
    try:
        resp = client.get(entry["check_url"], follow_redirects=True)
        resp.raise_for_status()
    except Exception:
        return result
    latest = _PARSERS[entry["parser"]](resp.text)
    if latest is None:
        return result
    result["latest_checked"] = latest
    if entry["latest_known"] == "live":
        # A live API is current whenever it answers.
        result["status"] = "live"
        return result
    cmp = compare_versions(entry["local"], latest)
    if cmp is None:
        # Local version not pinned (e.g. miRTarBase via Harmonizome).
        result["status"] = "unpinned"
    elif cmp < 0:
        result["status"] = "outdated"
    else:
        result["status"] = "up-to-date"
    return result


def _make_client(timeout):
    """Build an httpx client that actually reaches the internet from
    wherever KENRYU runs.

    Two traps, both learned the hard way on this VM:
    - httpx 0.26 with trust_env=True raises InvalidURL when NO_PROXY holds
      bracketed IPv6 (``[::1]``). Passing the proxy explicitly with
      trust_env=False skips that parsing entirely (our four hosts are all
      public, so no_proxy exceptions are irrelevant here).
    - Egress proxies that intercept TLS need their CA: honour SSL_CERT_FILE
      when set, else the system bundle, else the default bundle.
    """
    proxy = (
        os.environ.get("HTTPS_PROXY")
        or os.environ.get("https_proxy")
        or os.environ.get("HTTP_PROXY")
        or os.environ.get("http_proxy")
    )
    verify = os.environ.get("SSL_CERT_FILE") or (
        str(_SYSTEM_CA) if _SYSTEM_CA.exists() else True
    )
    if proxy:
        return httpx.Client(proxy=proxy, trust_env=False, verify=verify, timeout=timeout)
    return httpx.Client(trust_env=True, verify=verify, timeout=timeout)


def check_updates(client=None, timeout=15.0):
    """Check every source. Returns a report dict; never raises. When no
    client is given, one is created (see _make_client) and closed here."""
    own_client = client is None
    if own_client:
        client = _make_client(timeout)
    try:
        sources = [check_source(e, client) for e in DATA_SOURCES]
    finally:
        if own_client:
            client.close()
    return {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "updates_available": sum(1 for s in sources if s["status"] == "outdated"),
        "unpinned": sum(1 for s in sources if s["status"] == "unpinned"),
        "unknown": sum(1 for s in sources if s["status"] == "unknown"),
        "sources": sources,
    }


def summarize(report):
    """One honest line per source, for humans."""
    lines = []
    for s in report["sources"]:
        if s["status"] == "outdated":
            lines.append(
                f"{s['name']}: local {s['local']} -> {s['latest_checked']} "
                f"available ({s['download_url']})"
            )
        elif s["status"] == "up-to-date":
            lines.append(f"{s['name']}: {s['local']} up to date")
        elif s["status"] == "live":
            lines.append(f"{s['name']}: live API, current when online")
        elif s["status"] == "unpinned":
            lines.append(
                f"{s['name']}: version not pinned "
                f"(latest known {s['latest_known']})"
            )
        else:
            lines.append(f"{s['name']}: could not check (offline?)")
    return lines
