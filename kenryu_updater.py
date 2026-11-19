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
        "files": [
            {"url": "https://www.mirbase.org/download/CURRENT/mature.fa",
             "local": "local_db/mirbase_mature.fa", "verify": "fasta",
             "transform": "html_fasta"},
            {"url": "https://www.mirbase.org/download/CURRENT/hairpin.fa",
             "local": "local_db/mirbase_hairpin.fa", "verify": "fasta",
             "transform": "html_fasta"},
        ],
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
        "files": [],
        "note": "Local database (targetscan_full.json.zip, 8.3 MB). "
                "No v9.0 release found as of 2026-11-18. No direct file URL "
                "registered yet: downloads go through the data-download CGI "
                "(/cgi-bin/targetscan/data_download.vert80.cgi).",
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
        "files": [],
        "note": "Validated targets are fetched live from Harmonizome; the "
                "underlying miRTarBase version is not pinned. miRTarBase "
                "10.0 (Sep 2024) holds 3,817,550 validated interactions "
                "from 13,690 articles. No direct file URL registered yet "
                "(site unreachable from this machine).",
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
        "files": [],
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


def _client_kwargs(timeout):
    """Shared httpx client options: explicit proxy + egress CA when a
    proxy is configured (skips the broken NO_PROXY parsing in httpx 0.26),
    else the environment defaults."""
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
        return {"proxy": proxy, "trust_env": False, "verify": verify, "timeout": timeout}
    return {"trust_env": True, "verify": verify, "timeout": timeout}


def _make_client(timeout):
    """Build an httpx client that actually reaches the internet from
    wherever KENRYU runs. See _client_kwargs for the two traps handled."""
    return httpx.Client(**_client_kwargs(timeout))


def make_async_client(timeout):
    """Async twin of _make_client (the analysis path uses AsyncClient)."""
    return httpx.AsyncClient(**_client_kwargs(timeout))


def check_updates(client=None, timeout=15.0, local_versions=None):
    """Check every source. Returns a report dict; never raises. When no
    client is given, one is created (see _make_client) and closed here.
    The 'local' version compared is the pinned record in
    data_versions.json (or `local_versions` when given), not the registry
    default, so a source installed yesterday is reported current today."""
    own_client = client is None
    if own_client:
        client = _make_client(timeout)
    pinned = local_versions if local_versions is not None else load_local_versions()
    try:
        sources = [
            check_source({**e, "local": pinned.get(e["key"], e["local"])}, client)
            for e in DATA_SOURCES
        ]
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


# ── Download + install ───────────────────────────────────────────────────
# The second half of the updater: when check_updates() says a source is
# outdated, this downloads the new release, verifies it, swaps it in with
# the old file kept as a backup, and records the new version. Nothing is
# ever installed silently: the caller (API or daily job) asks for it by
# source key, and only the registry's own whitelisted file URLs are used.

import hashlib
import html as html_lib
import os
import shutil
import zipfile


def verify_fasta(path):
    """A FASTA file must start with '>' and hold at least one entry."""
    with open(path, encoding="utf-8", errors="strict") as fh:
        first = fh.readline()
        if not first.startswith(">"):
            raise ValueError(f"{path}: not a FASTA file (first line: {first[:40]!r})")
        for line in fh:
            if line.startswith(">"):
                return True
    raise ValueError(f"{path}: FASTA file has no entries")


def verify_zip(path):
    """A ZIP must open and pass its internal CRC check."""
    with zipfile.ZipFile(path) as zf:
        if zf.testzip() is not None:
            raise ValueError(f"{path}: ZIP failed its CRC check")
    return True


_VERIFIERS = {"fasta": verify_fasta, "zip": verify_zip}


def transform_html_fasta(path):
    """miRBase serves its .fa files wrapped in HTML (``<p>...</p>``,
    ``<br>`` line breaks, ``&gt;`` escapes). Unwrap to plain FASTA, in
    place. Raises ValueError if nothing FASTA-like remains."""
    text = Path(path).read_text(encoding="utf-8", errors="strict")
    text = re.sub(r"</?p>", "", text)
    text = text.replace("<br>", "\n").replace("<br/>", "\n").replace("<br />", "\n")
    text = html_lib.unescape(text)
    text = re.sub(r"<[^>]+>", "", text)
    text = text.strip() + "\n"
    if not text.startswith(">"):
        raise ValueError(f"{path}: no FASTA content after HTML unwrap")
    Path(path).write_text(text, encoding="utf-8")
    return True


_TRANSFORMS = {"html_fasta": transform_html_fasta}


def download_file(url, dest, client, timeout=300.0, min_bytes=1, verify=None,
                  transform=None):
    """Stream-download url to dest. Writes to a temp file first and only
    renames it into place after the optional transform, the size check and
    the optional verifier pass, so a failed download never leaves a partial
    file at dest. Returns (dest, sha256_hex, size_bytes). Raises on any
    failure."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    sha = hashlib.sha256()
    size = 0
    try:
        with client.stream("GET", url, follow_redirects=True, timeout=timeout) as resp:
            resp.raise_for_status()
            with open(tmp, "wb") as fh:
                for chunk in resp.iter_bytes(1 << 16):
                    fh.write(chunk)
                    sha.update(chunk)
                    size += len(chunk)
        if size < min_bytes:
            raise ValueError(f"{url}: downloaded {size} bytes, expected >= {min_bytes}")
        if transform:
            _TRANSFORMS[transform](tmp)
            # Hash/size describe the stored file, not the raw download.
            data = tmp.read_bytes()
            sha = hashlib.sha256(data)
            size = len(data)
        if verify:
            _VERIFIERS[verify](tmp)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, dest)
    return dest, sha.hexdigest(), size


def _record_version(key, new_version, versions_path=VERSIONS_FILE):
    """Update the pinned local version in data_versions.json."""
    try:
        with open(versions_path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        data = {}
    entry = data.get(key, {})
    if not isinstance(entry, dict):
        entry = {"local": entry}
    old = entry.get("local")
    entry["local"] = new_version
    entry["updated_at"] = datetime.now(timezone.utc).isoformat()
    data[key] = entry
    versions_path.parent.mkdir(parents=True, exist_ok=True)
    with open(versions_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
    return old


def install_source_update(key, client=None, base_dir=BASE_DIR, timeout=300.0):
    """Download + install the new release for one registered source.

    Steps: verify the latest release exists (check_source), download each
    registered file, verify it, back the old file up as
    ``<name>.bak-<old-version>``, swap it in, record the new version.
    Raises ValueError for an unknown source or one with no downloadable
    files, RuntimeError when the latest version cannot be verified.
    Returns a report dict with the installed files, hashes and backups.
    """
    entry = next((s for s in DATA_SOURCES if s["key"] == key), None)
    if entry is None:
        raise ValueError(f"unknown source: {key}")
    if not entry.get("files"):
        raise ValueError(
            f"{entry['name']}: no downloadable files registered yet "
            f"(see the note in the registry)"
        )
    own_client = client is None
    if own_client:
        client = _make_client(timeout)
    try:
        checked = check_source(entry, client)
        new_version = checked["latest_checked"]
        if new_version is None:
            raise RuntimeError(
                f"{entry['name']}: could not verify the latest release "
                f"({entry['check_url']}); refusing to install blind"
            )
        old_version = load_local_versions(
            base_dir / "local_db" / "data_versions.json").get(key)
        installed, backups = [], []
        for spec in entry["files"]:
            dest = base_dir / spec["local"]
            if dest.exists():
                backup = dest.with_name(
                    f"{dest.stem}.bak-{old_version or 'unknown'}{dest.suffix}"
                )
                shutil.copy2(dest, backup)
                backups.append(str(backup))
            _, sha, size = download_file(
                spec["url"], dest, client,
                timeout=timeout, min_bytes=1, verify=spec.get("verify"),
                transform=spec.get("transform"),
            )
            installed.append({
                "file": spec["local"], "url": spec["url"],
                "sha256": sha, "bytes": size,
            })
        _record_version(
            key, new_version,
            versions_path=base_dir / "local_db" / "data_versions.json")
        return {
            "source": key,
            "name": entry["name"],
            "old_version": old_version,
            "new_version": new_version,
            "installed": installed,
            "backups": backups,
            "notice": (
                f"{entry['name']} updated {old_version} -> {new_version}. "
                f"Previous files kept as backups; re-run your panels to see "
                f"what changed in the results."
            ),
        }
    finally:
        if own_client:
            client.close()
