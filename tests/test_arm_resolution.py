"""Arm resolution (-5p/-3p): a miRNA name with an arm must never silently
receive another arm's TargetScan targets.

Audit (2026-11-21, v1.56, scripts/audit_arm_resolution.py): 33 mature
human miRBase 23.1 names resolve to a TargetScan family whose key carries
no arm; in 14 hairpins (28 names) the -5p and -3p arms return exactly the
same gene set. With different seeds that is impossible in a seed predictor:
at least one of the two results is necessarily wrong (e.g. miR-107,
miR-133b, miR-137, miR-147b, miR-184, miR-217). The fix verifies the arm
by SEED: the TargetScan family seed (miR_Family_Info.txt) must occur in
the query's mature sequence (miRBase 23.1) starting at offset 0-3, which
tolerates 1-2 nt of 5' annotation drift between miRBase versions (e.g.
miR-144-3p, miR-802-5p) while still rejecting the wrong arm, whose seed
is absent from the sequence entirely.

These tests are written fail-first: they fail on v1.57 (the -5p arm
returns the -3p family's targets, there is no arm_unverified field).
"""
import re

import pytest

import kenryu_core as kc
import kenryu_engine as ke

# The 14 hairpins where both arms returned the same set in v1.56.
SWAPPED_HAIRPINS = [
    "107", "133b", "137", "147b", "184", "217", "326", "346",
    "375", "421", "429", "543", "665", "760",
]

LIPID_PANEL = [
    "hsa-miR-33a-5p", "hsa-miR-33b-5p", "hsa-miR-144-3p",
    "hsa-miR-758-3p", "hsa-miR-106b-5p",
]
LIPID_GENES = ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]


def _needs_arm_data():
    """The real miRBase FASTA, family info and TargetScan data are
    tracked files, so they exist in CI and in a fresh clone. Skip
    (rather than fail) if a deployment deleted them."""
    from pathlib import Path
    missing = [p for p in (
        "local_db/mirbase_mature.fa", "local_db/mir_family_info.txt",
        "data/targetscan_full.txt.zip",
    ) if not Path(p).exists()]
    return pytest.mark.skipif(bool(missing), reason=f"arm data missing: {missing}")


# ── Pure core functions (tiny fixtures, run in CI) ──────────────────────────

def test_mirna_seed_is_positions_2_to_8():
    assert kc.mirna_seed("AGCAGCAUUGUACAGGGCUAUCA") == "GCAGCAU"
    assert kc.mirna_seed("AGCUUCUUUACAGUGUUGCCUUGU") == "GCUUCUU"


def test_parse_mirbase_sequences_keeps_human_only(tmp_path):
    fa = tmp_path / "mature.fa"
    fa.write_text(
        ">hsa-miR-107-3p MIMAT0000104 Homo sapiens miR-107-3p\n"
        "AGCAGCAUUGUACAGGGCUAUCA\n"
        ">mmu-miR-107-3p MIMAT0000105 Mus musculus miR-107-3p\n"
        "AGCAGCAUUGUACAGGGCUAUCA\n"
        ">hsa-miR-107-5p MIMAT0050677 Homo sapiens miR-107-5p\n"
        "AGCUUCUUUACAGUGUUGCCUUGU\n"
    )
    seqs = kc.parse_mirbase_sequences(fa)
    assert seqs["mir-107-3p"] == "AGCAGCAUUGUACAGGGCUAUCA"
    assert seqs["mir-107-5p"] == "AGCUUCUUUACAGUGUUGCCUUGU"
    assert "mmu-miR-107-3p" not in seqs
    assert len(seqs) == 2


def test_parse_family_seeds_human_only_and_no_arm_keys(tmp_path):
    fi = tmp_path / "miR_Family_Info.txt"
    fi.write_text(
        "miR family\tSeed+m8\tSpecies ID\tMiRBase ID\tMature sequence\t"
        "Family Conservation?\tMiRBase Accession\n"
        "miR-103-3p/107\tGCAGCAU\t9606\thsa-miR-103a-3p\t"
        "AGCAGCAUUGUACAGGGCUAUGA\t2\tMIMAT0000101\n"
        "miR-103-3p/107\tGCAGCAU\t9606\thsa-miR-107\t"
        "AGCAGCAUUGUACAGGGCUAUCA\t2\tMIMAT0000104\n"
        "miR-103/107\tGCAGCAU\t9606\thsa-miR-103a-3p\t"
        "AGCAGCAUUGUACAGGGCUAUGA\t2\tMIMAT0000101\n"
        "let-7\tUAUACGA\t9615\tcfa-let-7d\tCUAUACGACCUGCUGCCUUUCUUAG\t"
        "-1\tMIMAT0034396\n"
    )
    seeds, conflicts = kc.parse_family_seeds(fi)
    # The map is keyed by family name (full and sub-parts) and by
    # MiRBase ID, all lowercased: the armless sub-part "107" and the
    # MiRBase ID without arm both map to the family seed; the full
    # family name and the armed MiRBase ID map to it too.
    assert seeds["mir-107"] == {"GCAGCAU"}
    assert seeds["mir-103-3p/107"] == {"GCAGCAU"}
    assert seeds["mir-103a-3p"] == {"GCAGCAU"}
    assert "let-7" not in seeds  # not human (Species ID 9615)
    assert conflicts == []


def test_verify_arm_seed_rule():
    # miR-107: the -3p arm carries the family seed at offset 1,
    # the -5p arm does not carry it at all.
    assert kc.verify_arm_seed("AGCAGCAUUGUACAGGGCUAUCA", {"GCAGCAU"}) is True
    assert kc.verify_arm_seed("AGCUUCUUUACAGUGUUGCCUUGU", {"GCAGCAU"}) is False
    # miR-802-5p: miRBase 23.1 dropped the 5' U TargetScan had, so the
    # family seed sits at offset 0, not 1. Still the same arm.
    assert kc.verify_arm_seed("CAGUAACAAAGAUUCAUCCUUGU", {"CAGUAAC"}) is True
    # miR-144-3p: miRBase 23.1 added 2 nt at the 5' end; the seed sits
    # at offset 3. Still the same arm.
    assert kc.verify_arm_seed("ACUACAGUAUAGAUGAUGUACU", {"ACAGUAU"}) is True
    # A seed starting at offset 4 is not annotation drift any more.
    assert kc.verify_arm_seed("UUUUACAGUAUAGAUGAUGUACU", {"ACAGUAU"}) is False
    # No seed data at all: cannot verify.
    assert kc.verify_arm_seed("AGCAGCAUUGUACAGGGCUAUCA", set()) is False
    assert kc.verify_arm_seed("", {"GCAGCAU"}) is False


# ── Engine behaviour with the real data ─────────────────────────────────────

@_needs_arm_data()
def test_mir107_3p_keeps_targets_5p_is_blocked():
    ke.TARGETSCAN_DB = ke.load_local_data()
    assert len(ke.get_targets_local("hsa-miR-107-3p")) == 836
    assert ke.get_targets_local("hsa-miR-107-5p") == set()
    check = ke.check_arm_resolution("hsa-miR-107-5p")
    assert check["status"] == "unverified"
    assert check["family_key"] == "mir-107"
    check3 = ke.check_arm_resolution("hsa-miR-107-3p")
    assert check3["status"] == "verified"


@_needs_arm_data()
def test_no_hairpin_returns_same_set_for_both_arms():
    ke.TARGETSCAN_DB = ke.load_local_data()
    for hairpin in SWAPPED_HAIRPINS:
        a = ke.get_targets_local(f"hsa-miR-{hairpin}-5p")
        b = ke.get_targets_local(f"hsa-miR-{hairpin}-3p")
        # The silent swap is gone: the two arms never return the same
        # non-empty set. And the family is not lost entirely: exactly
        # one arm keeps the modelled targets.
        assert not (a and a == b), f"miR-{hairpin}: both arms return the same set"
        assert a or b, f"miR-{hairpin}: both arms empty, family lost"


@_needs_arm_data()
def test_arm_unverified_reported_in_analyze():
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-107-3p", "hsa-miR-107-5p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    body = r.json()
    assert body["arm_verification"] == "active"
    flagged = [e["mirna"] for e in body["arm_unverified"]]
    assert flagged == ["hsa-miR-107-5p"]
    assert "hsa-miR-107-3p" in body["found_mirnas"]
    # The blocked arm falls back to "no data anywhere", never to the
    # other arm's targets.
    assert "hsa-miR-107-5p" in body["missing_mirnas"]


@_needs_arm_data()
def test_thesis_panel_unchanged():
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": LIPID_PANEL, "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    body = r.json()
    assert sorted(body["common_genes"]) == sorted(LIPID_GENES)
    assert body["arm_unverified"] == []


@_needs_arm_data()
def test_letter_drop_collapse_is_blocked():
    """Second bug class the seed check catches: the m_no_letter
    fallback drops the a/b/c letter and collapses distinct paralogs
    into the wrong family (hsa-miR-519e-3p into the miR-17 family,
    1,385 genes in v1.56). The seed does not match, so it is blocked."""
    ke.TARGETSCAN_DB = ke.load_local_data()
    assert ke.get_targets_local("hsa-miR-519e-3p") == set()
    assert ke.check_arm_resolution("hsa-miR-519e-3p")["status"] == "unverified"
    # And the letter-drop cases from the audit (122b, 151b, 219b, 302f,
    # 320a, 323b, 371b, 374c, 450b, 451b, 892a, 892b) are all blocked too.
    for m in ["hsa-miR-122b-5p", "hsa-miR-151b-3p", "hsa-miR-219b-5p",
              "hsa-miR-302f-3p", "hsa-miR-320a-5p", "hsa-miR-323b-3p",
              "hsa-miR-371b-5p", "hsa-miR-374c-5p", "hsa-miR-450b-5p",
              "hsa-miR-451b-3p", "hsa-miR-451b-5p", "hsa-miR-892a-3p",
              "hsa-miR-892b-3p"]:
        assert ke.get_targets_local(m) == set(), m


@_needs_arm_data()
def test_validated_fallback_other_arm_is_dropped():
    """The same swap lived in the miRTarBase lookup: hsa-miR-107-5p
    has no pairs of its own, so the bare-name fallback handed it the
    -3p arm's 298 pairs. The guard drops them; the -3p arm keeps them
    (its seed matches the bare key's)."""
    ke.load_mirtarbase_local()
    assert ke.get_validated_local_checked("hsa-miR-107-5p", ke.MIRTARBASE_DB) == set()
    assert len(ke.get_validated_local_checked("hsa-miR-107-3p", ke.MIRTARBASE_DB)) == 298


@_needs_arm_data()
def test_audit_zero_swaps_with_seed_criterion():
    """CI regression for the audit: over every armed human miRNA with
    local targets (332 names), a non-empty result must carry the
    family's seed in the query sequence (offset 0-3). The name-based
    audit in scripts/audit_arm_resolution.py reports the same count."""
    ke.TARGETSCAN_DB = ke.load_local_data()
    seqs, fam_seeds = ke.arm_seed_maps()
    assert seqs and fam_seeds
    names = []
    for line in open("local_db/mirbase_mature.fa"):
        if line.startswith(">hsa-"):
            names.append(line[1:].split()[0])
    checked = 0
    for n in sorted(set(names)):
        m = re.match(r"^hsa-(?:miR|let)-(\d+[a-z]?(?:-\d+)?)-(5p|3p)$", n, re.I)
        if not m:
            continue
        genes = ke.get_targets_local(n)
        if not genes:
            continue
        checked += 1
        check = ke.check_arm_resolution(n)
        assert check["status"] == "verified", (
            f"{n}: non-empty result with status {check['status']}: {check}")
    # 332 armed names had local targets in v1.56; the fix blocks 29
    # (docs/validation/arm_blocked.json): 15 unverified arms plus 14
    # letter-drop collapses (e.g. miR-519e-3p into the miR-17 family),
    # so 303 keep targets. Every kept result carries the family's seed.
    assert checked == 303
