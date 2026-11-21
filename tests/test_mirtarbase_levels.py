"""miRTarBase level policy for v2.0 (miRTarBase 10.0 adoption, step D).

v10.0 adds 1,411,229 weak pairs with NO PMID — all of them CLIP-only
(measured 2026-11-21: weak_no_pmid = 1,411,229, non-CLIP = 0, and
strong_no_pmid = 0 in both versions). They get their own level:

- "strong" (Functional MTI, with PMID): +30, unlocks "alta".
- "weak" (Functional MTI (Weak), with PMID): +15, never "alta".
- "clip_only" (Functional MTI (Weak), NO PMID, all experiments CLIP):
  visible in the dossier, +0 points, never unlocks "alta".
- Bare miRNA names (no -5p/-3p arm) are marked "brazo no
  especificado" wherever their evidence is shown.

These tests are written fail-first: they fail on v1.60 (there is no
"clip_only" level, the constructor maps every weak row to "weak"
regardless of PMID, and nothing marks bare-name evidence).
"""
import csv
import gzip
import json

import pytest

import kenryu_core as kc

COLUMNS = ["miRTarBase ID", "miRNA", "Species (miRNA)", "Target Gene",
           "Target Gene (Entrez ID)", "Species (Target Gene)",
           "Experiments", "Support Type", "References (PMID)"]


def _write_csv(tmp_path, rows):
    path = tmp_path / "hsa_MTI.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLUMNS})
    return path


def _row(mirna, gene, support, pmid="", experiments=""):
    return {"miRNA": mirna, "Target Gene": gene, "Support Type": support,
            "References (PMID)": pmid, "Experiments": experiments}


def _build(tmp_path, rows, out_name="out.json.gz"):
    import build_mirtarbase_official as bmo
    src = _write_csv(tmp_path, rows)
    out = tmp_path / out_name
    return bmo.build(str(src), str(out), version="miRTarBase 10.0 (synthetic)")


def test_constructor_assigns_clip_only(tmp_path):
    """Weak rows without a PMID become "clip_only" (+0), not "weak"."""
    payload = _build(tmp_path, [
        _row("hsa-miR-1-5p", "GENEA", "Functional MTI", "111",
             "Luciferase reporter assay"),
        _row("hsa-miR-1-5p", "GENEB", "Functional MTI (Weak)", "222",
             "HITS-CLIP"),
        _row("hsa-miR-1-5p", "GENEC", "Functional MTI (Weak)", "",
             "HITS-CLIP"),
        _row("hsa-miR-1-5p", "GENED", "Functional MTI (Weak)", "",
             "HITS-CLIP;PAR-CLIP"),
    ])
    mirnas = payload["mirnas"]["hsa-miR-1-5p"]
    assert mirnas["GENEA"]["level"] == "strong"
    assert mirnas["GENEA"]["pmids"] == [111]
    assert mirnas["GENEB"]["level"] == "weak"
    assert mirnas["GENEB"]["pmids"] == [222]
    assert mirnas["GENEC"]["level"] == "clip_only"
    assert mirnas["GENEC"]["pmids"] == []
    assert mirnas["GENEC"]["experiments"] == ["HITS-CLIP"]
    assert mirnas["GENED"]["level"] == "clip_only"
    assert payload["levels"] == {"strong": 1, "weak": 1, "clip_only": 2}
    assert payload["version"] == "miRTarBase 10.0 (synthetic)"


def test_constructor_refuses_weak_without_pmid_and_not_clip(tmp_path):
    """A weak pair without PMID whose experiments are not all CLIP is
    an unknown category: the build must fail, not guess a level."""
    with pytest.raises(ValueError, match="sin PMID"):
        _build(tmp_path, [
            _row("hsa-miR-2-5p", "GENEE", "Functional MTI (Weak)", "",
                 "Microarray"),
        ])


def test_constructor_refuses_strong_without_pmid(tmp_path):
    """A Functional MTI row without PMID cannot be called "strong"."""
    with pytest.raises(ValueError, match="sin PMID"):
        _build(tmp_path, [
            _row("hsa-miR-2-5p", "GENEF", "Functional MTI", "",
                 "Luciferase reporter assay"),
        ])


def test_constructor_refuses_unknown_support_type(tmp_path):
    """A Support Type the builder does not know stops the build."""
    with pytest.raises(ValueError, match="Support Type desconocido"):
        _build(tmp_path, [
            _row("hsa-miR-2-5p", "GENEG", "Functional MTI (Super)", "333",
                 "Luciferase reporter assay"),
        ])


def test_clip_only_scores_zero_and_never_unlocks_alta():
    """clip_only adds 0 points and a gene cannot reach "alta" on it,
    even when every other layer is maxed (25 + 15 + 10 + 10 + 10 = 70)."""
    research = {"pubmed": [1, 2, 3], "omim": ["MIM:1"],
                "clinvar": {"count": 1}, "trials": ["NCT1"]}
    result = kc.validation_score(
        "GENEX", research, predicted=True, validated="clip_only")
    assert result["score"] == 70
    assert result["level"] == "media"
    assert result["evidence_level"] == "clip_only"
    assert any("CLIP sin PMID" in layer for layer in result["layers"])
    assert "solo CLIP sin PMID" in result["explanation"]


def test_evidence_tiers_clip_only_is_visible_not_validated():
    """The dossier shows clip_only evidence as present-but-unpointed:
    visible, never called "validado"."""
    tiers = kc.build_evidence_tiers(True, "clip_only", {})
    assert tiers["experimental"]["present"] is True
    summary = tiers["experimental"]["summary"]
    assert "sin PMID" in summary
    assert "sin puntos" in summary
    assert "Validado" not in summary


def test_best_evidence_rank_clip_only_between_none_and_weak():
    """clip_only beats "none" (it is evidence) and loses to "weak"."""
    by_mirna = {
        "hsa-miR-1-5p": {"GENEA": {"level": "clip_only", "pmids": [],
                                   "experiments": ["HITS-CLIP"]}},
        "hsa-miR-2-5p": {"GENEA": {"level": "none", "pmids": [999],
                                   "experiments": []}},
    }
    assert kc._best_evidence("GENEA", by_mirna)["level"] == "clip_only"
    by_mirna["hsa-miR-2-5p"]["GENEA"] = {"level": "weak", "pmids": [999],
                                         "experiments": ["HITS-CLIP"]}
    assert kc._best_evidence("GENEA", by_mirna)["level"] == "weak"


def test_validated_local_excludes_clip_only():
    """clip_only pairs do not put a gene in the validated set."""
    db = {"hsa-miR-1-5p": {
        "GENEA": {"level": "clip_only", "pmids": [],
                  "experiments": ["HITS-CLIP"]},
        "GENEB": {"level": "weak", "pmids": [222],
                  "experiments": ["HITS-CLIP"]},
    }}
    assert kc.get_validated_local("hsa-miR-1-5p", db) == {"GENEB"}


def test_validated_key_source_marks_bare_names():
    """Evidence from a bare-name index key is marked; an armed key is
    not; a missing miRNA has no source."""
    db = {"hsa-miR-758": {"ABCA1": {"level": "weak", "pmids": [],
                                    "experiments": ["CLIP-seq"]}},
          "hsa-miR-144-3p": {"ABCA1": {"level": "weak", "pmids": [1],
                                       "experiments": ["HITS-CLIP"]}}}
    assert kc.validated_key_source("hsa-miR-758-3p", db) == "bare"
    assert kc.validated_key_source("hsa-miR-144-3p", db) == "armed"
    assert kc.validated_key_source("hsa-miR-999-5p", db) is None
    # A bare query against a bare key is bare evidence too: the entry
    # does not say which arm (v2.0, paso D).
    assert kc.validated_key_source("hsa-miR-375", db) is None
    db_bare = {"hsa-miR-375": {"NFIA": {"level": "weak", "pmids": [],
                                        "experiments": ["CLIP-seq"]}}}
    assert kc.validated_key_source("hsa-miR-375", db_bare) == "bare"


def test_engine_marks_bare_name_evidence():
    """The engine flags evidence that came from a bare-name key so the
    dossier can label it "brazo no especificado" — both when an armed
    query falls back to the bare key and when the query itself is bare."""
    import kenryu_engine as ke
    db = {"hsa-miR-758": {"ABCA1": {"level": "clip_only", "pmids": [],
                                    "experiments": ["CLIP-seq"]}}}
    ev = ke.get_validated_evidence_checked("hsa-miR-758-3p", db)
    assert ev["ABCA1"]["level"] == "clip_only"
    assert ev["ABCA1"]["bare_name"] is True
    db2 = {"hsa-miR-758-3p": {"ABCA1": {"level": "weak", "pmids": [1],
                                        "experiments": ["HITS-CLIP"]}}}
    ev2 = ke.get_validated_evidence_checked("hsa-miR-758-3p", db2)
    assert ev2["ABCA1"].get("bare_name") is not True
    db3 = {"hsa-miR-375": {"NFIA": {"level": "clip_only", "pmids": [],
                                    "experiments": ["CLIP-seq"]}}}
    ev3 = ke.get_validated_evidence_checked("hsa-miR-375", db3)
    assert ev3["NFIA"]["bare_name"] is True


def test_evidence_table_carries_bare_name_flag():
    """The per-miRNA evidence table in the dossier carries the flag."""
    by_mirna = {"hsa-miR-758-3p": {
        "ABCA1": {"level": "clip_only", "pmids": [],
                  "experiments": ["CLIP-seq"], "bare_name": True}}}
    table = kc._evidence_by_mirna_table(
        "ABCA1", by_mirna, ["hsa-miR-758-3p"])
    assert table["hsa-miR-758-3p"]["bare_name"] is True
    assert table["hsa-miR-758-3p"]["level"] == "clip_only"
