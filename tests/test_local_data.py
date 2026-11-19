"""Tests for the local-data core: miRTarBase GMT, Enrichr GMTs,
Fisher/BH enrichment, structured export, heatmap, reverse lookup,
isomiR normalization. Pure functions, no network (the real downloaded
files in local_db/ are the fixtures).
"""

import base64
from pathlib import Path

import kenryu_core as kc

REPO = Path(__file__).resolve().parent.parent
LOCAL_DB = REPO / "local_db"
MIRTARBASE_GMT = LOCAL_DB / "mirtarbase_gene_set_library.gmt.gz"
KEGG_GMT = LOCAL_DB / "gmt_KEGG_2021_Human.gmt"


def test_load_mirtarbase_gmt_real_file():
    db = kc.load_mirtarbase_gmt(MIRTARBASE_GMT)
    assert len([k for k in db if k.startswith("hsa-")]) >= 500
    mir21 = db["hsa-miR-21-5p"]
    for expected in ("PTEN", "PDCD4", "AGO4", "MEIS1", "TGFBR2"):
        assert expected in mir21
    assert db["hsa-let-7a-5p"]  # non-empty
    assert "hsa-miR-33a-5p" in db


def test_get_validated_local():
    db = {"hsa-miR-21-5p": {"PTEN", "PDCD4"}, "hsa-let-7a": {"LIN28A"}}
    assert kc.get_validated_local("hsa-miR-21-5p", db) == {"PTEN", "PDCD4"}
    # Suffix fallback: hsa-let-7a-5p resolves to the hsa-let-7a entry.
    assert kc.get_validated_local("hsa-let-7a-5p", db) == {"LIN28A"}
    assert kc.get_validated_local("hsa-miR-999-5p", db) == set()
    assert kc.get_validated_local("hsa-miR-21-5p", {}) == set()


def test_load_gmt_library_real_file():
    lib = kc.load_gmt_library(KEGG_GMT)
    assert len(lib) > 100
    assert any("TGF-beta" in t for t in lib)
    assert all(isinstance(v, set) and v for v in lib.values())


def test_benjamini_hochberg():
    adj = kc.benjamini_hochberg([0.01, 0.04, 0.03])
    assert adj == [0.03, 0.04, 0.04]
    assert kc.benjamini_hochberg([]) == []
    assert kc.benjamini_hochberg([0.5]) == [0.5]


def test_fisher_enrichment_synthetic():
    library = {
        "Path1": {"A", "B", "C", "D", "E"},
        "Path2": {"X", "Y", "Z"},
        "Path3": {"A", "Q"},
    }
    rows = kc.fisher_enrichment({"A", "B", "C"}, library,
                                background={f"G{i}" for i in range(100)} | {"A", "B", "C", "D", "E", "X", "Y", "Z", "Q"})
    terms = [r["Term"] for r in rows]
    assert "Path1" in terms
    assert "Path2" not in terms  # zero overlap never enriches
    assert rows[0]["Term"] == "Path1"
    assert rows[0]["Overlap"] == 3
    assert rows[0]["Pval"] < 0.05


def test_local_enrichment_real_gmt():
    lib = kc.load_gmt_library(KEGG_GMT)
    rows = kc.fisher_enrichment(
        {"PTEN", "PDCD4", "TGFBR2", "TGFBR1", "SMAD4", "SMAD7", "AKT2", "PIK3R1"},
        lib, top_n=50)
    assert rows
    assert any("TGF-beta" in r["Term"] and r["Pval"] < 0.05 for r in rows)
    # Sorted by adjusted p, ascending.
    pvals = [r["Pval"] for r in rows]
    assert pvals == sorted(pvals)
    combined = kc.local_enrichment(
        {"PTEN", "PDCD4", "TGFBR2", "TGFBR1", "SMAD4", "SMAD7", "AKT2", "PIK3R1"},
        {"KEGG": lib})
    assert combined and all(r["Source"] == "KEGG" for r in combined)


def test_build_export():
    export = kc.build_export(
        mirnas=["hsa-miR-21-5p", "hsa-miR-155-5p"],
        common_genes=["PTEN", "AGO4"],
        gene_details={"PTEN": {"full_name": "phosphatase and tensin homolog"}},
        enrichment=[{"Term": "TGF-beta signaling pathway", "Source": "KEGG",
                     "Pval": 0.028}],
        validation={"genes": {"PTEN": {"score": 80, "level": "alta",
                                       "layers": ["TargetScan (predicho)"]},
                             "AGO4": {"score": 55, "level": "media",
                                      "layers": []}},
                    "data_sources": {"mirbase": "23.1"},
                    "summary": {"n_genes": 2}},
        gene_sets=[{"PTEN", "AGO4"}, {"PTEN"}],
    )
    assert export["schema_version"] == "1.0"
    assert export["panel"] == ["hsa-miR-21-5p", "hsa-miR-155-5p"]
    assert export["genes"][0]["gene"] == "PTEN"
    assert export["genes"][0]["validation_score"] == 80
    assert {"mirna": "hsa-miR-21-5p", "gene": "PTEN",
            "evidence": "TargetScan"} in export["edges"]
    assert {"mirna": "hsa-miR-155-5p", "gene": "AGO4",
            "evidence": "TargetScan"} not in export["edges"]
    assert export["enrichment"][0]["p_adj"] == 0.028
    assert export["validation_summary"]["n_genes"] == 2


def test_heatmap_matrix_and_png():
    genes, mirnas, matrix = kc.heatmap_matrix(
        ["m1", "m2"], ["G1", "G2"], [{"G1"}, {"G1", "G2"}])
    assert genes == ["G1", "G2"]
    assert matrix == [[1, 1], [0, 1]]
    png = kc.create_heatmap(["m1", "m2"], ["G1", "G2"], [{"G1"}, {"G1", "G2"}])
    assert base64.b64decode(png)[:8] == b"\x89PNG\r\n\x1a\n"
    assert kc.create_heatmap([], [], []) is None


def test_reverse_lookup():
    db = {"hsa-miR-21-5p": {"PTEN", "PDCD4"}, "hsa-miR-155-5p": {"PTEN"}}
    assert kc.reverse_lookup("PTEN", db) == ["hsa-miR-155-5p", "hsa-miR-21-5p"]
    assert kc.reverse_lookup("PDCD4", db) == ["hsa-miR-21-5p"]
    assert kc.reverse_lookup("NADIE", db) == []
    assert kc.reverse_lookup("", db) == []


def test_normalize_isomir():
    assert kc.normalize_isomir("hsa-miR-21-5p+2") == "hsa-miR-21-5p"
    assert kc.normalize_isomir("hsa-miR-21-5p-1") == "hsa-miR-21-5p"
    assert kc.normalize_isomir("hsa-miR-21-5p|5p-shift") == "hsa-miR-21-5p"
    assert kc.normalize_isomir("hsa-let-7a-5p+1") == "hsa-let-7a-5p"
    assert kc.normalize_isomir("hsa-miR-21-5p") == "hsa-miR-21-5p"
    # Regression: bare numbers are part of the name, not a trimming mark.
    assert kc.normalize_isomir("hsa-miR-375") == "hsa-miR-375"
    assert kc.normalize_isomir("hsa-miR-33a-5p") == "hsa-miR-33a-5p"
    assert kc.normalize_isomir("hsa-miR-124-3p") == "hsa-miR-124-3p"


def test_deep_clean_species_prefix():
    # One species prefix is stripped; "mir-" itself is never stripped.
    assert kc.deep_clean_mirna("hsa-miR-33a-5p") == "mir-33a-5p"
    assert kc.deep_clean_mirna("mmu-miR-21a-5p") == "mir-21a-5p"
    assert kc.deep_clean_mirna("mir-33a-5p") == "mir-33a-5p"
    assert kc.deep_clean_mirna("hsa-let-7a-5p") == "let-7a-5p"
    assert kc.clean_mirna_name_cloud("mmu-miR-21a-5p", "mmu") == "mmu-miR-21a-5p"
    assert kc.normalize_isomir("mmu-miR-21a-5p", "mmu") == "mmu-miR-21a-5p"


def test_heatmap_union_varies():
    """The union heatmap must VARY: the old core-only matrix was all 1s by
    construction (a solid teal block that explained nothing)."""
    sets = [{"A", "B", "C"}, {"A", "B", "D"}, {"A", "E"}]
    genes, mirnas, matrix, union_size = kc.heatmap_matrix_union(
        ["m1", "m2", "m3"], sets)
    assert union_size == 5
    assert genes[0] == "A"  # targeted by all 3 -> on top
    assert any(0 in row for row in matrix)  # real variation
    assert matrix[genes.index("A")] == [1, 1, 1]
    assert matrix[genes.index("E")] == [0, 0, 1]
    png = kc.create_heatmap(["m1", "m2", "m3"], ["A"], sets)
    assert base64.b64decode(png)[:8] == b"\x89PNG\r\n\x1a\n"
