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
REACTOME_GMT = LOCAL_DB / "gmt_Reactome_Pathways_2024.gmt"


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
    lib = kc.load_gmt_library(REACTOME_GMT)
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
    # v1.45: la GMT de KEGG salió del repo por licencia (ver
    # DATA_LICENSES.md); la prueba usa Reactome, que sigue en local_db/.
    lib = kc.load_gmt_library(REACTOME_GMT)
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
        {"Reactome": lib})
    assert combined and all(r["Source"] == "Reactome" for r in combined)


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


def test_fisher_enrichment_bh_counts_zero_overlap_terms():
    """v1.45 (auditoría externa, H-05): BH se calculaba solo sobre los
    términos con solape, así que m era ~28× menor que el total en GO y los
    p ajustados salían demasiado optimistas. Añadir términos sin solape
    debe hacer el ajuste MÁS conservador (m mayor), no dejarlo igual."""
    universe = {f"G{i}" for i in range(1000)}
    query = {f"Q{i}" for i in range(10)}
    enriched = {"Enriched": set(query)}
    others = {f"Other{i}": {f"Q{i}", f"G{i}"} | {f"G{100 + i * 3 + j}" for j in range(3)}
              for i in range(9)}
    lib_a = {**enriched, **others}
    lib_b = dict(lib_a)
    for i in range(100):
        lib_b[f"Null{i}"] = {f"G{500 + i * 3 + j}" for j in range(5)}
    rows_a = kc.fisher_enrichment(query, lib_a, background=universe | query)
    rows_b = kc.fisher_enrichment(query, lib_b, background=universe | query)
    pa = next(r["Pval"] for r in rows_a if r["Term"] == "Enriched")
    pb = next(r["Pval"] for r in rows_b if r["Term"] == "Enriched")
    assert pb > pa  # m mayor → ajuste más conservador
    assert pa < 0.05 and pb < 0.05  # el término real sigue siendo significativo


def test_fisher_enrichment_ignores_genes_outside_universe():
    """v1.45 (H-05): genes que no están en la biblioteca GMT inflaban b y
    cambiaban el p ajustado (la auditoría vio 7.75e-49 → 1.20e-44 al añadir
    5 símbolos inexistentes). Ahora se recorta la lista al universo y el
    resultado es idéntico."""
    library = {
        "Path1": {"A", "B", "C", "D", "E"},
        "Path2": {"X", "Y", "Z"},
        "Path3": {"A", "Q"},
    }
    bg = {f"G{i}" for i in range(100)} | {"A", "B", "C", "D", "E", "X", "Y", "Z", "Q"}
    base = kc.fisher_enrichment({"A", "B", "C"}, library, background=bg)
    with_ghosts = kc.fisher_enrichment(
        {"A", "B", "C", "ZZ1", "ZZ2", "ZZ3", "ZZ4", "ZZ5"}, library, background=bg)
    assert [(r["Term"], r["Pval"]) for r in base] == [(r["Term"], r["Pval"]) for r in with_ghosts]


def test_fisher_enrichment_null_list_has_no_significant_terms():
    """v1.45 (H-05): 40/40 listas nulas daban ≥1 término con p_adj < 0.3
    con el BH viejo. Con m = total de términos y corte 0.05, una lista
    sin enriquecimiento real no devuelve ningún término significativo."""
    import random
    rng = random.Random(123)
    universe = {f"G{i}" for i in range(1000)}
    library = {}
    for i in range(200):
        library[f"Term{i}"] = {f"G{rng.randrange(1000)}" for _ in range(20)}
    query = {f"G{i}" for i in range(0, 1000, 25)}  # 40 genes repartidos
    rows = kc.fisher_enrichment(query, library, background=universe, top_n=50)
    assert rows == []


def test_prioritize_genes_orders_by_regulation_then_validation():
    """v1.46 (H-06): el orden no es alfabético: primero los genes que
    regulan más miRNAs del panel y, en empate, los validados en miRTarBase
    por más miRNAs. El nombre solo desempata (determinismo)."""
    gene_sets = [{"AGO4", "FOXP1", "BNC2"}, {"AGO4", "FOXP1"}, {"AGO4"}]
    remote = {"m1": {"AGO4"}, "m2": {"AGO4", "BNC2"}, "m3": set()}
    ordered = kc.prioritize_genes(["FOXP1", "BNC2", "AGO4"], gene_sets, remote)
    assert ordered == ["AGO4", "FOXP1", "BNC2"]


def test_build_validation_report_marks_unresearched_as_not_evaluated():
    """v1.46 (H-06): un gen fuera del tope de investigación no es 'baja'
    por orden alfabético: es 'no evaluado', sin puntaje, y no entra en la
    media ni en los conteos."""
    local = {"hsa-miR-21-5p": {"AGO4", "FOXP1"}, "hsa-miR-155-5p": {"AGO4", "FOXP1"}}
    remote = {"hsa-miR-21-5p": {"AGO4"}, "hsa-miR-155-5p": {"AGO4"}}
    research = {"AGO4": {"pubmed": [{"pmid": "12345678", "term": "x", "source": "KEGG",
                                     "editorial_note": "", "year_window": "10y"}],
                         "omim": [], "clinvar": {"count": 0, "url": ""},
                         "trials": [], "diseases": []}}
    # FOXP1 no tiene entrada: no investigado.
    rep = kc.build_validation_report(
        ["hsa-miR-21-5p", "hsa-miR-155-5p"], ["AGO4", "FOXP1"], research, local, remote, {})
    assert rep["genes"]["FOXP1"]["level"] == "not_evaluated"
    assert rep["genes"]["FOXP1"]["score"] is None
    assert rep["genes"]["AGO4"]["level"] == "media"
    assert rep["genes"]["AGO4"]["score"] == 60
    assert rep["summary"]["n_not_evaluated"] == 1
    assert rep["summary"]["not_evaluated"] == ["FOXP1"]
    assert "FOXP1" not in rep["summary"]["not_validated"]
    assert rep["summary"]["mean_score"] == 60.0


def test_build_validation_report_labels_mean_as_evaluated_only():
    """v1.46 (ajuste pre-publicación): la media es solo de los genes
    evaluados y la etiqueta lo dice con el n, para que nadie lea
    'media 80' como 'este panel tiene evidencia fuerte'."""
    local = {"hsa-miR-21-5p": {"AGO4", "FOXP1"}, "hsa-miR-155-5p": {"AGO4", "FOXP1"}}
    remote = {"hsa-miR-21-5p": {"AGO4"}, "hsa-miR-155-5p": {"AGO4"}}
    research = {"AGO4": {"pubmed": [{"pmid": "12345678", "term": "x", "source": "KEGG",
                                     "editorial_note": "", "year_window": "10y"}],
                         "omim": [], "clinvar": {"count": 0, "url": ""},
                         "trials": [], "diseases": []}}
    rep = kc.build_validation_report(
        ["hsa-miR-21-5p", "hsa-miR-155-5p"], ["AGO4", "FOXP1"], research, local, remote, {})
    assert rep["summary"]["mean_score"] == 60.0
    assert rep["summary"]["mean_score_label"] == "Media de los genes evaluados (1 de 2)"
    assert "priorizados" in rep["method"] or "prioridad" in rep["method"]
    assert "miRTarBase" in rep["method"]
