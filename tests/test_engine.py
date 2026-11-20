"""Offline test suite for the KENRYU engine (kenryu_engine.py).

Every test runs without network access, using only the data files
shipped in the repository (data/hsa-miR-*.txt). Nothing here modifies
the engine: tests assert the behaviour the thesis documents.

Run:  pytest -v
"""
import base64
import json

import pytest
from fastapi.testclient import TestClient

import kenryu_engine as ke


# ── Name cleaning ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("raw,expected", [
    ("hsa-miR-33a-5p", "mir-33a-5p"),
    ("HSA-MIR-21-5P", "mir-21-5p"),
    ("miR-144-3p", "mir-144-3p"),
    ("hsa-hsa-miR-375", "mir-375"),
    ("miR\u201133a\u20115p", "mir-33a-5p"),  # unicode dashes
])
def test_deep_clean_mirna(raw, expected):
    assert ke.deep_clean_mirna(raw) == expected


@pytest.mark.parametrize("raw,expected", [
    ("mir-33a-5p", "hsa-miR-33a-5p"),
    ("hsa-miR-21-5p", "hsa-miR-21-5p"),
    ("miR-144-3p", "hsa-miR-144-3p"),
    ("mir-106b-5p", "hsa-miR-106b-5p"),
    # let-7 family keeps its name: hsa-let-7a-5p, never hsa-miR-7a-5p.
    ("let-7a-5p", "hsa-let-7a-5p"),
    ("hsa-let-7a-5p", "hsa-let-7a-5p"),
    ("mir-let-7a-5p", "hsa-let-7a-5p"),
    ("let-7a", "hsa-let-7a"),
])
def test_clean_mirna_name_cloud(raw, expected):
    assert ke.clean_mirna_name_cloud(raw) == expected


# ── Local TargetScan lookup (repo data files, no network) ────────────────────

PANEL = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
         "hsa-miR-106b-5p", "hsa-miR-33b-5p"]

# Gene counts per file, as shipped in data/ (thesis §4.1: miR-33-5p -> 545).
EXPECTED_COUNTS = {
    "hsa-miR-33a-5p": 545,
    "hsa-miR-33b-5p": 545,
    "hsa-miR-758-3p": 271,
    "hsa-miR-144-3p": 1048,
    "hsa-miR-106b-5p": 1385,
}


@pytest.mark.parametrize("mirna,count", list(EXPECTED_COUNTS.items()))
def test_get_targets_local_counts(mirna, count):
    targets = ke.get_targets_local(mirna)
    assert len(targets) == count


def test_get_targets_local_contains_abca1():
    # ABCA1 is a validated target of every miRNA in the lipid panel.
    for mirna in PANEL:
        assert "ABCA1" in ke.get_targets_local(mirna)


def test_strict_intersection_reproduces_thesis():
    """The Strict-mode core of the lipid panel must be exactly the five
    convergent genes of the UNIR thesis (Table 4.1 / §4.2)."""
    sets = [ke.get_targets_local(m) for m in PANEL]
    assert set.intersection(*sets) == {"ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"}


def test_get_targets_local_unknown_returns_empty():
    assert ke.get_targets_local("hsa-miR-9999-5p") == set()


# ── Request model ────────────────────────────────────────────────────────────

def test_analysis_request_defaults():
    req = ke.AnalysisRequest(mirnas=PANEL)
    assert req.years == 10
    assert req.mode == "strict"
    assert req.month is None


# ── Academic synthesis (pure function) ───────────────────────────────────────

def test_build_synthesis_structure():
    common = ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]
    gene_details = {
        "ABCA1": {
            "full_name": "ATP Binding Cassette Subfamily A Member 1",
            "pathology": "Regulador maestro del eflujo de colesterol.",
            "clinical_evidence": [
                {"source": "KEGG", "term": "Metabolismo del colesterol",
                 "desc": "miR-33 suprime ABCA1.", "pmid": None},
            ],
        }
    }
    enrichment = [{"Term": "Metabolismo del colesterol", "Source": "KEGG",
                   "Pval": 0.016, "ScientificDesc": "Ruta identificada con p-adj=0.016"}]
    references = [{"id": 1, "authors": "Ramirez et al.", "title": "miR-33 and ABCA1",
                   "journal": "Nature", "year": "2020", "url": "https://pubmed.ncbi.nlm.nih.gov/30212345/",
                   "pmid": "30212345", "type": "pubmed", "source": "KEGG"}]
    out = ke.build_synthesis(PANEL, common, gene_details, enrichment, references)
    assert set(out) == {"academic", "functional", "references"}
    assert "5 genes" in out["academic"]
    assert "ABCA1" in out["academic"]
    assert "30212345" in out["references"]


# ── Visualizations (pure, no network) ────────────────────────────────────────

def test_create_visuals_real_data_only():
    """Venn is always real (sets). Volcano needs real PubMed counts and PPI
    needs real STRING edges; without them the plot is None, never simulated."""
    sets = [{"ABCA1", "KPNA3"}, {"ABCA1", "SCN1A"}, {"ABCA1", "SNTB2"}]
    counts = {"ABCA1": 1200, "KPNA3": 45, "SCN1A": 900, "SNTB2": 30}
    edges = [("ABCA1", "KPNA3", 0.7), ("SCN1A", "SNTB2", 0.9)]
    venn, volcano, ppi = ke.create_visuals(
        sets, PANEL[:3], ["ABCA1"], pubmed_counts=counts, ppi_edges=edges)
    for img in (venn, volcano, ppi):
        raw = base64.b64decode(img)
        assert raw[:8] == b"\x89PNG\r\n\x1a\n"  # PNG magic bytes
    # Without real data: volcano and PPI are None (honest, not simulated).
    venn2, volcano2, ppi2 = ke.create_visuals(sets, PANEL[:3], ["ABCA1"])
    assert base64.b64decode(venn2)[:8] == b"\x89PNG\r\n\x1a\n"
    assert volcano2 is None
    assert ppi2 is None
    # The volcano depends on the data: different counts -> different image.
    counts2 = dict(counts, ABCA1=5)
    _, volcano_b, _ = ke.create_visuals(
        sets, PANEL[:3], ["ABCA1"], pubmed_counts=counts2, ppi_edges=edges)
    assert volcano_b != volcano


# ── Consensus modes (pure function, extracted from analyze()) ───────────────

def test_consensus_genes_strict():
    sets = [{"A", "B"}, {"A", "C"}, {"A", "D"}]
    assert ke.consensus_genes(sets, "strict") == {"A"}


def test_consensus_genes_n1():
    sets = [{"A", "B"}, {"A", "C"}, {"B", "C"}]
    assert ke.consensus_genes(sets, "n-1") == {"A", "B", "C"}


def test_consensus_genes_n2():
    sets = [{"A", "B"}, {"A", "C"}, {"A", "D"}]
    assert ke.consensus_genes(sets, "n-2") == {"A", "B", "C", "D"}


def test_consensus_genes_single_panel_behaves_as_strict():
    assert ke.consensus_genes([{"A", "B"}], "n-1") == {"A", "B"}


# ── API validation ───────────────────────────────────────────────────────────

def test_analyze_rejects_empty_panel():
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={"mirnas": []})
    assert r.status_code == 422


# ── HTTP health + HEAD endpoints ───────────────────────────────────────────────

def test_head_root_returns_200():
    # Backported from the Space: the root serves HEAD so uptime checks and
    # the Hugging Face platform can probe / without a body.
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.head("/")
    assert r.status_code == 200


# ── HTTP health endpoint ─────────────────────────────────────────────────────

def test_health_endpoint():
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.get("/health")
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "OK"
    assert body["version"] == "1.55"
    assert body["db_families"] >= 5  # fallback: the 5 panel .txt files
    assert body["raw_files"] == 5


# ── Protein relevance (UniProt) + STRING neighbours ──────────────────────────

def test_parse_uniprot_entry():
    entry = {
        "primaryAccession": "O95477",
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Phospholipid-transporting ATPase ABCA1"}}},
        "comments": [
            {"commentType": "FUNCTION", "texts": [{"value": "Catalyzes the translocation of phospholipids."}]},
            {"commentType": "SUBCELLULAR LOCATION", "subcellularLocations": [{"location": {"value": "Cell membrane"}}]},
            {"commentType": "DISEASE", "disease": {"diseaseId": "Tangier disease", "description": "Deficiency causes Tangier."}},
        ],
    }
    p = ke.parse_uniprot_entry(entry, "ABCA1")
    assert p["uniprot_id"] == "O95477"
    assert p["protein_name"] == "Phospholipid-transporting ATPase ABCA1"
    assert "translocation" in p["function"]
    assert p["locations"] == ["Cell membrane"]
    assert p["diseases"][0].startswith("Tangier disease")
    assert p["url"] == "https://rest.uniprot.org/uniprotkb/O95477"
    assert ke.parse_uniprot_entry({}, "X") is None


def test_fetch_string_network_add_nodes_keeps_neighbours():
    import asyncio
    import httpx

    def _mock(request):
        assert "add_nodes=5" in str(request.url)
        return httpx.Response(200, json=[
            {"preferredName_A": "ABCA1", "preferredName_B": "SNTB2", "score": 0.836},
            {"preferredName_A": "SNTA1", "preferredName_B": "ABCA1", "score": 0.568},
            {"preferredName_A": "ZZZ", "preferredName_B": "YYY", "score": 0.9},
        ])

    async def _run():
        transport = httpx.MockTransport(_mock)
        async with httpx.AsyncClient(transport=transport) as client:
            return await ke.fetch_string_network(["ABCA1", "SNTB2"], client, add_nodes=5)

    edges = asyncio.run(_run())
    assert ("ABCA1", "SNTB2", 0.836) in edges
    assert ("ABCA1", "SNTA1", 0.568) in edges  # neighbour kept
    assert all("ZZZ" not in e and "YYY" not in e for e in edges)  # unrelated dropped


def test_fetch_uniprot_for_gene_mocked():
    import asyncio
    import httpx

    def _mock(request):
        assert "ABCA1" in str(request.url)
        return httpx.Response(200, json={"results": [{
            "primaryAccession": "O95477",
            "proteinDescription": {"recommendedName": {"fullName": {"value": "Phospholipid-transporting ATPase ABCA1"}}},
            "comments": [{"commentType": "FUNCTION", "texts": [{"value": "Moves phospholipids."}]}],
        }]})

    async def _run():
        transport = httpx.MockTransport(_mock)
        async with httpx.AsyncClient(transport=transport) as client:
            return await ke.fetch_uniprot_for_gene("ABCA1", client)

    p = asyncio.run(_run())
    assert p["uniprot_id"] == "O95477"
    assert "Moves phospholipids" in p["function"]


def test_ppi_isolated_core_gene_and_neighbour():
    """Regression: a core gene with NO STRING edges (KPNA3 in the thesis
    panel) must not crash the PPI — it is drawn as an isolated node, and
    neighbour proteins (first shell) are drawn grey."""
    import base64
    sets = [{"ABCA1", "KPNA3"}, {"ABCA1", "SCN1A"}, {"ABCA1", "SNTB2"}]
    counts = {"ABCA1": 1200, "KPNA3": 45, "SCN1A": 900, "SNTB2": 30}
    edges = [("ABCA1", "SNTB2", 0.836), ("ABCA1", "SNTA1", 0.568)]
    venn, volcano, ppi = ke.create_visuals(
        sets, PANEL[:3], ["ABCA1", "KPNA3", "SCN1A", "SNTB2"],
        pubmed_counts=counts, ppi_edges=edges)
    for img in (venn, volcano, ppi):
        assert base64.b64decode(img)[:8] == b"\x89PNG\r\n\x1a\n"


def test_dedupe_pubmed_refs():
    """Regression (auditoría externa, 2026-11-18): dos términos de evidencia
    pueden devolver el mismo PMID y cada repetición sumaba 5 puntos.
    KPNA3, SCN1A, SNTB2 y TSC22D2 del panel de la tesis tenían el mismo
    PMID dos veces en vivo."""
    refs = [
        {"pmid": "32492291", "term": "term A"},
        {"pmid": "32492291", "term": "term B"},  # mismo artículo, otro término
        {"pmid": "36153432", "term": "term C"},
        {"pmid": None, "term": "sin pmid"},      # sin PMID: se conserva
    ]
    out = ke.dedupe_pubmed_refs(refs)
    assert [r["pmid"] for r in out] == ["32492291", "36153432", None]
    assert out[0]["term"] == "term A"  # se conserva la primera ocurrencia
    assert ke.dedupe_pubmed_refs([]) == []
    assert ke.dedupe_pubmed_refs(None) == []


def test_preset_gene_details_have_no_handwritten_pmids():
    """Regression (auditoría externa, 2026-11-19): the 7 preset genes carried
    11 hand-written PMIDs and EVERY one resolved in PubMed to an unrelated
    paper (e.g. PMID 25648732 is pregnancy counselling, cited for ABCA1 HDL
    assembly and for GMFB MAPK/ERK). A preset PMID must come from a live
    query or be None — never written by hand. The gene panel only shows the
    'verify evidence' link when a real PMID exists."""
    import asyncio
    import httpx

    def _mock(request):
        # UniProt returns nothing; get_gene_details must still work offline.
        return httpx.Response(200, json={"results": []})

    async def _run(gene):
        transport = httpx.MockTransport(_mock)
        async with httpx.AsyncClient(transport=transport) as client:
            return await ke.get_gene_details(gene, client)

    for gene in ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2", "GMFB", "KLF12"]:
        info = asyncio.run(_run(gene))
        assert info["pmid"] is None, gene
        assert all(ev.get("pmid") is None for ev in info["clinical_evidence"]), gene


def test_build_synthesis_separates_editorial_note_from_live_pmid():
    """v1.44 (auditoría externa): la síntesis pegaba la descripción escrita
    a mano junto a un PMID — primero uno inventado, y tras v1.43 el literal
    'PMID: None'. La nota editorial y el artículo en vivo van separados y
    se dice así: la nota sin cita, el PMID como primer resultado de la
    búsqueda por gen + término."""
    common = ["ABCA1"]
    gene_details = {
        "ABCA1": {
            "full_name": "ATP Binding Cassette Subfamily A Member 1",
            "pathology": "Regulador maestro del eflujo de colesterol.",
            "clinical_evidence": [
                {"source": "KEGG", "term": "Metabolismo del colesterol",
                 "desc": "miR-33 suprime ABCA1.", "pmid": None},
            ],
        }
    }
    gene_research = {
        "ABCA1": {
            "pubmed": [{"pmid": "12345678", "term": "Metabolismo del colesterol",
                        "source": "KEGG", "editorial_note": "miR-33 suprime ABCA1.",
                        "year_window": "10y"}],
            "omim": [], "clinvar": {"count": 0, "url": ""}, "trials": [],
        }
    }
    out = ke.build_synthesis(PANEL, common, gene_details, [], [], gene_research)
    text = out["academic"]
    assert "PMID: None" not in text
    assert "Notas editoriales" in text
    assert "Primeros resultados en PubMed" in text
    assert "12345678" in text
    # Ninguna línea junta la frase editorial con un PMID.
    for line in text.splitlines():
        if "miR-33 suprime ABCA1" in line:
            assert "PMID" not in line


def test_research_gene_complete_labels_desc_as_editorial_note():
    """v1.44: el registro de PubMed ya no lleva la descripción del preset
    en el campo 'desc' junto al PMID (eso implicaba que el artículo la
    respalda). Ahora va en 'editorial_note' y el PMID es solo el primer
    resultado de la búsqueda por gen + término."""
    import asyncio
    import httpx

    def _mock(request):
        url = str(request.url)
        params = request.url.params
        if "clinicaltrials.gov" in url:
            return httpx.Response(200, json={"studies": []})
        if "esearch.fcgi" in url:
            db = params.get("db", "")
            if db == "pubmed":
                return httpx.Response(200, json={"esearchresult": {"idlist": ["12345678"], "count": "1"}})
            return httpx.Response(200, json={"esearchresult": {"idlist": [], "count": "0"}})
        return httpx.Response(200, json={"result": {}})

    preset = {"clinical_evidence": [
        {"term": "Metabolismo del colesterol", "source": "KEGG",
         "desc": "miR-33 suprime ABCA1."},
    ]}

    async def _run():
        transport = httpx.MockTransport(_mock)
        async with httpx.AsyncClient(transport=transport) as client:
            return await ke.research_gene_complete("TESTGENE", preset, 10, client)

    res = asyncio.run(_run())
    assert len(res["pubmed"]) == 1
    ref = res["pubmed"][0]
    assert ref["pmid"] == "12345678"
    assert ref["editorial_note"] == "miR-33 suprime ABCA1."
    assert "desc" not in ref


def test_get_targets_local_fuzzy_match_is_token_bounded():
    """v1.45 (auditoría externa, H-03): la subcadena simple asignaba dianas
    de familias ajenas — '1-5p' casa dentro de 'mir-211-5p' y devolvía 793
    genes de miR-204/211 para un hsa-miR-1-5p inexistente. El número solo
    vale en límite de token: '873-5p' sí encuentra 'mir-873-5p.2'
    (variante real sin clave exacta), pero '1-5p' ya no casa con 'mir-211-5p'."""
    db = {
        "mir-204-5p/211-5p": {"GENE1", "GENE2"},
        "mir-211-5p": {"GENE3"},
        "mir-873-5p.2": {"GENE4"},
        "mir-17-5p/20-5p/93-5p/106-5p/519-3p": {"GENE5"},
    }
    # Casos de la auditoría: miRNAs inexistentes → conjunto vacío.
    assert ke.get_targets_local("hsa-miR-1-5p", db) == set()
    assert ke.get_targets_local("hsa-miR-3-5p", db) == set()
    assert ke.get_targets_local("hsa-miR-9-3p", db) == set()
    assert ke.get_targets_local("hsa-miR-2-5p", db) == set()
    # Variante real .2: sigue encontrándose por límite de token.
    assert ke.get_targets_local("hsa-miR-873-5p", db) == {"GENE4"}


def test_analyze_validates_years_mode_and_panel_size():
    """v1.45 (auditoría externa, H-07/H-10): la UI enviaba el año
    calendario (2026) como `years` y el backend lo aceptaba sin validar
    (2026-2026=0 anulaba el filtro de PubMed); `mode` inválido caía
    silenciosamente a n-2. Ahora years es 1-50, mode es un literal y el
    panel tiene tope de 50 miRNAs."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        base = {"mirnas": ["hsa-miR-33a-5p"], "offline": True}
        assert client.post("/api/v1/analyze", json={**base, "years": 2026}).status_code == 422
        assert client.post("/api/v1/analyze", json={**base, "years": 0}).status_code == 422
        assert client.post("/api/v1/analyze", json={**base, "years": 51}).status_code == 422
        assert client.post("/api/v1/analyze", json={**base, "mode": "foo"}).status_code == 422
        assert client.post("/api/v1/analyze", json={"mirnas": ["hsa-miR-33a-5p"] * 51, "offline": True}).status_code == 422
        assert client.post("/api/v1/analyze", json={**base, "years": 10}).status_code == 200


def test_compare_validates_panels_and_types():
    """v1.45 (H-10): /compare devolvía 500 con `mirnas: [1, 2]` (int sin
    .strip) y aceptaba 300 paneles sin límite. Ahora los tipos y los
    límites se validan (422)."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        good = {"panels": [{"name": "a", "mirnas": ["hsa-miR-33a-5p"]}], "mode": "strict"}
        assert client.post("/api/v1/compare", json=good).status_code == 200
        assert client.post("/api/v1/compare", json={"panels": [{"name": "x", "mirnas": [1, 2]}]}).status_code == 422
        too_many = {"panels": [{"name": f"p{i}", "mirnas": ["hsa-miR-33a-5p"]} for i in range(21)]}
        assert client.post("/api/v1/compare", json=too_many).status_code == 422
        assert client.post("/api/v1/compare", json={"panels": []}).status_code == 422
        assert client.post("/api/v1/compare", json={"panels": [{"name": "x", "mirnas": ["hsa-miR-33a-5p"] * 31}]}).status_code == 422


def test_batch_rejects_oversize_file_and_too_many_rows():
    """v1.45 (H-10): /batch leía el archivo completo sin límite (200k
    filas → 22 s de CPU). Ahora 400 con >256 KB o >500 filas."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        big = b"panel,mirna\n" + b"p1,hsa-miR-33a-5p\n" * 300000
        assert len(big) > 256 * 1024
        r = client.post("/api/v1/batch", files={"file": ("big.csv", big, "text/csv")})
        assert r.status_code == 400
        many_rows = b"panel,mirna\n" + b"p1,hsa-miR-33a-5p\n" * 501
        r2 = client.post("/api/v1/batch", files={"file": ("rows.csv", many_rows, "text/csv")})
        assert r2.status_code == 400


def test_install_endpoint_requires_admin_token():
    """v1.45 (auditoría externa, H-01): /api/v1/updates/{source}/install
    permitía a cualquier cliente anónimo descargar y sobrescribir archivos
    del servidor. Sin KENRYU_ADMIN_TOKEN configurado → 403; con token,
    un token distinto → 403; con el token correcto y fuente desconocida →
    400 (llega a la lógica, sin descargar nada)."""
    import os
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        assert client.post("/api/v1/updates/foo/install").status_code == 403
        os.environ["KENRYU_ADMIN_TOKEN"] = "test-token"
        try:
            assert client.post("/api/v1/updates/foo/install").status_code == 403
            r = client.post("/api/v1/updates/foo/install",
                            headers={"X-Admin-Token": "test-token"})
            assert r.status_code == 400
        finally:
            del os.environ["KENRYU_ADMIN_TOKEN"]


def test_analyze_offline_marks_genes_beyond_research_cap_not_evaluated():
    """v1.46 (H-06): el panel de control (6 genes) supera el tope de 5
    investigaciones: el gen de menor prioridad sale como 'no evaluado'
    (sin puntaje), nunca 'baja', y la respuesta dice el tope."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    d = r.json()
    assert sorted(d["common_genes"]) == ["GMFB", "KLF12", "NFIA", "PDE3A", "SKI", "SOX6"]
    val = d["validation"]
    not_ev = [g for g, s in val["genes"].items() if s["level"] == "not_evaluated"]
    assert len(not_ev) == 1
    assert val["genes"][not_ev[0]]["score"] is None
    assert val["summary"]["n_not_evaluated"] == 1
    scope = d["research_scope"]
    assert scope["cap"] == 5
    assert len(scope["researched_genes"]) == 5
    assert scope["not_evaluated_genes"] == not_ev


def test_analyze_offline_cancer_panel_researches_prioritized_genes():
    """v1.46 (H-06): en el panel de cáncer (18 genes, strict) los 5
    investigados deben ser los de mayor prioridad por validación en
    miRTarBase (AGO4, MEIS1, PURB, RPS6KA3, WNK3), no los 5 primeros
    alfabéticos (AAK1, ACVR2A, AGO4, BNC2, BRWD3). Los 13 restantes son
    'no evaluado', nunca 'baja' por no haberse investigado. Excepción
    v1.52: sin ningún miRTarBase local (clon fresco) no hay validación
    que priorice, así que los 5 son los alfabéticos, con 25 ("baja")."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    d = r.json()
    assert len(d["common_genes"]) == 18
    val = d["validation"]
    researched = {g for g, s in val["genes"].items() if s["level"] != "not_evaluated"}
    not_ev = {g for g, s in val["genes"].items() if s["level"] == "not_evaluated"}
    # v1.48 (Ronda 3): con el índice oficial de miRTarBase, AAK1/CHD7/IL6R
    # tienen validación en 2 miRNAs del panel (antes 0 en el snapshot), así
    # que entran en el tope de 5 por delante de PURB/RPS6KA3/WNK3 (1 miRNA).
    # Con el GMT anterior, los 5 investigados eran los únicos validados.
    # v1.52: sin ningún miRTarBase local (clon fresco: el GMT ya no está
    # en el repo), todos los genes empatan a 0 validación y el tope de 5
    # cae a los 5 primeros alfabéticos, con 25 puntos ("baja", solo
    # TargetScan); los 13 restantes quedan sin evaluar.
    if _OFFICIAL.exists():
        expected = {"AAK1", "AGO4", "CHD7", "IL6R", "MEIS1"}
    elif _GMT.exists():
        expected = {"AGO4", "MEIS1", "PURB", "RPS6KA3", "WNK3"}
    else:
        expected = {"AAK1", "ACVR2A", "AGO4", "BNC2", "BRWD3"}
    assert researched == expected
    assert len(not_ev) == 13
    assert all(val["genes"][g]["score"] is None for g in not_ev)
    # Con validación local los 5 investigados son media/alta (low 0); sin
    # ningún miRTarBase local los 5 son "baja" (25, solo TargetScan).
    assert val["summary"]["low"] == (0 if (_OFFICIAL.exists() or _GMT.exists()) else 5)


def test_analyze_offline_n2_truncates_enrichment_by_priority():
    """v1.46 (H-06): en N-2 el panel de cáncer da 1.941 genes; el
    enriquecimiento usa los 300 priorizados (no los A–C alfabéticos) y
    la respuesta lo dice."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"],
            "years": 10, "mode": "n-2", "offline": True})
    assert r.status_code == 200
    d = r.json()
    assert len(d["common_genes"]) > 300
    assert d["enrichment_truncated"] is True
    assert "300" in (d.get("enrichment_note") or "")


class _FakeEnrichrResp:
    def __init__(self, status, payload):
        self.status_code = status
        self._payload = payload

    def json(self):
        return self._payload


class _FakeEnrichrClient:
    """Enrichr falso: un término significativo (p-adj 0,01) y uno que no
    lo es (0,0584) en cada biblioteca."""

    async def post(self, url, files=None, timeout=None):
        return _FakeEnrichrResp(200, {"userListId": 1})

    async def get(self, url, params=None, timeout=None):
        return _FakeEnrichrResp(200, {params["backgroundType"]: [
            [1, "Term A", 0.001, 5.0, 10.0, "GENE", 0.01],
            [2, "Term B", 0.002, 4.0, 9.0, "GENE", 0.0584],
        ]})


def test_enrich_genes_direct_ignores_stale_cache_and_filters_padj():
    """v1.46: la caché de Enrichr guardaba términos con el corte viejo
    (p-adj < 0,3); la clave no llevaba versión y se seguían sirviendo
    términos de 0,0584 como significativos. La clave versionada ignora
    esas entradas y el filtro < 0,05 se aplica a los resultados frescos."""
    import asyncio
    import kenryu_engine as ke
    genes = ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]
    old_key = ",".join(sorted(genes))
    ke.PERSISTENT_CACHE["enrichr"][old_key] = [
        {"Term": "Stale", "Source": "GO", "Pval": 0.0584, "RawTerm": "Stale"}]
    res = asyncio.run(ke.enrich_genes_direct(genes, _FakeEnrichrClient()))
    assert all(r["Pval"] < 0.05 for r in res)
    assert not any(r["Term"] == "Stale" for r in res)
    assert any(r["Term"] == "Term A" for r in res)
    assert ("enr-v2|" + old_key) in ke.PERSISTENT_CACHE["enrichr"]


def test_research_scope_note_states_priority_basis_and_snapshot():
    """v1.46 (ajuste pre-publicación): la nota del tope dice por qué se
    eligieron esos genes (nº de miRNAs del panel y validación en el
    snapshot local de miRTarBase) y con qué versión del snapshot."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    note = r.json()["research_scope"]["note"]
    assert "miRNAs del panel" in note
    assert "miRTarBase" in note
    assert "Harmonizome" in note


def test_analyze_offline_reports_source_per_mirna():
    """v1.46 (ajuste pre-publicación): la respuesta dice, por miRNA, qué
    fuente se usó — TargetScan local o miRTarBase como respaldo — porque
    un consenso strict puede mezclar tipos de evidencia distintos."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-33a-5p", "hsa-miR-106b-3p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    src = r.json()["source_by_mirna"]
    assert src["hsa-miR-33a-5p"].startswith("TargetScan")
    # hsa-miR-106b-3p has no TargetScan data: it falls back to miRTarBase
    # when a local copy exists (official index or GMT); with neither (fresh
    # clone since v1.52), it has no data anywhere and is reported as such.
    if _OFFICIAL.exists() or _GMT.exists():
        assert "miRTarBase" in src["hsa-miR-106b-3p"]
        assert "respaldo" in src["hsa-miR-106b-3p"]
    else:
        assert src["hsa-miR-106b-3p"] == "sin datos en ninguna fuente"
        assert "hsa-miR-106b-3p" in r.json()["missing_mirnas"]
    # El panel de control usa solo TargetScan.
    r2 = client.post("/api/v1/analyze", json={
        "mirnas": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
        "years": 10, "mode": "strict", "offline": True})
    assert r2.status_code == 200
    assert all(v.startswith("TargetScan") for v in r2.json()["source_by_mirna"].values())


# ── Ronda 2, Commit A (v1.47) — H-17, H-18, H-12, H-21 ─────────────────────
# Método: cada prueba se escribió primero y se confirmó que FALLA con el
# código de v1.46; luego se corrige. H-04 (paper.md / caso de estudio /
# core:265-266) NO está aquí: es decisión del usuario, pendiente.


def test_consensus_genes_empty_returns_empty_set():
    """H-17: consensus_genes([], ...) no debe lanzar TypeError."""
    assert ke.consensus_genes([], "strict") == set()
    assert ke.consensus_genes([], "n-1") == set()
    assert ke.consensus_genes([], "n-2") == set()


def test_consensus_dedupes_duplicate_mirnas():
    """H-17: los duplicados no deben sesgar el umbral de N-1/N-2.
    [21,21,21,155] en n-1 debe dar lo mismo que [21,155]."""
    s21 = {"A", "B", "C", "D"}
    s155 = {"C", "D", "E"}
    assert ke.consensus_genes([s21, s21, s21, s155], "n-1") == ke.consensus_genes([s21, s155], "n-1")


def test_analyze_reports_missing_mirnas():
    """H-17: un miRNA sin datos en ninguna fuente aparece en
    'missing_mirnas' en vez de descartarse en silencio."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-21-5p", "hsa-miR-9999-5p", "hsa-miR-155-5p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    d = r.json()
    assert "hsa-miR-9999-5p" in d["missing_mirnas"]
    assert "hsa-miR-21-5p" not in d["missing_mirnas"]
    assert "hsa-miR-155-5p" not in d["missing_mirnas"]


def test_analyze_single_mirna_warns_no_convergence():
    """H-17: un panel de 1 miRNA no tiene convergencia multi-miRNA;
    la respuesta debe advertirlo en 'single_mirna_warning'."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-21-5p"], "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    assert r.json()["single_mirna_warning"]
    assert "convergencia" in r.json()["single_mirna_warning"]
    # Un panel de varios miRNAs no lleva la advertencia.
    r2 = client.post("/api/v1/analyze", json={
        "mirnas": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
        "years": 10, "mode": "strict", "offline": True})
    assert r2.json()["single_mirna_warning"] is None


def test_analyze_invalid_mode_returns_422():
    """H-17 (regresión, ya en v1.45): mode fuera de {strict, n-1, n-2}
    devuelve 422, no cae en silencio a n-2."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-21-5p"], "years": 10, "mode": "foo", "offline": True})
    assert r.status_code == 422


def test_gene_details_has_no_fixed_confidence():
    """H-18: get_gene_details no asigna un 'confidence' fijo (85/88/92/95).
    El único puntaje que existe es el del dossier de validación."""
    import asyncio
    import httpx

    def _mock(request):
        return httpx.Response(200, json={"results": []})

    async def _run(gene):
        transport = httpx.MockTransport(_mock)
        async with httpx.AsyncClient(transport=transport) as client:
            return await ke.get_gene_details(gene, client)

    for gene in ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2", "GMFB", "KLF12", "AGO4"]:
        info = asyncio.run(_run(gene))
        assert "confidence" not in info, gene


def test_script_js_rigor_bar_uses_validation_score_only():
    """H-18: la barra 'Rigor Científico' usa solo el puntaje real del
    dossier; sin él muestra 'Sin evaluar', nunca un 85 fijo."""
    from pathlib import Path
    js = (Path(__file__).parent.parent / "static" / "script.js").read_text()
    assert "confidence?.score" not in js
    assert "|| 85" not in js
    assert "Sin evaluar" in js


def test_build_synthesis_no_unsupported_causal_claim():
    """H-12: sin clinical_evidence, la síntesis no afirma consecuencias
    causales ni cita con un índice [i+1] que no corresponde a nada."""
    syn = ke.build_synthesis(
        mirnas=["hsa-miR-21-5p", "hsa-miR-155-5p"],
        common=["GENE1"],
        gene_details={"GENE1": {"full_name": "Gene 1", "pathology": "Nodo regulador."}},
        enrichment=[],
        references=[],
        gene_research={})
    text = syn["academic"]
    assert "consecuencias directas en la estabilidad tisular" not in text
    assert "[1]" not in text
    assert "[2]" not in text


def test_build_synthesis_significance_label_only_below_005():
    """H-12: un término con p_adj 0,04 no es 'significancia funcional'
    ni 'equilibrio fisiológico'; la etiqueta dice significancia
    estadística y nada más."""
    syn = ke.build_synthesis(
        mirnas=["hsa-miR-21-5p"],
        common=["GENE1"],
        gene_details={"GENE1": {"full_name": "Gene 1", "pathology": "Nodo."}},
        enrichment=[{"Term": "Ruta X", "Source": "Reactome", "Pval": 0.04}],
        references=[],
        gene_research={})
    text = syn["functional"]
    assert "significancia funcional" not in text
    assert "mantener el equilibrio fisiológico" not in text
    assert "significancia estadística" in text


def test_gene_details_conclusion_is_hypothesis_not_intervention():
    """H-12: la conclusión genérica de un gen sin preset no afirma que
    el eje sea 'un punto crítico de intervención'; lo presenta como
    candidato para validación experimental."""
    import asyncio
    import httpx

    def _mock(request):
        # MyGene devuelve un hit para que se llegue a la conclusión
        # genérica; el resto (UniProt, traducción) devuelve vacío.
        if "mygene.info" in str(request.url):
            return httpx.Response(200, json={"hits": [
                {"name": "Argonaute 4", "summary": "involved in RNA silencing",
                 "pathway": {}, "go": {}}]})
        return httpx.Response(200, json={"results": []})

    async def _run(gene):
        transport = httpx.MockTransport(_mock)
        async with httpx.AsyncClient(transport=transport) as client:
            return await ke.get_gene_details(gene, client)

    info = asyncio.run(_run("AGO4"))
    conclusion = info.get("conclusion", "")
    assert "punto crítico de intervención" not in conclusion
    assert "candidato" in conclusion
    assert "validación experimental" in conclusion


def test_script_js_no_hardcoded_version():
    """H-21: la versión no está fija como 'v1.38' en el frontend; la
    exportación usa JSON.stringify para el YAML y lee la versión del
    backend."""
    from pathlib import Path
    js = (Path(__file__).parent.parent / "static" / "script.js").read_text()
    assert "v1.38" not in js
    assert "JSON.stringify(m)" in js
    assert "engine_version" in js


def test_analyze_response_carries_engine_version():
    """H-21: la respuesta de /analyze lleva la versión del motor para
    que la exportación no dependa de un literal fijo."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
            "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    assert r.json()["engine_version"] == ke.APP_VERSION


# ── Ronda 2, Commit B (v1.47) — H-08, H-14 ────────────────────────────────
# Método: cada prueba se escribió primero y se confirmó que FALLA con el
# código de v1.47-Commit-A; luego se corrige.


def test_cors_rejects_other_hf_spaces():
    """H-14: un Space ajeno (https://evil.hf.space) no recibe
    access-control-allow-origin ni allow-credentials."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.options("/api/v1/analyze", headers={
            "Origin": "https://evil.hf.space",
            "Access-Control-Request-Method": "POST"})
    assert r.headers.get("access-control-allow-origin") != "https://evil.hf.space"
    assert r.headers.get("access-control-allow-credentials") != "true"


def test_cors_allows_exact_space_origin():
    """H-14: el origen exacto del Space sí recibe allow-origin (sin
    credenciales: la app no usa cookies de sesión)."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.options("/api/v1/analyze", headers={
            "Origin": "https://kenryu007-bioinformatica.hf.space",
            "Access-Control-Request-Method": "POST"})
    assert r.headers.get("access-control-allow-origin") == "https://kenryu007-bioinformatica.hf.space"
    assert r.headers.get("access-control-allow-credentials") != "true"


def test_security_headers_present():
    """H-14: la respuesta lleva cabeceras de seguridad básicas
    (X-Content-Type-Options: nosniff como mínimo)."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.get("/health")
    assert r.headers.get("x-content-type-options") == "nosniff"


def test_script_js_escapes_external_data_in_innerhtml():
    """H-08: los datos externos y de usuario (gen buscado, genes del
    consenso, términos de enriquecimiento, texto de síntesis) se escapan
    antes de insertarse con innerHTML."""
    from pathlib import Path
    js = (Path(__file__).parent.parent / "static" / "script.js").read_text()
    assert "function escapeHtml" in js
    # Búsqueda inversa: el gen que escribe el usuario se refleja en la página.
    assert "escapeHtml(data.gene)" in js
    # Genes del consenso y tabla de validación.
    assert "escapeHtml(g)" in js
    # Términos de enriquecimiento (vienen de Enrichr).
    assert "escapeHtml" in js and js.count("escapeHtml(") >= 6


# ── Ronda 2, H-04 (v1.47) — el texto describe lo que el código hace ────────
# Decisión del usuario (2026-11-19, por recomendación externa): corregir el
# texto, no cambiar el código. El consenso se construye sobre TargetScan;
# miRTarBase es respaldo cuando TargetScan no tiene el miRNA, y verificación
# posterior en el dossier. La evidencia dual se mostrará por niveles
# (Ronda 3: miRTarBase oficial + dossier de tres niveles).


def test_paper_md_describes_real_algorithm():
    """H-04: paper.md no afirma evidencia dual obligatoria; describe la
    intersección TargetScan + verificación posterior."""
    from pathlib import Path
    text = (Path(__file__).parent.parent / "paper.md").read_text()
    assert "requires dual evidence" not in text
    assert "dual-evidence convergence" not in text
    assert "dual-evidence consensus" not in text
    assert "post-hoc" in text


def test_case_study_describes_real_strict_mode():
    """H-04: el caso de estudio no dice que strict exija validación en
    miRTarBase para cada miRNA, ni que cada gen lleve evidencia
    PubMed/OMIM (solo los investigados)."""
    from pathlib import Path
    text = (Path(__file__).parent.parent / "docs" / "case-study-cancer-panel.md").read_text()
    assert "validated in miRTarBase for every miRNA" not in text
    assert "Each gene in the JSON carries its PubMed" not in text
    assert "post-hoc" in text


def test_build_synthesis_intro_describes_real_algorithm():
    """H-04: la introducción de la síntesis no afirma intersección de
    TargetScan y miRTarBase."""
    syn = ke.build_synthesis(
        mirnas=["hsa-miR-21-5p", "hsa-miR-155-5p"],
        common=["GENE1"],
        gene_details={"GENE1": {"full_name": "Gene 1", "pathology": "Nodo."}},
        enrichment=[],
        references=[],
        gene_research={})
    assert "intersección de TargetScan y miRTarBase" not in syn["academic"]
    assert "TargetScan" in syn["academic"]
    assert "miRTarBase" in syn["academic"]


# ── Ronda 3 (v1.48) — interruptor "exigir evidencia experimental" ───────────
import pytest as _pytest

_THESIS = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
           "hsa-miR-106b-5p", "hsa-miR-33b-5p"]
_OFFICIAL = ke.BASE_DIR / "local_db" / "mirtarbase_official.json.gz"
# v1.52: the Harmonizome GMT is download-on-demand (not in the repo), so
# a fresh clone has NO local miRTarBase at all; tests that depend on the
# fallback branch on whichever local source exists (official, GMT, none).
_GMT = ke.BASE_DIR / "local_db" / "mirtarbase_gene_set_library.gmt.gz"




def test_analyze_default_keeps_all_thesis_genes():
    """v1.48: sin el interruptor (por defecto), el consenso no cambia:
    los 5 genes de la tesis, como siempre."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": _THESIS, "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    assert sorted(r.json()["common_genes"]) == [
        "ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]


def test_analyze_legacy_mode_warns_and_caps_at_media(tmp_path, monkeypatch):
    """v1.54: sin el índice oficial (solo el GMT legacy), la respuesta
    dice mirtarbase_mode="legacy" y trae un aviso visible; la evidencia
    de ABCA1 es "legacy" (sin nivel verificado), así que nunca llega a
    "alta" aunque el puntaje lo permitiría (ver test_legacy_capped_at_media)."""
    import gzip
    gmt = tmp_path / "mirtarbase_gene_set_library.gmt.gz"
    gmt.write_bytes(gzip.compress(
        b"hsa-miR-33a-5p\tMIRT000001\tABCA1\tSNTB2\n"
        b"hsa-miR-33b-5p\tMIRT000002\tABCA1\n"))
    monkeypatch.setattr(ke, "MIRTARBASE_OFFICIAL", str(tmp_path / "nope.json.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_GMT", str(gmt))
    monkeypatch.setattr(ke, "MIRTARBASE_DB", {})
    monkeypatch.setattr(ke, "MIRTARBASE_FORMAT", None)
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": _THESIS, "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    d = r.json()
    assert d["mirtarbase_mode"] == "legacy"
    assert d["warnings"], "el modo legacy debe traer un aviso visible"
    assert "legacy" in d["warnings"][0].lower()
    assert "oficial" in d["warnings"][0].lower()
    assert d["validation"]["genes"]["ABCA1"]["evidence_level"] == "legacy"
    assert d["validation"]["summary"]["high"] == 0




# ── v1.49 — fondo del enriquecimiento: universo predicho ───────────────────
# Método: prueba escrita primero; falló con el fondo de la unión de GMTs
# (el panel de control devolvía 1 término offline, artefacto del fondo
# equivocado: las dianas del panel ya están enriquecidas en esa ruta,
# así que el fondo correcto la descuenta).


def test_analyze_offline_enrichment_uses_predicted_universe():
    """v1.49: el enriquecimiento local usa como fondo la unión de las
    dianas predichas del panel, no la unión de cada GMT. El panel de
    control lo demuestra: con la unión de GMTs daba 1 término offline;
    con el universo predicho da 0, y la respuesta coincide con el
    cálculo directo con ese fondo."""
    from fastapi.testclient import TestClient
    from kenryu_core import local_enrichment

    mirnas = ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"]
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": mirnas, "years": 10, "mode": "strict", "offline": True})
    assert r.status_code == 200
    body = r.json()
    api_terms = body["enrichment"]
    ke.load_local_data()
    ke.load_gmt_libraries()
    predicted = set()
    for m in mirnas:
        predicted |= ke.get_targets_local(m)
    assert len(predicted) > 1000  # universo predicho real, no vacío
    direct_p = local_enrichment(
        body["common_genes"], ke.GMT_LIBRARIES, background=predicted)
    direct_u = local_enrichment(body["common_genes"], ke.GMT_LIBRARIES)
    # La API usa el universo predicho, no la unión de GMTs.
    assert [(t["Term"], t["Source"]) for t in api_terms] == [
        (t["Term"], t["Source"]) for t in direct_p]
    for a, d in zip(api_terms, direct_p):
        assert a["Pval"] == pytest.approx(d["Pval"])
    # Y el fondo cambia el resultado de verdad: el término GO que daba
    # la unión de GMTs era un artefacto del fondo equivocado.
    assert len(direct_u) == 1 and direct_u[0]["Source"] == "GO"
    assert direct_p != direct_u


def test_analyze_offline_enrichment_uses_targetscan_universe():
    """v1.54: the offline/local enrichment background is the full
    TargetScan universe (all genes with any predicted site), not the
    panel's own union. Discriminating case: the cancer panel in N-1
    mode (188 genes) gives 0 local terms with the panel-union
    background (v1.49) and 16 with the TargetScan universe. The
    response must match the direct per-library computation with the
    intersected background."""
    import kenryu_core as kc
    from fastapi.testclient import TestClient
    cancer = ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"]
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": cancer, "years": 10, "mode": "n-1", "offline": True})
    assert r.status_code == 200
    d = r.json()
    genes = set(d["common_genes"][:300])
    ts_universe = set().union(*[set(g) for g in ke.TARGETSCAN_DB.values()])
    expected = []
    for label, lib in ke.GMT_LIBRARIES.items():
        bg = ts_universe & set().union(*lib.values())
        for row in kc.fisher_enrichment(genes, lib, background=bg):
            expected.append({"Term": row["Term"], "Source": label,
                             "Pval": row["Pval"]})
    expected.sort(key=lambda x: x["Pval"])
    expected = expected[:24]
    got = [(e["Term"], e["Source"], e["Pval"]) for e in d["enrichment"]]
    assert got == [(e["Term"], e["Source"], e["Pval"]) for e in expected]
    assert len(got) > 0, "the TargetScan background must return terms here"


def test_analyze_offline_enrichment_has_specificity():
    """v1.54: cada término de enriquecimiento trae 'specificity' — la
    fracción de paneles aleatorios del mismo tamaño en los que el
    término aparece por azar. El panel de cáncer N-1 (188 genes) da 16
    términos locales; todos llevan el campo, entre 0 y 1."""
    from fastapi.testclient import TestClient
    cancer = ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"]
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": cancer, "years": 10, "mode": "n-1", "offline": True})
    assert r.status_code == 200
    terms = r.json()["enrichment"]
    assert len(terms) > 0
    for t in terms:
        assert "specificity" in t
        assert t["specificity"] is None or 0.0 <= t["specificity"] <= 1.0


def test_annotate_specificity_uses_raw_term_for_translated_terms():
    """v1.54: the live path translates Term to Spanish; the English name
    survives in RawTerm, and that is what the GMT lookup must use.
    The 3 Reactome terms of the cancer case study are specific
    (0.0, 0.0, 0.01): they essentially never appear in random panels
    of the same size."""
    import kenryu_engine as ke
    # load_local_data() returns the dict but does not assign
    # TARGETSCAN_DB (the app startup does that); assign it here.
    ke.TARGETSCAN_DB = ke.load_local_data()
    ke.load_gmt_libraries()
    cancer = ["AAK1", "ACVR2A", "AGO4", "BNC2", "BRWD3", "CHD7", "DCUN1D3",
              "FOXP1", "IL6R", "KPNA4", "MEIS1", "PURB", "RDX", "RPS6KA3",
              "SESN3", "UBN2", "WNK3", "XKR4"]
    terms = [
        {"Term": "Biología del desarrollo", "RawTerm": "Developmental Biology",
         "Source": "Reactome", "Pval": 0.0176},
        {"Term": "Vía de reciclaje de L1", "RawTerm": "Recycling Pathway of L1",
         "Source": "Reactome", "Pval": 0.0424},
        {"Term": "Señalización mediante TGFBR3", "RawTerm": "Signaling by TGFBR3",
         "Source": "Reactome", "Pval": 0.0424},
    ]
    ke.annotate_enrichment_specificity(terms, cancer)
    assert [t["specificity"] for t in terms] == [0.0, 0.0, 0.01]


def test_health_reports_commit_hash_from_build_info(tmp_path, monkeypatch):
    # Unreleased: /health identifica el código exacto con el hash corto
    # de GitHub, leído de local_db/build_info.json (generado al sincronizar).
    info = tmp_path / "build_info.json"
    info.write_text(json.dumps({"commit": "abc1234", "version": "1.54"}))
    monkeypatch.setattr(ke, "BUILD_INFO_PATH", info)
    with TestClient(ke.app) as client:
        r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["commit"] == "abc1234"


def test_health_reports_unknown_commit_without_build_info(tmp_path, monkeypatch):
    # Sin build_info.json no se inventa nada: "unknown".
    monkeypatch.setattr(ke, "BUILD_INFO_PATH", tmp_path / "missing.json")
    with TestClient(ke.app) as client:
        r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["commit"] == "unknown"


def test_analyze_includes_commit_hash(tmp_path, monkeypatch):
    info = tmp_path / "build_info.json"
    info.write_text(json.dumps({"commit": "deadbee"}))
    monkeypatch.setattr(ke, "BUILD_INFO_PATH", info)
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-33a-5p"], "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    assert r.json()["commit"] == "deadbee"


def test_analyze_unknown_commit_without_build_info(tmp_path, monkeypatch):
    monkeypatch.setattr(ke, "BUILD_INFO_PATH", tmp_path / "missing.json")
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-33a-5p"], "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    assert r.json()["commit"] == "unknown"
