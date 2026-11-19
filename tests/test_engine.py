"""Offline test suite for the KENRYU engine (kenryu_engine.py).

Every test runs without network access, using only the data files
shipped in the repository (data/hsa-miR-*.txt). Nothing here modifies
the engine: tests assert the behaviour the thesis documents.

Run:  pytest -v
"""
import base64

import pytest

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
                 "desc": "miR-33 suprime ABCA1.", "pmid": "30212345"},
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
    assert body["version"] == "1.41"
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
