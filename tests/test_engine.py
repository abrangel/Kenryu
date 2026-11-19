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

def test_create_visuals_returns_three_pngs():
    sets = [{"ABCA1", "KPNA3"}, {"ABCA1", "SCN1A"}, {"ABCA1", "SNTB2"}]
    venn, volcano, ppi = ke.create_visuals(sets, PANEL[:3], ["ABCA1"])
    for img in (venn, volcano, ppi):
        raw = base64.b64decode(img)
        assert raw[:8] == b"\x89PNG\r\n\x1a\n"  # PNG magic bytes


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
    assert body["version"] == "1.39"
    assert body["db_families"] >= 5  # fallback: the 5 panel .txt files
    assert body["raw_files"] == 5
