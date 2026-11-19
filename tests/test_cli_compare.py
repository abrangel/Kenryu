"""Tests for the CLI and the compare/batch endpoints (local data only)."""

import json
from pathlib import Path

import kenryu_cli as cli


def test_read_panel(tmp_path):
    p = tmp_path / "panel.txt"
    p.write_text("# comment\n\nhsa-miR-21-5p, hsa-miR-155-5p\nhsa-let-7a-5p\n")
    assert cli.read_panel(p) == ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"]


def test_cli_offline(tmp_path):
    p = tmp_path / "panel.txt"
    p.write_text("hsa-miR-375\nhsa-miR-29a-3p\nhsa-miR-21-5p\n")
    out = tmp_path / "out"
    rc = cli.main([str(p), "--offline", "--out", str(out)])
    assert rc == 0
    data = json.loads((out / "analysis.json").read_text())
    assert data["common_genes"] == ["GMFB", "KLF12", "NFIA", "PDE3A", "SKI", "SOX6"]
    assert (out / "report.md").read_text().startswith("# KENRYU")
    assert (out / "venn.png").stat().st_size > 1000
    assert (out / "heatmap.png").stat().st_size > 1000


def test_compare_endpoint():
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    lipid = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
             "hsa-miR-106b-5p", "hsa-miR-33b-5p"]
    cancer = ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"]
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/compare",
                        json={"panels": [{"name": "lipid", "mirnas": lipid},
                                         {"name": "cancer", "mirnas": cancer}]})
    assert r.status_code == 200
    d = r.json()
    by_name = {p["name"]: p for p in d["panels"]}
    assert by_name["lipid"]["genes"] == ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]
    assert len(by_name["cancer"]["genes"]) == 18
    assert d["n_union"] == len(set(by_name["lipid"]["genes"]) | set(by_name["cancer"]["genes"]))
    assert d["pairwise"][0]["a"] == "lipid"


def test_batch_endpoint_csv_and_text():
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    csv_body = ("panel,mirna\n"
                "lipid,hsa-miR-33a-5p\nlipid,hsa-miR-144-3p\n"
                "cancer,hsa-miR-21-5p\ncancer,hsa-miR-155-5p\ncancer,hsa-let-7a-5p\n")
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/batch",
                        files={"file": ("panels.csv", csv_body, "text/csv")})
    assert r.status_code == 200
    assert {p["name"] for p in r.json()["panels"]} == {"lipid", "cancer"}

    text_body = "lipid: hsa-miR-33a-5p, hsa-miR-144-3p\ncancer: hsa-miR-21-5p, hsa-miR-155-5p, hsa-let-7a-5p\n"
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/batch",
                        files={"file": ("panels.txt", text_body, "text/plain")})
    assert r.status_code == 200
    assert len(r.json()["panels"]) == 2


def test_disgenet_mocked():
    import asyncio
    import httpx
    import kenryu_engine as ke

    def handler(request):
        assert request.headers["Authorization"] == "KEY123"
        return httpx.Response(200, json=[
            {"disease_name": "Tangier disease", "score": 0.9},
            {"disease_name": "Atherosclerosis", "score": 0.7},
        ])

    async def run():
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        out = await ke.fetch_disgenet("ABCA1", client, api_key="KEY123")
        assert out["diseases"][0]["name"] == "Tangier disease"
        assert len(out["diseases"]) == 2
        # Without a key it is skipped, never faked.
        assert (await ke.fetch_disgenet("ABCA1", client, api_key=None)) == {"diseases": []}
        await client.aclose()

    asyncio.run(run())


def test_mouse_db_loads_and_analyze_offline():
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    mouse = ke.load_mouse_data()
    assert len(mouse) > 100
    # TargetScanMouse names the family mir-21-5p (mouse miR-21a lives there).
    assert "mir-21-5p" in mouse
    assert len(mouse["mir-21-5p"]) > 50
    # Mouse gene symbols are title-case (Acat1, Ago2), not human caps.
    assert "Acat1" in mouse["mir-21-5p"] and "Ago2" in mouse["mir-21-5p"]

    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["mmu-miR-21-5p", "mmu-miR-155-5p", "mmu-let-7a-5p"],
            "mode": "strict", "years": 10, "offline": True, "species": "mouse"})
    assert r.status_code == 200
    d = r.json()
    assert d["species"] == "mouse"
    assert d["found_mirnas"] == ["mmu-miR-21-5p", "mmu-miR-155-5p", "mmu-let-7a-5p"]
    assert isinstance(d["common_genes"], list)


def test_validation_ui_hooks():
    """The web page must have a place to SHOW the validation dossier
    (table + heatmap) and the gene panel must use the real score."""
    import re
    html = (Path(__file__).resolve().parent.parent / "static" / "index.html").read_text()
    js = (Path(__file__).resolve().parent.parent / "static" / "script.js").read_text()
    assert "validation-container" in html
    assert "heatmap-container" in html
    assert "validation-badge" in html
    assert "data.validation" in js
    assert "heatmap_plot" in js
    assert "validation?.genes?.[gene]?.score" in js
    # No hardcoded gene scores left in the panel.
    assert "gene === 'ABCA1'" not in js


def test_reverse_lookup_ui_and_explanations():
    """The page must EXPLAIN every graphic in plain language and offer the
    reverse lookup (gene -> miRNAs) as a visible control."""
    import re
    html = (Path(__file__).resolve().parent.parent / "static" / "index.html").read_text()
    js = (Path(__file__).resolve().parent.parent / "static" / "script.js").read_text()
    assert "Búsqueda inversa" in html
    assert 'id="reverse-gene"' in html
    assert "reverseLookup" in js
    assert "/api/v1/reverse" in js
    assert "card-sub" in html
    assert "Cada círculo es un miARN" in html
    assert "Cada punto es un gen real" in html
    assert "Cada fila es un gen" in html
    assert "STRING-DB" in html
