"""Gene-symbol normalisation (HGNC) — v1.59, policy tightened 2026-11-21.

The audit (scripts/audit_gene_symbols.py) found 2,496 of 13,035
TargetScan genes in no local enrichment library, and HGNC recovers
487 of them with a unique current symbol in the libraries (+11 only
via KEGG). The map (local_db/hgnc_symbol_map.json, 498 entries) is
applied at load time to the TargetScan sets and to the miRTarBase
index keys, so both layers speak current symbols. Ambiguous genes
are never mapped; unmapped genes stay as they are. GUARD: a gene
that is itself an Approved HGNC symbol is never mapped — it denotes
that gene, not the prev/alias target (TargetScan's BCAM is
ENSG00000187244 = HGNC BCAM, not BCAT2). 1,504 current genes the
libraries lack stay untouched.
"""

import json
from pathlib import Path

import kenryu_engine as ke

REPO = Path(__file__).resolve().parent.parent
MAP_PATH = REPO / "local_db" / "hgnc_symbol_map.json"


def _map():
    return json.loads(MAP_PATH.read_text(encoding="utf-8"))


def test_symbol_map_file():
    """The map is exactly the audited set: 498 entries (487 + 11),
    3 ambiguous, AARS -> AARS1. POLICY (2026-11-21): the previous
    symbol is preferred over an alias — an ambiguous previous symbol
    is discarded, never rescued by an alias. That recovers 5 genes
    the merged logic left ambiguous (AGPAT9, LOR, MPP6, NOV, SEPT2)
    and leaves 3 discarded (B3GNT1, MUM1, QARS). GUARD: a gene that
    is itself an Approved symbol is never mapped — 1,504 current
    genes the libraries lack stay untouched (e.g. BCAM stays BCAM,
    not BCAT2)."""
    doc = _map()
    assert len(doc["mapping"]) == 498
    assert doc["counts"]["in_three_libraries"] == 487
    assert doc["counts"]["kegg_only"] == 11
    assert len(doc["ambiguous_not_mapped"]) == 3
    assert doc["counts"]["approved_current_skipped"] == 1504
    assert doc["mapping"]["AARS"] == "AARS1"
    assert doc["mapping"]["AGPAT9"] == "GPAT3"
    assert doc["mapping"]["LOR"] == "LORICRIN"
    assert doc["mapping"]["MPP6"] == "PALS2"
    assert doc["mapping"]["NOV"] == "CCN3"
    assert doc["mapping"]["SEPT2"] == "SEPTIN2"
    assert "BCAM" not in doc["mapping"]
    assert "FADS6" not in doc["mapping"]
    assert "CHD5" not in doc["mapping"]
    assert sorted(doc["ambiguous_not_mapped"]) == [
        "B3GNT1", "MUM1", "QARS"]
    assert doc["hgnc_licence"].startswith("CC0")
    assert doc["hgnc_downloaded"] == "2026-11-21"


def test_targetscan_loader_maps_symbols():
    """AARS is a predicted target of mir-145-5p under its old symbol;
    after loading, the set carries AARS1, not AARS. The universe
    shrinks from 13,035 to 13,016 (15 merges into a symbol already
    present under its current name + 4 collision groups). No gene
    is lost: every dropped count is two symbols of one gene."""
    db = ke.load_local_data()
    assert "AARS1" in db["mir-145-5p"]
    assert "AARS" not in db["mir-145-5p"]
    universe = set().union(*db.values())
    assert len(universe) == 13016
    assert "AARS" not in universe
    assert "ACN9" not in universe and "SDHAF3" in universe
    assert "BCAM" in universe and "BCAT2" in universe


def test_mirtarbase_loader_maps_symbols():
    """The validated index carries AARS for hsa-miR-124-3p (and four
    more miRNAs) under its old symbol; after loading, the evidence
    key is AARS1. The fixture stands in for the real index, which is
    gitignored and absent in CI; the real index is only asserted when
    it is actually loaded. A second fixture checks the merge: two old
    symbols falling on one current symbol keep the best level and the
    union of PMIDs."""
    fixture = {
        "hsa-miR-124-3p": {
            "AARS": {"level": "weak", "pmids": [12345678],
                     "experiments": ["CLIP-seq"]},
        },
        "hsa-miR-21-5p": {
            "FAM21A": {"level": "weak", "pmids": [111],
                       "experiments": ["HITS-CLIP"]},
            "FAM21B": {"level": "strong", "pmids": [222],
                       "experiments": ["Luciferase"]},
        },
    }
    mapped = ke._apply_symbol_map_to_mirtarbase(fixture)
    ev = ke.get_validated_evidence("hsa-miR-124-3p", mapped)
    assert "AARS1" in ev
    assert "AARS" not in ev
    ev21 = ke.get_validated_evidence("hsa-miR-21-5p", mapped)
    assert "WASHC2A" in ev21
    assert "FAM21A" not in ev21 and "FAM21B" not in ev21
    assert ev21["WASHC2A"]["level"] == "strong"
    assert sorted(ev21["WASHC2A"]["pmids"]) == [111, 222]
    # Real index, only when it is loaded (gitignored, absent in CI).
    ke.load_mirtarbase_local()
    if ke.MIRTARBASE_DB:
        ev_real = ke.get_validated_evidence("hsa-miR-124-3p",
                                            ke.MIRTARBASE_DB)
        assert "AARS1" in ev_real
        assert "AARS" not in ev_real


def test_ambiguous_and_current_symbols_untouched():
    """Ambiguous genes are never mapped, and current symbols pass
    through unchanged — including Approved symbols that share a
    string with another gene's prev/alias (CHD5 is an Approved
    symbol in its own right; BCAM and FADS6 likewise)."""
    db = ke.load_local_data()
    assert ke.map_gene_symbol("B3GNT1") == "B3GNT1"
    assert ke.map_gene_symbol("CHD5") == "CHD5"
    assert ke.map_gene_symbol("MUM1") == "MUM1"
    assert ke.map_gene_symbol("QARS") == "QARS"
    assert ke.map_gene_symbol("BCAM") == "BCAM"
    assert ke.map_gene_symbol("FADS6") == "FADS6"
    assert ke.map_gene_symbol("ABCA1") == "ABCA1"
    assert ke.map_gene_symbol("AARS") == "AARS1"
    assert ke.map_gene_symbol("AGPAT9") == "GPAT3"
    # Ambiguous genes stay in the universe under their own symbol:
    # B3GNT1 is a TargetScan gene and is NOT replaced by B3GNT2/B4GAT1.
    universe = set().union(*db.values())
    assert "B3GNT1" in universe
