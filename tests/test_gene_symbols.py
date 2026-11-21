"""Gene-symbol normalisation (HGNC) — v1.59.

The audit (scripts/audit_gene_symbols.py) found 2,496 of 13,035
TargetScan genes in no local enrichment library, and HGNC recovers
527 of them with a unique current symbol in the libraries (+11 only
via KEGG). The map (local_db/hgnc_symbol_map.json, 538 entries) is
applied at load time to the TargetScan sets and to the miRTarBase
index keys, so both layers speak current symbols. Ambiguous genes
are never mapped; unmapped genes stay as they are.
"""

import json
from pathlib import Path

import kenryu_engine as ke

REPO = Path(__file__).resolve().parent.parent
MAP_PATH = REPO / "local_db" / "hgnc_symbol_map.json"


def _map():
    return json.loads(MAP_PATH.read_text(encoding="utf-8"))


def test_symbol_map_file():
    """The map is exactly the audited set: 538 entries (527 + 11),
    9 ambiguous, AARS -> AARS1."""
    doc = _map()
    assert len(doc["mapping"]) == 538
    assert doc["counts"]["in_three_libraries"] == 527
    assert doc["counts"]["kegg_only"] == 11
    assert len(doc["ambiguous_not_mapped"]) == 9
    assert doc["mapping"]["AARS"] == "AARS1"
    assert "AGPAT9" not in doc["mapping"]
    assert "AGPAT9" in doc["ambiguous_not_mapped"]
    assert doc["hgnc_licence"].startswith("CC0")


def test_targetscan_loader_maps_symbols():
    """AARS is a predicted target of mir-145-5p under its old symbol;
    after loading, the set carries AARS1, not AARS. The universe
    shrinks from 13,035 to 12,985 (5 merges + 45 targets already
    present under their current name)."""
    db = ke.load_local_data()
    assert "AARS1" in db["mir-145-5p"]
    assert "AARS" not in db["mir-145-5p"]
    universe = set().union(*db.values())
    assert len(universe) == 12985
    assert "AARS" not in universe
    assert "ACN9" not in universe and "SDHAF3" in universe


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
    through unchanged."""
    db = ke.load_local_data()
    assert ke.map_gene_symbol("AGPAT9") == "AGPAT9"
    assert ke.map_gene_symbol("ABCA1") == "ABCA1"
    assert ke.map_gene_symbol("AARS") == "AARS1"
    # Ambiguous genes stay in the universe under their own symbol:
    # AGPAT9 is a TargetScan gene and is NOT replaced by GPAT3/LPCAT1.
    universe = set().union(*db.values())
    assert "AGPAT9" in universe
