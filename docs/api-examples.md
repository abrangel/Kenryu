# API examples: 3 calls that run

Base URL: the public Space `https://kenryu007-bioinformatica.hf.space`
or a local run at `http://localhost:7860`. The versioned schema is
`docs/openapi.json` (also live at `/openapi.json`; interactive docs at
`/docs`). All three payloads live in `docs/examples/` and are run by
`tests/test_openapi.py` — they are verified, not just written.

`offline: true` means local data only: no calls to Enrichr, PubMed or
any other service. The gene order in `common_genes` depends on the
validation data available: with the official miRTarBase index it is
prioritized by evidence (ABCA1 first); without it, alphabetical.

## 1. Analyze a panel (curl)

```bash
curl -s -X POST https://kenryu007-bioinformatica.hf.space/api/v1/analyze \
  -H 'Content-Type: application/json' \
  -d @docs/examples/analyze-lipid.json
```

Payload (`docs/examples/analyze-lipid.json`): the five-miRNA lipid
panel, Strict mode, offline. Abridged response (public Space, v1.57):

```json
{
  "common_genes": ["ABCA1", "SNTB2", "TSC22D2", "KPNA3", "SCN1A"],
  "convergence": {
    "genes_observed": 5, "k_effective": 4,
    "p_by_tolerance": {"0.10": 0.126, "0.25": 0.123, "0.50": 0.097},
    "seed": 2547548954
  },
  "execution_manifest": {
    "engine_version": "1.57", "commit": "1d720ee",
    "parameters": {"mirnas": ["hsa-miR-33a-5p", "hsa-miR-144-3p",
      "hsa-miR-758-3p", "hsa-miR-106b-5p", "hsa-miR-33b-5p"],
      "mode": "strict", "offline": true},
    "seeds": {"convergence": 2547548954},
    "reproduce": "curl -s -X POST http://localhost:7860/api/v1/analyze ..."
  }
}
```

## 2. Analyze a panel (Python) and read the manifest

```python
import json, requests

BASE = "https://kenryu007-bioinformatica.hf.space"
payload = json.load(open("docs/examples/analyze-cancer-n1.json"))
r = requests.post(f"{BASE}/api/v1/analyze", json=payload, timeout=120)
r.raise_for_status()
d = r.json()
print(len(d["common_genes"]), "genes")          # 188 (n-1 mode, cancer panel)
print(d["convergence"]["p_by_tolerance"])       # {'0.10': 0.505, '0.25': 0.501, '0.50': 0.254}
print(d["execution_manifest"]["reproduce"])     # the curl command that repeats this run
```

Payload (`docs/examples/analyze-cancer-n1.json`): the three-miRNA
cancer panel, N-1 mode, offline. N-1 admits genes missed by one
miRNA, so the panel grows from 18 genes (Strict) to 188 — and the
convergence p rises to ≈0.5, which is what the test is for: the count
is not unusual. The manifest's `reproduce` command is the exact call
to repeat the run.

## 3. Compare two panels (curl)

```bash
curl -s -X POST https://kenryu007-bioinformatica.hf.space/api/v1/compare \
  -H 'Content-Type: application/json' \
  -d @docs/examples/compare-lipid-cancer.json
```

Payload (`docs/examples/compare-lipid-cancer.json`): the lipid panel
and the cancer panel, both Strict. Abridged response:

```json
{
  "panels": [
    {"name": "lipid", "genes": ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]},
    {"name": "cancer", "genes": ["AGO4", "MEIS1", "... 18 genes ..."]}
  ],
  "n_union": 23
}
```

The two panels share no gene (`n_union` 23 = 5 + 18): the lipid and
cancer convergences are disjoint sets.
