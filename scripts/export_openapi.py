#!/usr/bin/env python3
"""Export the live OpenAPI schema to docs/openapi.json.

Run from the repo root:  python scripts/export_openapi.py
The export must stay in sync with the app — tests/test_openapi.py
fails if it drifts. Run this after any endpoint change and commit
the result with the change.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import kenryu_engine as ke  # noqa: E402
from fastapi.testclient import TestClient  # noqa: E402

OUT = ROOT / "docs" / "openapi.json"

with TestClient(ke.app) as client:
    schema = client.get("/openapi.json").json()

OUT.write_text(json.dumps(schema, indent=2, ensure_ascii=False) + "\n",
               encoding="utf-8")
print(f"openapi.json: {schema['openapi']}, "
      f"version {schema['info']['version']}, "
      f"{len(schema['paths'])} paths -> {OUT}")
