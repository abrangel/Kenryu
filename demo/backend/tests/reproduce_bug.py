import pytest
from fastapi.testclient import TestClient
from demo.backend.main import app

client = TestClient(app)

def test_analyze_single_mirna_index_error():
    # Este test debería fallar con IndexError en la versión actual si len(mirnas) == 1
    response = client.post("/api/v1/analyze", json={
        "mirnas": ["hsa-miR-33a-5p"],
        "years": 10
    })
    # Si hay un error 500 y un IndexError en el servidor, esto lo atrapará
    assert response.status_code == 200
