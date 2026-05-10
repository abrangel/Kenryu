import pytest
from demo.backend.main import clean_mirna_name

def test_clean_mirna_name_basic():
    # Test normalización básica
    assert clean_mirna_name("mir-33a") == "hsa-miR-33a"
    assert clean_mirna_name("hsa-mir-33a-5p") == "hsa-miR-33a-5p"

def test_clean_mirna_name_edge_cases():
    # Test casos complejos (guiones especiales, espacios)
    assert clean_mirna_name(" mir‑33a ") == "hsa-miR-33a"
    assert clean_mirna_name("LET-7") == "hsa-miR-7" # Mapeo simplificado

def test_pipeline_request_validation():
    # Aquí podríamos añadir tests de integración con la API si fuera necesario
    pass
