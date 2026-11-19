"""Rate limiting on the public API (auditoría externa, H-10).

H-10 found no rate limiting on any endpoint: a single client could
flood /api/v1/analyze (each call burns 30-90 s of live NCBI/Enrichr
calls) and exhaust the Space's quota for those APIs. The fix is a
per-IP sliding-window middleware on /api/v1/*, with /health exempt
(the Hugging Face health checks must never fail).

These tests FAIL on v1.49 (no limiter exists yet) and pass once the
middleware lands. The limiter is OFF by default inside the test
suite (conftest sets KENRYU_RATE_LIMIT=0) so the other 160+ tests
never trip it; each test here enables it explicitly and resets its
state.
"""
import pytest
from fastapi.testclient import TestClient

import kenryu_engine as ke


@pytest.fixture
def limited(monkeypatch):
    """Enable the limiter with a tiny budget (3/min) for one test."""
    monkeypatch.setattr(ke, "RATE_LIMIT_ENABLED", True)
    monkeypatch.setattr(ke.RATE_LIMITER, "per_minute", 3)
    ke.RATE_LIMITER.reset()
    yield
    ke.RATE_LIMITER.reset()


def test_rate_limit_blocks_after_limit(limited):
    """v1.50 (H-10): the 4th request inside the window gets 429."""
    with TestClient(ke.app) as client:
        for _ in range(3):
            assert client.get("/api/v1/reverse?gene=ABCA1").status_code == 200
        r = client.get("/api/v1/reverse?gene=ABCA1")
        assert r.status_code == 429
        assert "Retry-After" in r.headers


def test_rate_limit_exempts_health(limited):
    """v1.50 (H-10): /health is never limited (HF health checks)."""
    with TestClient(ke.app) as client:
        for _ in range(5):
            assert client.get("/health").status_code == 200


def test_rate_limit_per_ip_via_x_forwarded_for(limited):
    """v1.50 (H-10): behind the HF proxy every client arrives with the
    proxy's IP, so the bucket key is X-Forwarded-For (first hop).
    Two different IPs get independent budgets."""
    with TestClient(ke.app) as client:
        a = {"X-Forwarded-For": "1.2.3.4"}
        b = {"X-Forwarded-For": "5.6.7.8"}
        assert client.get("/api/v1/reverse?gene=ABCA1", headers=a).status_code == 200
        assert client.get("/api/v1/reverse?gene=ABCA1", headers=a).status_code == 200
        assert client.get("/api/v1/reverse?gene=ABCA1", headers=a).status_code == 200
        assert client.get("/api/v1/reverse?gene=ABCA1", headers=a).status_code == 429
        # The other IP still has its full budget.
        assert client.get("/api/v1/reverse?gene=ABCA1", headers=b).status_code == 200


def test_rate_limit_off_in_test_suite():
    """The suite runs with KENRYU_RATE_LIMIT=0 (conftest): a burst of
    API calls must never 429, or the whole suite would be flaky."""
    assert ke.RATE_LIMIT_ENABLED is False
    with TestClient(ke.app) as client:
        for _ in range(10):
            assert client.get("/api/v1/reverse?gene=ABCA1").status_code == 200
