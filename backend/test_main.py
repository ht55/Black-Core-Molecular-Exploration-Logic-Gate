# test_main.py

import pytest
from fastapi.testclient import TestClient
from main import app
from graph_engine import (
    extract_json,
    calculate_metrics,
    design_router,
    physics_router,
)

client = TestClient(app)


# ========================================
# Existing Test
# ========================================

def test_catchall_returns_something():
    response = client.get("/")
    assert response.status_code in [200, 404]


# ========================================
# Phase 1-A: extract_json
# ========================================

def test_extract_json_normal():
    """Valid JSON is parsed correctly."""
    result = extract_json('{"candidates": [{"smiles": "C"}]}')
    assert result is not None
    assert result["candidates"][0]["smiles"] == "C"

def test_extract_json_broken_recovers():
    """Broken JSON with missing closing tags is recovered."""
    broken = '{"candidates": [{"smiles": "C", "description": "test"'
    result = extract_json(broken)
    assert result is not None

def test_extract_json_empty_returns_none():
    """Empty string returns None."""
    assert extract_json("") is None

def test_extract_json_no_braces_returns_none():
    """Plain text with no JSON structure returns None."""
    assert extract_json("plain text, no json here") is None

def test_extract_json_with_markdown_fences():
    """JSON wrapped in markdown fences is extracted after stripping (common LLM output)."""
    wrapped = '```json\n{"candidates": [{"smiles": "CC"}]}\n```'
    cleaned = wrapped.replace("```json", "").replace("```", "").strip()
    result = extract_json(cleaned)
    assert result is not None
    assert result["candidates"][0]["smiles"] == "CC"


# ========================================
# Phase 1-B: calculate_metrics
# ========================================

def test_calculate_metrics_valid_smiles():
    """Benzene SMILES returns valid physical descriptors."""
    result = calculate_metrics("C1=CC=CC=C1")
    assert result is not None
    assert result["mw"] > 0
    assert result["rings"] == 1

def test_calculate_metrics_complex_smiles():
    """Naphthalene (2-ring PAH) ring count is correctly detected."""
    result = calculate_metrics("c1ccc2ccccc2c1")
    assert result is not None
    assert result["rings"] == 2

def test_calculate_metrics_invalid_smiles():
    """Invalid SMILES returns None."""
    result = calculate_metrics("INVALID_SMILES")
    assert result is None

def test_calculate_metrics_returns_expected_keys():
    """Return value contains all required descriptor keys."""
    result = calculate_metrics("C1=CC=CC=C1")
    assert result is not None
    for key in ["mw", "logp", "tpsa", "rings", "heavy_atoms"]:
        assert key in result


# ========================================
# Phase 1-C: API Endpoint Validation
# ========================================

def test_stream_missing_api_key_returns_401():
    """Request without API key returns 401."""
    response = client.post(
        "/api/mission/stream/m1",
        json={"philosophy": "STANDARD", "entropy_level": 0.1, "language": "en"}
    )
    assert response.status_code == 401

def test_stream_failed_status_returns_429():
    """Payload with status=FAILED returns 429."""
    response = client.post(
        "/api/mission/stream/m1?key=dummy_key_12345",
        json={"status": "FAILED", "language": "en"}
    )
    assert response.status_code == 429


# ========================================
# Phase 3: Router Logic
# ========================================

def test_design_router_serendipity_goes_to_mutate():
    """SERENDIPITY philosophy routes to mutate node."""
    state = {"philosophy": "SERENDIPITY", "status": "RUNNING"}
    assert design_router(state) == "mutate"

def test_design_router_standard_goes_to_physics():
    """STANDARD philosophy routes to physics node."""
    state = {"philosophy": "STANDARD", "status": "RUNNING"}
    assert design_router(state) == "physics"

def test_design_router_taboo_goes_to_physics():
    """TABOO philosophy routes to physics node."""
    state = {"philosophy": "TABOO", "status": "RUNNING"}
    assert design_router(state) == "physics"

def test_design_router_failed_goes_to_end():
    """status=FAILED routes to END to prevent zombie loops."""
    state = {"philosophy": "STANDARD", "status": "FAILED"}
    assert design_router(state) == "end"

def test_physics_router_valid_result_goes_to_sherry():
    """At least one valid candidate routes to sherry."""
    state = {
        "status": "RUNNING",
        "philosophy": "STANDARD",
        "physics_fail_count": 1,
        "metrics_results": [{"valid": True, "smiles": "C"}],
    }
    assert physics_router(state) == "sherry"

def test_physics_router_all_invalid_goes_to_vermouth():
    """All invalid candidates with fail_count < 3 routes back to vermouth."""
    state = {
        "status": "RUNNING",
        "philosophy": "STANDARD",
        "physics_fail_count": 1,
        "metrics_results": [{"valid": False}],
    }
    assert physics_router(state) == "vermouth"

def test_physics_router_third_failure_goes_to_sherry():
    """Third consecutive failure forces routing to sherry for deep analysis."""
    state = {
        "status": "RUNNING",
        "philosophy": "STANDARD",
        "physics_fail_count": 3,
        "metrics_results": [{"valid": False}],
    }
    assert physics_router(state) == "sherry"

def test_physics_router_serendipity_always_goes_to_sherry():
    """SERENDIPITY philosophy always routes to sherry regardless of validity."""
    state = {
        "status": "RUNNING",
        "philosophy": "SERENDIPITY",
        "physics_fail_count": 1,
        "metrics_results": [{"valid": False}],
    }
    assert physics_router(state) == "sherry"

def test_physics_router_failed_status_goes_to_end():
    """status=FAILED immediately routes to END."""
    state = {
        "status": "FAILED",
        "philosophy": "STANDARD",
        "physics_fail_count": 0,
        "metrics_results": [],
    }
    assert physics_router(state) == "end"