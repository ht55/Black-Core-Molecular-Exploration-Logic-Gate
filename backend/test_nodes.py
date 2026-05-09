# test_nodes.py
# Phase 2: LLM Mock Tests for design_node, mutate_node, analysis_node

import json
import pytest
from unittest.mock import AsyncMock, MagicMock
from graph_engine import (
    design_node,
    mutate_node,
    analysis_node,
    MASTER_ASSETS,
)


# ========================================
# Shared Fixtures
# ========================================

def make_state(overrides: dict = {}) -> dict:
    """Base state factory for all node tests."""
    base = {
        "mission_id": "m1",
        "language": "en",
        "philosophy": "STANDARD",
        "entropy_level": 0.3,
        "iteration": 0,
        "status": "RUNNING",
        "physics_fail_count": 0,
        "dialogue": [],
        "candidates": [],
        "metrics_results": [],
        "transformation_metrics": {},
        "scientific_report": "",
        "next_instruction": "",
        "distilled_lesson": "",
        "current_smiles": "",
        "base_scaffold_data": MASTER_ASSETS["m1"]["scaffolds"][0],
    }
    base.update(overrides)
    return base

def make_config(mock_llm) -> dict:
    """Config factory with a mock LLM."""
    return {"configurable": {"llm": mock_llm, "language": "en"}}

def make_llm(content: str) -> AsyncMock:
    """Create a mock LLM that returns the given content string."""
    mock_llm = AsyncMock()
    mock_llm.ainvoke.return_value.content = content
    return mock_llm


# ========================================
# design_node: Happy Path
# ========================================

@pytest.mark.asyncio
async def test_design_node_returns_two_candidates():
    """design_node returns two candidates on valid LLM response."""
    llm = make_llm(json.dumps({
        "candidates": [
            {"smiles": "C1=CC=CC=C1", "description": "Benzene derivative"},
            {"smiles": "c1ccc2ccccc2c1", "description": "Naphthalene derivative"}
        ],
        "text": "How about these?"
    }))
    result = await design_node(make_state(), make_config(llm))
    assert len(result["candidates"]) == 2
    assert result["candidates"][0]["smiles"] == "C1=CC=CC=C1"

@pytest.mark.asyncio
async def test_design_node_increments_iteration():
    """design_node increments iteration counter on success."""
    llm = make_llm(json.dumps({
        "candidates": [
            {"smiles": "C", "description": "Methane"},
            {"smiles": "CC", "description": "Ethane"}
        ],
        "text": "Done."
    }))
    result = await design_node(make_state({"iteration": 0}), make_config(llm))
    assert result["iteration"] == 1

@pytest.mark.asyncio
async def test_design_node_resets_next_instruction_on_success():
    """design_node clears next_instruction on successful design."""
    llm = make_llm(json.dumps({
        "candidates": [
            {"smiles": "C", "description": "Test"},
            {"smiles": "CC", "description": "Test2"}
        ],
        "text": "Done."
    }))
    result = await design_node(make_state(), make_config(llm))
    assert result["next_instruction"] == ""

@pytest.mark.asyncio
async def test_design_node_adds_strategy_field():
    """design_node adds a strategy field to each candidate for downstream compatibility."""
    llm = make_llm(json.dumps({
        "candidates": [
            {"smiles": "C", "description": "Test"},
            {"smiles": "CC", "description": "Test2"}
        ],
        "text": "Done."
    }))
    result = await design_node(make_state(), make_config(llm))
    for c in result["candidates"]:
        assert "strategy" in c


# ========================================
# design_node: Failure Cases
# ========================================

@pytest.mark.asyncio
async def test_design_node_json_parse_failure_increments_fail_count():
    """design_node increments physics_fail_count when LLM returns unparseable JSON."""
    llm = make_llm("This is not JSON at all.")
    result = await design_node(make_state(), make_config(llm))
    assert result["physics_fail_count"] == 1
    assert result["candidates"] == []

@pytest.mark.asyncio
async def test_design_node_429_sets_status_failed():
    """design_node sets status=FAILED on 429 / RESOURCE_EXHAUSTED error."""
    llm = AsyncMock()
    llm.ainvoke.side_effect = Exception("429 RESOURCE_EXHAUSTED: quota exceeded")
    result = await design_node(make_state(), make_config(llm))
    assert result["status"] == "FAILED"

@pytest.mark.asyncio
async def test_design_node_three_consecutive_failures_sets_status_failed():
    """design_node sets status=FAILED after 3 consecutive failures."""
    llm = make_llm("not json")
    state = make_state({"physics_fail_count": 2})
    result = await design_node(state, make_config(llm))
    assert result["status"] == "FAILED"

@pytest.mark.asyncio
async def test_design_node_already_failed_raises():
    """design_node raises RuntimeError immediately if status is already FAILED."""
    llm = make_llm("")
    state = make_state({"status": "FAILED"})
    with pytest.raises(RuntimeError):
        await design_node(state, make_config(llm))


# ========================================
# mutate_node: Happy Path
# ========================================

@pytest.mark.asyncio
async def test_mutate_node_returns_candidates():
    """mutate_node returns mutated candidates on valid LLM response."""
    llm = make_llm(json.dumps({
        "mutated_candidates": [
            {"smiles": "C1=CC=CC=C1F", "description": "Fluorinated accident"},
            {"smiles": "c1ccc2ccccc2c1Cl", "description": "Chlorinated accident"}
        ]
    }))
    state = make_state({
        "philosophy": "SERENDIPITY",
        "entropy_level": 0.8,
        "candidates": [
            {"smiles": "C1=CC=CC=C1", "description": "Original", "strategy": ""},
            {"smiles": "c1ccc2ccccc2c1", "description": "Original2", "strategy": ""}
        ]
    })
    result = await mutate_node(state, make_config(llm))
    assert len(result["candidates"]) == 2

@pytest.mark.asyncio
async def test_mutate_node_marks_serendipity_in_strategy():
    """mutate_node tags mutated candidates with SERENDIPITY_ACCIDENT in strategy."""
    llm = make_llm(json.dumps({
        "mutated_candidates": [
            {"smiles": "C1=CC=CC=C1F", "description": "Accident 1"},
            {"smiles": "C1CCCCC1", "description": "Accident 2"}
        ]
    }))
    state = make_state({
        "philosophy": "SERENDIPITY",
        "entropy_level": 1.0,
        "candidates": [
            {"smiles": "C", "description": "Base", "strategy": ""},
            {"smiles": "CC", "description": "Base2", "strategy": ""}
        ]
    })
    result = await mutate_node(state, make_config(llm))
    assert any("SERENDIPITY_ACCIDENT" in c.get("strategy", "") for c in result["candidates"])


# ========================================
# mutate_node: Failure Cases
# ========================================

@pytest.mark.asyncio
async def test_mutate_node_llm_error_preserves_original_candidates():
    """mutate_node preserves original candidates when LLM call fails."""
    llm = AsyncMock()
    llm.ainvoke.side_effect = Exception("LLM connection error")
    original_candidates = [
        {"smiles": "C1=CC=CC=C1", "description": "Original", "strategy": ""},
        {"smiles": "CC", "description": "Original2", "strategy": ""}
    ]
    state = make_state({
        "philosophy": "SERENDIPITY",
        "entropy_level": 0.8,
        "candidates": original_candidates
    })
    result = await mutate_node(state, make_config(llm))
    assert result["candidates"][0]["smiles"] == "C1=CC=CC=C1"
    assert result["candidates"][1]["smiles"] == "CC"


# ========================================
# analysis_node: Happy Path (with survivors)
# ========================================

@pytest.mark.asyncio
async def test_analysis_node_returns_scientific_report():
    """analysis_node returns a scientific report when valid candidates exist."""
    llm = make_llm(json.dumps({
        "analysis": "The mutation rate is 45.0%. Fitness score: 0.95.",
        "text": "The mutation rate is 45.0%. Structural analysis looks promising."
    }))
    state = make_state({
        "metrics_results": [{
            "valid": True,
            "smiles": "C1=CC=CC=C1",
            "strategy": "test",
            "descriptors": {"mw": 78.11, "psa": 0.0, "rings": 1, "rot_bonds": 0, "h_bond": 0}
        }]
    })
    result = await analysis_node(state, make_config(llm))
    assert result["scientific_report"] != ""

@pytest.mark.asyncio
async def test_analysis_node_sets_current_smiles():
    """analysis_node sets current_smiles to the best candidate's SMILES."""
    llm = make_llm(json.dumps({
        "analysis": "Mutation rate 45.0%. Fitness: 0.95.",
        "text": "The mutation rate is 45.0%. Looks good."
    }))
    state = make_state({
        "metrics_results": [{
            "valid": True,
            "smiles": "C1=CC=CC=C1",
            "strategy": "test",
            "descriptors": {"mw": 78.11, "psa": 0.0, "rings": 1, "rot_bonds": 0, "h_bond": 0}
        }]
    })
    result = await analysis_node(state, make_config(llm))
    assert result["current_smiles"] == "C1=CC=CC=C1"

@pytest.mark.asyncio
async def test_analysis_node_resets_fail_count_on_success():
    """analysis_node resets physics_fail_count to 0 when survivors exist."""
    llm = make_llm(json.dumps({
        "analysis": "Mutation rate 45.0%.",
        "text": "The mutation rate is 45.0%. Proceed."
    }))
    state = make_state({
        "physics_fail_count": 2,
        "metrics_results": [{
            "valid": True,
            "smiles": "C1=CC=CC=C1",
            "strategy": "test",
            "descriptors": {"mw": 78.11, "psa": 0.0, "rings": 1, "rot_bonds": 0, "h_bond": 0}
        }]
    })
    result = await analysis_node(state, make_config(llm))
    assert result["physics_fail_count"] == 0


# ========================================
# analysis_node: No Survivors (Physical Collapse)
# ========================================

@pytest.mark.asyncio
async def test_analysis_node_no_survivors_returns_report():
    """analysis_node returns a collapse report when no valid candidates exist."""
    llm = make_llm(json.dumps({
        "analysis": "All structures violated valence rules.",
        "text": "Vermouth, these structures are physically impossible."
    }))
    state = make_state({
        "metrics_results": [{"valid": False, "smiles": "INVALID"}]
    })
    result = await analysis_node(state, make_config(llm))
    assert result["scientific_report"] != ""

@pytest.mark.asyncio
async def test_analysis_node_no_survivors_status_remains_running():
    """analysis_node keeps status=RUNNING when no survivors (design loop continues)."""
    llm = make_llm(json.dumps({
        "analysis": "Collapse detected.",
        "text": "Fix the structure."
    }))
    state = make_state({
        "metrics_results": [{"valid": False, "smiles": "INVALID"}]
    })
    result = await analysis_node(state, make_config(llm))
    assert result["status"] == "RUNNING"

@pytest.mark.asyncio
async def test_analysis_node_llm_error_falls_back_gracefully():
    """analysis_node returns fallback data without crashing on LLM error."""
    llm = AsyncMock()
    llm.ainvoke.side_effect = Exception("LLM timeout")
    state = make_state({
        "metrics_results": [{
            "valid": True,
            "smiles": "C1=CC=CC=C1",
            "strategy": "test",
            "descriptors": {"mw": 78.11, "psa": 0.0, "rings": 1, "rot_bonds": 0, "h_bond": 0}
        }]
    })
    result = await analysis_node(state, make_config(llm))
    assert "scientific_report" in result
    assert result["scientific_report"] != ""