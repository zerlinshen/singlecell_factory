"""Tests for workflow.modular._densify_policy decision logic."""
from __future__ import annotations

import os
from unittest.mock import patch

import numpy as np
import pytest

from workflow.modular._densify_policy import (
    DensifyDecision,
    plan_densify,
    planned_bytes,
)

_FLOAT32_ITEMSIZE = np.dtype(np.float32).itemsize  # 4 bytes


def test_planned_bytes_small():
    assert planned_bytes((100, 200), np.float32) == 100 * 200 * _FLOAT32_ITEMSIZE


def test_planned_bytes_float64():
    assert planned_bytes((10, 10), np.float64) == 10 * 10 * 8


def test_plan_densify_go_for_small_array():
    """Small array well within caps → GO."""
    # 1000 x 100 float32 = 400 KB, tiny
    decision = plan_densify((1000, 100), np.float32, reason="test small array")
    assert decision == DensifyDecision.GO


def test_plan_densify_abort_above_hard_cap():
    """Array >= hard cap → ABORT regardless of free RAM."""
    hard_cap = 12 * 1024**3  # 12 GiB
    # Use env override to lower hard cap to 1 byte for predictable test
    with patch.dict(os.environ, {"SC_DENSIFY_HARD_CAP_BYTES": "1",
                                  "SC_DENSIFY_SOFT_CAP_BYTES": "1"}):
        decision = plan_densify((100, 100), np.float32, reason="test abort")
    assert decision == DensifyDecision.ABORT


def test_plan_densify_chunk_between_soft_and_hard():
    """Array between soft and hard caps → CHUNK."""
    # Set soft=1 byte, hard=very large so nbytes falls between them
    nbytes = planned_bytes((100, 100), np.float32)  # 40000 bytes
    with patch.dict(os.environ, {
        "SC_DENSIFY_SOFT_CAP_BYTES": "1",
        "SC_DENSIFY_HARD_CAP_BYTES": str(nbytes + 1),
    }):
        decision = plan_densify((100, 100), np.float32, reason="test chunk")
    assert decision == DensifyDecision.CHUNK


def test_plan_densify_chunk_when_exceeds_80pct_free_ram():
    """When planned > 80% of free RAM but < hard cap → CHUNK."""
    # Free RAM = 1000 bytes; planned = 900 bytes (90% of free) → CHUNK
    with patch.dict(os.environ, {
        "SC_DENSIFY_SOFT_CAP_BYTES": str(10 * 1024**3),  # 10 GiB soft (well above)
        "SC_DENSIFY_HARD_CAP_BYTES": str(12 * 1024**3),  # 12 GiB hard
    }):
        with patch("workflow.modular._densify_policy._free_bytes", return_value=1000):
            # 1 x 225 float32 = 900 bytes > 80% of 1000
            decision = plan_densify((1, 225), np.float32, reason="test ram guard")
    assert decision == DensifyDecision.CHUNK


def test_plan_densify_go_when_free_ram_is_sufficient():
    """When planned is well under 80% of free RAM and under soft cap → GO."""
    with patch.dict(os.environ, {
        "SC_DENSIFY_SOFT_CAP_BYTES": str(10 * 1024**3),
        "SC_DENSIFY_HARD_CAP_BYTES": str(12 * 1024**3),
    }):
        with patch("workflow.modular._densify_policy._free_bytes", return_value=10 * 1024**3):
            # Very small array: 100 x 100 float32 = 40000 bytes, << 80% of 10 GiB
            decision = plan_densify((100, 100), np.float32, reason="test go")
    assert decision == DensifyDecision.GO


def test_densify_decision_enum_values():
    assert DensifyDecision.GO == "go"
    assert DensifyDecision.CHUNK == "chunk"
    assert DensifyDecision.ABORT == "abort"
