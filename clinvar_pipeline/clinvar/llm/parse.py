"""Shared LLM response parsing and validation."""

from __future__ import annotations

import json
import re
from typing import Any, Optional


def parse_llm_json_response(response_text: str) -> Optional[dict[str, Any]]:
    if not response_text:
        return None

    try:
        return json.loads(response_text)
    except json.JSONDecodeError:
        pass

    code_block_pattern = r"```(?:json)?\s*\n?(.*?)\n?```"
    for match in re.findall(code_block_pattern, response_text, re.DOTALL):
        try:
            return json.loads(match.strip())
        except json.JSONDecodeError:
            continue

    json_pattern = r"\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}"
    for match in re.findall(json_pattern, response_text, re.DOTALL):
        try:
            return json.loads(match)
        except json.JSONDecodeError:
            continue

    return None


def validate_evaluation_response(
    parsed: dict[str, Any],
    variant_type: str = "snp",
) -> dict[str, Any]:
    expected_fields: dict[str, Any] = {
        "concordance": "PARSE_ERROR",
        "explanation": "",
        "signal_category": "UNKNOWN",
        "primary_mechanism": "",
        "key_signals": "",
        "rationale_mechanism": "",
        "nt_missed": None,
        "notes": "",
    }
    if variant_type == "indel":
        expected_fields["embedding_impact"] = "UNKNOWN"

    valid_concordance = {"CONCORDANT", "PARTIAL", "DISCORDANT", "NOT_APPLICABLE"}
    valid_signal_category = {"STRONG", "MODERATE", "WEAK", "ABSENT", "UNKNOWN"}
    valid_embedding_impact = {"HIGH", "MODERATE", "LOW", "NONE", "UNKNOWN"}

    result = {field: parsed.get(field, default) for field, default in expected_fields.items()}

    if result["concordance"] not in valid_concordance:
        result["concordance"] = "PARSE_ERROR"
    if result["signal_category"] not in valid_signal_category:
        result["signal_category"] = "UNKNOWN"
    if variant_type == "indel" and result.get("embedding_impact") not in valid_embedding_impact:
        result["embedding_impact"] = "UNKNOWN"
    if result["nt_missed"] is not None:
        result["nt_missed"] = bool(result["nt_missed"])

    return result
