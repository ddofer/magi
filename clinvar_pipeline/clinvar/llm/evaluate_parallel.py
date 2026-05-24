"""Threaded Gemini evaluation with checkpointing."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Optional

import pandas as pd
import threading

from clinvar.llm.parse import parse_llm_json_response, validate_evaluation_response
from clinvar.prompts import get_system_prompt

_thread_local = threading.local()


def get_thread_client(api_key: str):
    if not hasattr(_thread_local, "client"):
        from google import genai

        _thread_local.client = genai.Client(api_key=api_key)
    return _thread_local.client


def call_llm_api_robust_threaded(
    user_content: str,
    system_content: str,
    model_name: str,
    api_key: str,
    max_retries: int = 10,
) -> Optional[str]:
    from google.api_core import exceptions

    client = get_thread_client(api_key)

    for attempt in range(max_retries):
        try:
            resp = client.models.generate_content(
                model=model_name,
                contents=user_content,
                config={"system_instruction": system_content},
            )
            return resp.text.strip() if resp and resp.text else ""
        except exceptions.ResourceExhausted:
            wait_time = (2**attempt) * 5
            print(f"\n[Rate Limit] sleep {wait_time}s (attempt {attempt + 1}/{max_retries})")
            time.sleep(wait_time)
        except exceptions.InvalidArgument as e:
            print(f"\n[Invalid Request] {e}")
            return None
        except Exception as e:
            print(f"\n[Error] {e}")
            time.sleep(5)

    return None


def _error_result(variation_id, concordance: str, explanation: str, notes: str, raw=None):
    return {
        "#VariationID": variation_id,
        "concordance": concordance,
        "explanation": explanation,
        "signal_category": "UNKNOWN",
        "primary_mechanism": "",
        "key_signals": "",
        "rationale_mechanism": "",
        "nt_missed": None,
        "notes": notes,
        "raw_response": raw,
        "parse_success": False,
    }


def evaluate_single_variant(row, system_prompt, variant_type, model_name, api_key):
    variation_id = row.get("#VariationID", row.name)
    user_prompt = row.get("llm_prompt", "")
    if not user_prompt:
        return _error_result(variation_id, "ERROR", "No prompt available", "Missing llm_prompt column")

    raw_response = call_llm_api_robust_threaded(
        user_content=user_prompt,
        system_content=system_prompt,
        model_name=model_name,
        api_key=api_key,
    )
    if raw_response is None:
        return _error_result(variation_id, "API_ERROR", "API call failed", "API returned None")

    parsed = parse_llm_json_response(raw_response)
    if parsed is None:
        return _error_result(
            variation_id,
            "PARSE_ERROR",
            "Failed to parse JSON",
            f"Raw response: {raw_response[:500]}",
            raw=raw_response,
        )

    validated = validate_evaluation_response(parsed, variant_type)
    validated["#VariationID"] = variation_id
    validated["raw_response"] = raw_response
    validated["parse_success"] = True
    return validated


def evaluate_variants_parallel(
    df: pd.DataFrame,
    system_prompt: str,
    variant_type: str,
    model_name: str,
    api_key: str,
    max_workers: int = 10,
    save_checkpoint_every: int = 50,
    checkpoint_path: str | None = None,
    resume_from_checkpoint: bool = True,
) -> pd.DataFrame:
    results: list[dict] = []
    evaluated_ids: set[str] = set()

    if checkpoint_path and resume_from_checkpoint:
        try:
            ck = pd.read_parquet(checkpoint_path)
            results = ck.to_dict("records")
            evaluated_ids = set(ck["#VariationID"].astype(str))
            print(f"Resumed {len(evaluated_ids)} already done")
        except FileNotFoundError:
            pass

    remaining = []
    for row_idx, row in df.iterrows():
        vid = str(row.get("#VariationID", row_idx))
        if vid not in evaluated_ids:
            remaining.append((row_idx, row))

    total_target = len(evaluated_ids) + len(remaining)
    done_counter = 0
    start = time.time()

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        future_to_vid = {}
        for row_idx, row in remaining:
            vid = str(row.get("#VariationID", row_idx))
            fut = ex.submit(
                evaluate_single_variant,
                row,
                system_prompt,
                variant_type,
                model_name,
                api_key,
            )
            future_to_vid[fut] = vid

        for fut in as_completed(future_to_vid):
            vid = future_to_vid[fut]
            try:
                res = fut.result()
            except Exception as e:
                res = _error_result(vid, "ERROR", f"Worker exception: {e}", "Worker crashed")

            results.append(res)
            done_counter += 1
            elapsed = time.time() - start
            rate_per_min = (done_counter / elapsed) * 60 if elapsed > 0 else 0
            status_char = "✓" if res.get("parse_success") else "✗"
            print(
                f"\r[{len(results)}/{total_target}] {status_char} {vid} | {rate_per_min:.1f}/min",
                end="",
                flush=True,
            )

            if checkpoint_path and (len(results) % save_checkpoint_every == 0):
                pd.DataFrame(results).to_parquet(checkpoint_path)
                print(f"\n[Checkpoint] saved {len(results)} to {checkpoint_path}")

    print(f"\nDone. Total results: {len(results)}")
    out = pd.DataFrame(results)
    if checkpoint_path:
        out.to_parquet(checkpoint_path)
    return out


def summarize_evaluation_results(results_df: pd.DataFrame) -> dict[str, Any]:
    total = len(results_df)
    summary: dict[str, Any] = {
        "total_variants": total,
        "parse_success_rate": (
            results_df["parse_success"].mean() if "parse_success" in results_df.columns else None
        ),
        "concordance_distribution": results_df["concordance"].value_counts().to_dict(),
        "signal_category_distribution": results_df["signal_category"].value_counts().to_dict(),
        "nt_missed_rate": (
            results_df["nt_missed"].mean() if "nt_missed" in results_df.columns else None
        ),
    }
    concordance_counts = results_df["concordance"].value_counts()
    for cat in ["CONCORDANT", "PARTIAL", "DISCORDANT", "NOT_APPLICABLE"]:
        count = concordance_counts.get(cat, 0)
        summary[f"pct_{cat.lower()}"] = count / total * 100 if total > 0 else 0
    return summary


def print_evaluation_summary(results_df: pd.DataFrame) -> None:
    summary = summarize_evaluation_results(results_df)
    print("=" * 60)
    print("EVALUATION SUMMARY")
    print("=" * 60)
    print(f"Total variants evaluated: {summary['total_variants']}")
    if summary["parse_success_rate"] is not None:
        print(f"Parse success rate: {summary['parse_success_rate'] * 100:.1f}%")
    print("\nConcordance Distribution:")
    print("-" * 30)
    for cat in ["CONCORDANT", "PARTIAL", "DISCORDANT", "NOT_APPLICABLE"]:
        pct = summary.get(f"pct_{cat.lower()}", 0)
        count = summary["concordance_distribution"].get(cat, 0)
        bar = "█" * int(pct / 2)
        print(f"  {cat:16} {count:4} ({pct:5.1f}%) {bar}")
    error_cats = ["API_ERROR", "PARSE_ERROR", "ERROR"]
    error_count = sum(summary["concordance_distribution"].get(c, 0) for c in error_cats)
    if error_count:
        print(f"\n  Errors: {error_count}")
    print("\nSignal Category Distribution:")
    print("-" * 30)
    for cat, count in summary["signal_category_distribution"].items():
        pct = count / summary["total_variants"] * 100
        print(f"  {cat:12} {count:4} ({pct:5.1f}%)")
    if summary["nt_missed_rate"] is not None:
        print(f"\nNT Missed Rate: {summary['nt_missed_rate'] * 100:.1f}%")
    print("=" * 60)


def run_evaluation_pipeline(
    df_with_prompts: pd.DataFrame,
    variant_type: str = "snp",
    output_path: str | None = None,
    checkpoint_path: str | None = None,
    model_name: str = "gemini-3-flash-preview",
    api_key: str | None = None,
    max_workers: int = 10,
    sample_size: int | None = None,
) -> pd.DataFrame:
    if not api_key:
        raise ValueError("api_key is required (set GEMINI_API_KEY env var)")

    system_prompt = get_system_prompt(variant_type)
    df_eval = (
        df_with_prompts.sample(n=sample_size, random_state=42)
        if sample_size and sample_size < len(df_with_prompts)
        else df_with_prompts
    )
    if sample_size and sample_size < len(df_with_prompts):
        print(f"Sampled {sample_size} variants for evaluation")

    print(f"Starting evaluation of {len(df_eval)} {variant_type} variants")
    print(f"Model: {model_name}")
    print()

    results_df = evaluate_variants_parallel(
        df=df_eval,
        system_prompt=system_prompt,
        variant_type=variant_type,
        model_name=model_name,
        api_key=api_key,
        max_workers=max_workers,
        checkpoint_path=checkpoint_path,
    )
    print_evaluation_summary(results_df)

    if output_path:
        results_df.to_parquet(output_path)
        print(f"\nResults saved to: {output_path}")
        csv_path = output_path.replace(".parquet", ".csv")
        results_df.to_csv(csv_path, index=False)
        print(f"CSV saved to: {csv_path}")

    return results_df


def merge_evaluation_results(
    original_df: pd.DataFrame,
    results_df: pd.DataFrame,
    on: str = "#VariationID",
) -> pd.DataFrame:
    merge_cols = [
        on,
        "concordance",
        "explanation",
        "signal_category",
        "primary_mechanism",
        "key_signals",
        "rationale_mechanism",
        "nt_missed",
        "parse_success",
    ]
    if "embedding_impact" in results_df.columns:
        merge_cols.append("embedding_impact")

    results_slim = results_df[[c for c in merge_cols if c in results_df.columns]].copy()
    original_df = original_df.copy()
    original_df[on] = original_df[on].astype(str)
    results_slim[on] = results_slim[on].astype(str)
    return original_df.merge(results_slim, on=on, how="left", suffixes=("", "_eval"))
