"""Gemini batch API evaluation."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from google import genai as genai_client
from google.genai import types

from clinvar.llm.evaluate_parallel import print_evaluation_summary
from clinvar.llm.parse import parse_llm_json_response, validate_evaluation_response

# BATCH FILE PREPARATION
# =============================================================================

def prepare_batch_jsonl(
    df_with_prompts: pd.DataFrame,
    system_prompt: str,
    output_path: str,
    variant_type: str = 'snp',
    model: str = 'gemini-3-flash-preview',
) -> str:
    """
    Prepare JSONL file for batch API submission.

    Each line follows the format required by Gemini Batch API:
    {"request": {"contents": [...], "system_instruction": {...}}, "custom_id": "..."}

    Args:
        df_with_prompts: DataFrame with 'llm_prompt' and '#VariationID' columns
        system_prompt: System instruction text
        output_path: Path for output JSONL file
        variant_type: 'snp' or 'indel'

    Returns:
        Path to created JSONL file
    """
    batch_requests = []

    for idx, row in df_with_prompts.iterrows():
        variation_id = str(row.get('#VariationID', row.name))
        user_prompt = row.get('llm_prompt', '')

        if not user_prompt:
            continue

        # Build request in Gemini batch format
        request = {
            "custom_id": variation_id,
            "request": {
                "model": f'models/{model}',
                "contents": [
                    {
                        "role": "user",
                        "parts": [{"text": user_prompt}]
                    }
                ],
                "system_instruction": {
                    "parts": [{"text": system_prompt}]
                },
                "generation_config": {
                    "temperature": 0.1,  # Low temperature for consistent JSON
                    "max_output_tokens": 4096,
                    "response_mime_type": "application/json"  # Request JSON output
                }
            }
        }

        batch_requests.append(request)

    # Write JSONL
    with open(output_path, 'w') as f:
        for req in batch_requests:
            f.write(json.dumps(req) + '\n')

    print(f"Created batch file: {output_path}")
    print(f"Total requests: {len(batch_requests)}")

    return output_path


def prepare_batch_jsonl_simple(
    df_with_prompts: pd.DataFrame,
    system_prompt: str,
    output_path: str
) -> str:
    """
    Simpler JSONL format - just the essentials.
    Adjust based on actual Gemini Batch API requirements.
    """
    with open(output_path, 'w') as f:
        for idx, row in df_with_prompts.iterrows():
            variation_id = str(row.get('#VariationID', row.name))
            user_prompt = row.get('llm_prompt', '')

            if not user_prompt:
                continue

            entry = {
                "custom_id": variation_id,
                "body": {
                    "contents": [
                        {"role": "user", "parts": [{"text": user_prompt}]}
                    ],
                    "systemInstruction": {
                        "parts": [{"text": system_prompt}]
                    },
                    "generationConfig": {
                        "temperature": 0.1,
                        "maxOutputTokens": 1024
                    }
                }
            }

            f.write(json.dumps(entry) + '\n')

    print(f"Created batch file: {output_path}")
    return output_path


# =============================================================================
# BATCH JOB MANAGEMENT
# =============================================================================

def upload_and_submit_batch(
    jsonl_path: str,
    display_name: str,
    api_key: str = None,
    model: str = "gemini-3-flash-preview"
) -> Dict[str, Any]:
    """
    Upload JSONL file and submit batch job.

    Args:
        jsonl_path: Path to prepared JSONL file
        display_name: Name for the batch job
        api_key: API key
        model: Model to use

    Returns:
        Dict with batch job info
    """
    client = genai_client.Client(api_key=api_key)

    # Upload the file
    print(f"Uploading {jsonl_path}...")
    uploaded_file = client.files.upload(
    file=jsonl_path,
    config=types.UploadFileConfig(mime_type="application/jsonl")
    )
    print(f"Uploaded: {uploaded_file.name}")

    # Create batch job
    print(f"Creating batch job...")
    batch_job = client.batches.create(
        model=model,
        src=uploaded_file.name,
        config={'display_name': display_name}
    )

    print(f"Batch job created: {batch_job.name}")
    print(f"Status: {batch_job.state}")

    return {
        'job_name': batch_job.name,
        'file_name': uploaded_file.name,
        'display_name': display_name,
        'state': str(batch_job.state),
        'created_time': time.time()
    }


def check_batch_status(
    job_name: str,
    api_key: str = None
) -> Dict[str, Any]:
    """
    Check status of a batch job.

    Args:
        job_name: Name of the batch job
        api_key: API key

    Returns:
        Dict with status info
    """
    client = genai_client.Client(api_key=api_key)

    batch_job = client.batches.get(name=job_name)

    status = {
        'name': batch_job.name,
        'state': str(batch_job.state),
        'display_name': getattr(batch_job, 'display_name', ''),
    }

    # Add progress info if available
    if hasattr(batch_job, 'request_counts'):
        counts = batch_job.request_counts
        status['total_requests'] = getattr(counts, 'total', 0)
        status['succeeded'] = getattr(counts, 'succeeded', 0)
        status['failed'] = getattr(counts, 'failed', 0)
        status['pending'] = getattr(counts, 'pending', 0)

    return status


def wait_for_batch_completion(
    job_name: str,
    api_key: str = None,
    poll_interval: int = 60,
    timeout: int = 7200  # 2 hours
) -> Dict[str, Any]:
    """
    Wait for batch job to complete, polling periodically.

    Args:
        job_name: Name of the batch job
        api_key: API key
        poll_interval: Seconds between status checks
        timeout: Maximum seconds to wait

    Returns:
        Final status dict
    """
    start_time = time.time()

    print(f"Waiting for batch job: {job_name}")
    print(f"Polling every {poll_interval}s, timeout {timeout}s")

    while True:
        status = check_batch_status(job_name, api_key)
        # state = status['state']
        state = str(status['state'])
        print(state)

        elapsed = time.time() - start_time

        # Progress display
        if 'total_requests' in status:
            succeeded = status.get('succeeded', 0)
            total = status.get('total_requests', 0)
            pct = (succeeded / total * 100) if total > 0 else 0
            print(f"\r[{elapsed/60:.1f}min] State: {state} | "
                  f"Progress: {succeeded}/{total} ({pct:.1f}%)", end='', flush=True)
        else:
            print(f"\r[{elapsed/60:.1f}min] State: {state}", end='', flush=True)

        # Check completion states
        if state in ['JobState.JOB_STATE_SUCCEEDED', 'JOB_STATE_SUCCEEDED', 'SUCCEEDED', 'STATE_SUCCEEDED']:
            print(f"\n✓ Batch job completed successfully!")
            return status

        if state in ['JOB_STATE_FAILED', 'FAILED', 'STATE_FAILED']:
            print(f"\n✗ Batch job failed!")
            return status

        if state in ['JOB_STATE_CANCELLED', 'CANCELLED', 'STATE_CANCELLED']:
            print(f"\n✗ Batch job was cancelled!")
            return status

        # Check timeout
        if elapsed > timeout:
            print(f"\n⚠ Timeout reached ({timeout}s)")
            return status

        time.sleep(poll_interval)


def download_batch_results(job_name: str, output_path: str, api_key: str = None) -> str:
    client = genai_client.Client(api_key=api_key)

    # 1. Retrieve the job
    batch_job = client.batches.get(name=job_name)

    # 2. Extract the filename from the .dest object
    # Based on your diagnostic, the attribute is 'file_name'
    result_file_name = None

    if hasattr(batch_job, 'dest') and batch_job.dest:
        result_file_name = batch_job.dest.file_name

    # Safety check if it's still None
    if not result_file_name:
        raise ValueError(f"Could not find 'file_name' in batch_job.dest. Object dump: {batch_job}")

    print(f"Downloading results from: {result_file_name}")

    # 3. Download using the 'file' keyword argument
    # This returns the raw bytes
    result_bytes = client.files.download(file=result_file_name)

    # 4. Save to disk
    with open(output_path, 'wb') as f:
        f.write(result_bytes)

    print(f"✓ Results saved to: {output_path}")
    return output_path


# =============================================================================
# RESULTS PARSING
# =============================================================================

def parse_batch_results(results_path: str, variant_type: str = "snp") -> pd.DataFrame:
    results = []

    with open(results_path, "r") as f:
        for line_num, line in enumerate(f):
            if not line.strip():
                continue

            try:
                entry = json.loads(line)
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse line {line_num}: {e}")
                continue

            custom_id = entry.get("custom_id", f"unknown_{line_num}")
            response = entry.get("response", {})

            if "error" in response:
                results.append(
                    {
                        "#VariationID": custom_id,
                        "concordance": "API_ERROR",
                        "explanation": str(response.get("error", "")),
                        "signal_category": "UNKNOWN",
                        "primary_mechanism": "",
                        "key_signals": "",
                        "rationale_mechanism": "",
                        "nt_missed": None,
                        "raw_response": json.dumps(response),
                        "parse_success": False,
                    }
                )
                continue

            try:
                if "candidates" in response:
                    text = response["candidates"][0]["content"]["parts"][0]["text"]
                elif "content" in response:
                    text = response["content"]["parts"][0]["text"]
                elif "text" in response:
                    text = response["text"]
                elif "body" in entry:
                    body = entry["body"]
                    if "candidates" in body:
                        text = body["candidates"][0]["content"]["parts"][0]["text"]
                    else:
                        text = str(body)
                else:
                    text = json.dumps(response)
            except (KeyError, IndexError, TypeError) as e:
                results.append(
                    {
                        "#VariationID": custom_id,
                        "concordance": "PARSE_ERROR",
                        "explanation": f"Could not extract text: {e}",
                        "signal_category": "UNKNOWN",
                        "primary_mechanism": "",
                        "key_signals": "",
                        "rationale_mechanism": "",
                        "nt_missed": None,
                        "raw_response": json.dumps(entry),
                        "parse_success": False,
                    }
                )
                continue

            parsed = parse_llm_json_response(text)
            if parsed is None:
                results.append(
                    {
                        "#VariationID": custom_id,
                        "concordance": "PARSE_ERROR",
                        "explanation": "Failed to parse JSON from response",
                        "signal_category": "UNKNOWN",
                        "primary_mechanism": "",
                        "key_signals": "",
                        "rationale_mechanism": "",
                        "nt_missed": None,
                        "raw_response": text,
                        "parse_success": False,
                    }
                )
                continue

            result = validate_evaluation_response(parsed, variant_type)
            result["#VariationID"] = custom_id
            result["raw_response"] = text
            result["parse_success"] = True
            results.append(result)

    df = pd.DataFrame(results)
    print(f"Parsed {len(df)} results from {results_path}")
    if "parse_success" in df.columns:
        print(f"Parse success rate: {df['parse_success'].mean() * 100:.1f}%")
    if "concordance" in df.columns:
        print("\nConcordance distribution:")
        print(df["concordance"].value_counts())
    return df


# =============================================================================
# COMPLETE BATCH PIPELINE
# =============================================================================

def run_batch_evaluation_pipeline(
    df_with_prompts: pd.DataFrame,
    system_prompt: str,
    variant_type: str = 'snp',
    batch_name: str = 'variant_evaluation',
    work_dir: str = './batch_eval',
    api_key: str = None,
    model: str = "gemini-3-flash-preview",
    wait_for_completion: bool = True,
    poll_interval: int = 60
) -> Optional[pd.DataFrame]:
    """
    Run complete batch evaluation pipeline.

    Args:
        df_with_prompts: DataFrame with 'llm_prompt' column
        system_prompt: System instruction
        variant_type: 'snp' or 'indel'
        batch_name: Name for the batch job
        work_dir: Directory for intermediate files
        api_key: API key
        model: Model to use
        wait_for_completion: If True, wait and return results
        poll_interval: Seconds between status checks

    Returns:
        DataFrame with results if wait_for_completion, else None
    """
    # Create work directory
    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)

    timestamp = time.strftime('%Y%m%d_%H%M%S')

    # Paths
    jsonl_path = work_path / f'{batch_name}_{timestamp}_requests.jsonl'
    results_path = work_path / f'{batch_name}_{timestamp}_results.jsonl'
    job_info_path = work_path / f'{batch_name}_{timestamp}_job_info.json'

    # Step 1: Prepare JSONL
    print("=" * 60)
    print("STEP 1: Preparing batch request file")
    print("=" * 60)
    prepare_batch_jsonl(
        df_with_prompts=df_with_prompts,
        system_prompt=system_prompt,
        output_path=str(jsonl_path),
        variant_type=variant_type,
        model=model,
    )

    # Step 2: Upload and submit
    print("\n" + "=" * 60)
    print("STEP 2: Uploading and submitting batch job")
    print("=" * 60)
    job_info = upload_and_submit_batch(
        jsonl_path=str(jsonl_path),
        display_name=f'{batch_name}_{timestamp}',
        api_key=api_key,
        model=model
    )

    # Save job info for later reference
    with open(job_info_path, 'w') as f:
        json.dump(job_info, f, indent=2)
    print(f"Job info saved to: {job_info_path}")

    if not wait_for_completion:
        print("\nBatch job submitted. Use check_batch_status() to monitor progress.")
        print(f"Job name: {job_info['job_name']}")
        return None

    # Step 3: Wait for completion
    print("\n" + "=" * 60)
    print("STEP 3: Waiting for batch completion")
    print("=" * 60)
    final_status = wait_for_batch_completion(
        job_name=job_info['job_name'],
        api_key=api_key,
        poll_interval=poll_interval
    )

    if final_status['state'] not in ['JobState.JOB_STATE_SUCCEEDED','JOB_STATE_SUCCEEDED', 'SUCCEEDED', 'STATE_SUCCEEDED']:
        print(f"\nBatch job did not complete successfully: {final_status['state']}")
        return None

    # Step 4: Download results
    print("\n" + "=" * 60)
    print("STEP 4: Downloading results")
    print("=" * 60)
    download_batch_results(
        job_name=job_info['job_name'],
        output_path=str(results_path),
        api_key=api_key
    )

    # Step 5: Parse results
    print("\n" + "=" * 60)
    print("STEP 5: Parsing results")
    print("=" * 60)
    results_df = parse_batch_results(
        results_path=str(results_path),
        variant_type=variant_type
    )

    # Save parsed results
    parquet_path = work_path / f'{batch_name}_{timestamp}_parsed.parquet'
    results_df.to_parquet(parquet_path)
    print(f"\nParsed results saved to: {parquet_path}")

    # Print summary
    print_evaluation_summary(results_df)

    return results_df


# =============================================================================
# UTILITY: Resume from job name
# =============================================================================

def resume_batch_evaluation(
    job_name: str,
    output_path: str,
    variant_type: str = 'snp',
    api_key: str = None,
    wait_if_pending: bool = True
) -> pd.DataFrame:
    """
    Resume/retrieve results from an existing batch job.

    Useful if you submitted a job and closed your session.

    Args:
        job_name: Name of the batch job
        output_path: Path to save results
        variant_type: 'snp' or 'indel'
        api_key: API key
        wait_if_pending: If True, wait for completion if job is still running

    Returns:
        DataFrame with results
    """
    # Check status
    status = check_batch_status(job_name, api_key)
    print(f"Job: {job_name}")
    print(f"State: {status['state']}\n")

    if status['state'] not in ['JOB_STATE_SUCCEEDED', 'SUCCEEDED', 'STATE_SUCCEEDED']:
        if wait_if_pending:
            print("Job not complete. Waiting...")
            status = wait_for_batch_completion(job_name, api_key)
        else:
            print("Job not complete yet.")
            return None

    # Download and parse
    results_jsonl = output_path.replace('.parquet', '_raw.jsonl')
    download_batch_results(job_name, results_jsonl, api_key)

    results_df = parse_batch_results(results_jsonl, variant_type)
    results_df.to_parquet(output_path)

    return results_df

