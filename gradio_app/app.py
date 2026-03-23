#!/usr/bin/env python3
"""
MAGI Variant Interpreter - Gradio Web Interface
===============================================
Interactive web app for scoring and summarizing the predicted effects of
human genomic variants.

The app uses the external NTv3 foundation model together with MAGI scoring,
annotation, and interpretation layers.

Usage:
    python app.py

For HuggingFace Spaces deployment with ZeroGPU:
    Uses @spaces.GPU decorator for GPU bursting
"""

import os
import warnings
from pathlib import Path

import torch
import gradio as gr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import our modules
from inference import (
    CONTEXT_LEN,
    predict_variants,
    load_model_and_resources,
    _MODEL_CACHE,
    SUPPORTED_UI_SPECIES,
)
from annotation import annotate_dataframe, load_mane_data, _MANE_CACHE
from analysis import (
    rank_top_disrupted_tracks,
    build_top_track_table,
    compute_impact_scores,
    extract_top_summary_signals,
    summarize_variant,
    get_top_bigwig_descriptions,
    make_fingerprint_plot,
    identify_delta_columns,
    format_summary_table,
)
from tracks import generate_region_tracks_plot, get_track_view_bounds
from interpretation import build_signal_interpretation

# Suppress warnings for cleaner UI
warnings.filterwarnings("ignore")

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
EXAMPLES_DIR = BASE_DIR / "examples"
METADATA_FILE = DATA_DIR / "functional_tracks_metadata_human.csv"
DEFAULT_ZOOM_BP = 128
MIN_ZOOM_BP = 8
INITIAL_ZOOM_MAX_BP = CONTEXT_LEN // 2
ANNOTATION_FLAGS = [
    "gene",
    "mRNA",
    "mRNA_promoter",
    "mRNA_exon",
    "coding_sequence",
    "start_codon",
    "stop_codon",
    "five_prime_UTR",
    "three_prime_UTR",
    "mRNA_intron",
    "mRNA_splice",
    "lncRNA",
    "lncRNA_promoter",
    "lncRNA_exon",
    "snRNA",
    "snRNA_promoter",
    "snRNA_exon",
    "antisenseRNA",
    "antisenseRNA_promoter",
    "antisenseRNA_exon",
    "telomeraseRNA",
    "telomeraseRNA_promoter",
    "telomeraseRNA_exon",
    "RNaseMRPRNA",
    "RNaseMRPRNA_promoter",
    "RNaseMRPRNA_exon",
    "snoRNA",
    "snoRNA_promoter",
    "snoRNA_exon",
    "other",
]


def _export_figure_png(fig, file_name: str, dpi: int = 300):
    """Persist a Matplotlib figure for figure-ready download."""
    if fig is None:
        return None

    output_path = BASE_DIR / file_name
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    return str(output_path)


def _is_human_species(species: str) -> bool:
    return str(species).strip() == "human"


def _species_label(species: str) -> str:
    return SUPPORTED_UI_SPECIES.get(
        str(species).strip(), str(species).replace("_", " ").title()
    )


def _apply_species_annotations(results_df: pd.DataFrame, species: str) -> pd.DataFrame:
    if _is_human_species(species):
        return annotate_dataframe(results_df)

    annotated = results_df.copy()
    for col in ANNOTATION_FLAGS:
        annotated[col] = 0
    annotated["transcript_set"] = [set() for _ in range(len(annotated))]
    annotated["promoter_transcript_set"] = [set() for _ in range(len(annotated))]
    annotated["region"] = "non_human_unavailable"
    annotated["region_class"] = "NON_HUMAN"
    annotated["gene_name"] = "N/A (non-human)"
    return annotated


def _build_zoom_slider_update(requested_zoom_bp):
    """Update the zoom slider to the current variant's available track-view span."""
    try:
        requested_zoom_bp = max(MIN_ZOOM_BP, int(requested_zoom_bp))
    except (TypeError, ValueError):
        requested_zoom_bp = DEFAULT_ZOOM_BP

    view_bounds = get_track_view_bounds()
    max_radius = view_bounds.get("max_radius")
    if max_radius is None:
        return gr.update()

    return gr.update(
        value=min(requested_zoom_bp, max_radius),
        maximum=max_radius,
        info=f"Current available radius for this prediction: up to {max_radius:,} bp",
    )


def predict_single_variant_ui(species, chrom, pos, ref, alt, zoom_bp=DEFAULT_ZOOM_BP):
    """UI wrapper: run prediction and prepare user-facing outputs."""
    result = predict_single_variant(species, chrom, pos, ref, alt, zoom_bp=zoom_bp)
    # Engine returns 9 values
    (
        summary_md,
        interpretation_md,
        top_table_df,
        fingerprint_fig,
        region_tracks_fig,
        csv_path,
        bed_table_df,
        mlm_md,
        ranked,
    ) = result

    tracks_png_path = None
    try:
        tracks_png_path = _export_figure_png(
            region_tracks_fig, "last_variant_tracks.png"
        )
    except Exception as exc:
        print(f"\u26a0\ufe0f  Track PNG export failed: {exc}")

    zoom_update = gr.update()
    if isinstance(ranked, list):
        zoom_update = _build_zoom_slider_update(zoom_bp)

    return (
        summary_md,
        interpretation_md,
        top_table_df,
        fingerprint_fig,
        region_tracks_fig,
        csv_path,
        tracks_png_path,
        bed_table_df,
        mlm_md,
        zoom_update,
        ranked[:50]
        if isinstance(ranked, list)
        else ranked,  # → gr.State for zoom slider
    )


def _update_region_zoom(zoom_bp, ranked, track_filter="All"):
    """Re-render region tracks plot when the zoom slider or filter changes."""
    if not ranked:
        return None, None
    try:
        filtered = ranked
        if track_filter and track_filter != "All":
            filtered = [r for r in ranked if r["track_type"] == track_filter]
        if not filtered:
            filtered = ranked
        fig = generate_region_tracks_plot(
            ranked_tracks=filtered,
            max_ranked_tracks=10,
            visible_radius_bp=max(int(zoom_bp), 8),
        )
        png_path = _export_figure_png(fig, "last_variant_tracks.png")
        return fig, png_path
    except Exception as exc:
        print(f"\u26a0\ufe0f  Zoom re-render failed: {exc}")
        return None, None


def predict_single_variant_wrapper(
    species, chrom, pos, ref, alt, zoom_bp=DEFAULT_ZOOM_BP
):
    """Wrapper with ZeroGPU decorator if available."""
    return predict_single_variant_ui(species, chrom, pos, ref, alt, zoom_bp)


# Check if running on HuggingFace Spaces with ZeroGPU
try:
    import spaces

    HAS_SPACES = True
    print("\u2705 Running with HuggingFace Spaces ZeroGPU support")
except ImportError:
    spaces = None
    HAS_SPACES = False
    print("\u26a0\ufe0f  Not running on HuggingFace Spaces - using local GPU/CPU")


# Apply decorators if available
if HAS_SPACES and spaces is not None:
    predict_single_variant_wrapper = spaces.GPU(predict_single_variant_wrapper)

# Device selection: prefer CUDA if available, fall back to CPU
if os.environ.get("SPACES_GPU") or (HAS_SPACES and torch.cuda.is_available()):
    DEVICE = "cuda"
elif torch.cuda.is_available():
    DEVICE = "cuda"
else:
    DEVICE = "cpu"
print(f"⚙️  Using device: {DEVICE}")

# Load metadata for BigWig track descriptions
TRACK_METADATA = pd.read_csv(METADATA_FILE)

# Pre-build metadata lookup dictionary for O(1) access
TRACK_META_DICT = {
    row["file_id"]: {
        "tissue": str(row.get("tissue", "") or ""),
        "assay": str(row.get("assay", "") or ""),
        "target": str(row.get("experiment_target", "") or ""),
    }
    for _, row in TRACK_METADATA.iterrows()
}


# ============================================================================
# STARTUP - LOAD MODEL AND ANNOTATION DATA
# ============================================================================
def initialize_app():
    """Load model and annotation data at startup."""
    print("🚀 Initializing MAGI Variant Interpreter...")
    try:
        load_model_and_resources(device=DEVICE)
        load_mane_data()
        print("✅ App ready!")
        return True
    except Exception as e:
        print(f"❌ Model initialization failed: {e}")
        import traceback

        traceback.print_exc()
        return False


# Initialize on import
_APP_READY = initialize_app()


# ============================================================================
# UTILITY HELPERS
# ============================================================================
def _fmt(val, fmt: str = ".4f", na: str = "N/A") -> str:
    """Safely format a numeric value; return 'na' string for NaN/None."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return na
    try:
        return format(float(val), fmt)
    except (TypeError, ValueError):
        return na


def _render_track_signal_block(title: str, items) -> str:
    """Render compact BED/BigWig bullets for the summary card."""
    if not items:
        return f"**{title}**\n- none above 3% absolute change"

    lines = [f"**{title}**"]
    for item in items:
        lines.append(
            "- "
            f"{item['label']}: {item['delta']:+.1%} "
            f"({_fmt(item['ref_val'])} → {_fmt(item['alt_val'])})"
        )
    return "\n".join(lines)


def _render_mlm_signal_block(items) -> str:
    """Render compact sequence-model bullets for the summary card."""
    if not items:
        return "**Top sequence-model signals**\n- none above 0.03 absolute magnitude"

    lines = ["**Top sequence-model signals**"]
    for item in items:
        value_fmt = (
            f"{item['value']:+.4f}"
            if item.get("signed", False)
            else f"{item['value']:.4f}"
        )
        lines.append(f"- {item['label']}: {value_fmt}")
    return "\n".join(lines)


def _build_top_signal_summary(signals, include_bigwig: bool = True) -> str:
    """Build the top-signal markdown block for the summary card."""
    blocks = [
        "#### Top Signals",
        _render_track_signal_block("Top BED signals", signals.get("bed", [])),
        _render_mlm_signal_block(signals.get("mlm", [])),
    ]
    if include_bigwig:
        blocks.insert(
            2,
            _render_track_signal_block("Top BigWig signals", signals.get("bigwig", [])),
        )
    return "\n\n".join(blocks)


def _validate_allele(allele: str) -> bool:
    """Return True if allele contains only valid nucleotide characters."""
    return bool(allele) and all(c in "ACGTNacgtn" for c in allele.strip())


# ============================================================================
# SINGLE VARIANT PREDICTION
# ============================================================================
def predict_single_variant(
    species: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    zoom_bp: int = DEFAULT_ZOOM_BP,
):
    """
    Predict functional impact for a single variant.

    Returns:
        Tuple of (summary_md, interpretation_md, top_table_df,
                ranked_tracks=ranked,
                max_ranked_tracks=15,
                  bed_table_df, mlm_md, ranked_tracks)
    """
    _none9 = (None,) * 9

    try:
        # Validate inputs
        if not chrom or not ref or not alt:
            return (
                "❌ Error: Please provide all required fields (Chromosome, Position, Ref, Alt)",
                *_none9[1:],
            )

        if pos <= 0:
            return (
                "❌ Error: Position must be a positive integer",
                *_none9[1:],
            )

        ref_clean = ref.upper().strip()
        alt_clean = alt.upper().strip()

        if not _validate_allele(ref_clean):
            return (
                f"❌ Error: REF allele '{ref}' contains invalid characters. Use A/C/G/T/N only.",
                *_none9[1:],
            )
        if not _validate_allele(alt_clean):
            return (
                f"❌ Error: ALT allele '{alt}' contains invalid characters. Use A/C/G/T/N only.",
                *_none9[1:],
            )

        # Create input DataFrame
        input_df = pd.DataFrame(
            [
                {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref_clean,
                    "alt": alt_clean,
                }
            ]
        )

        # Run inference
        species = str(species).strip()
        print(f"🔬 Predicting {species} {chrom}:{pos} {ref_clean}>{alt_clean}...")
        results_df = predict_variants(input_df, device=DEVICE, species=species)

        # Annotate
        results_df = _apply_species_annotations(results_df, species)

        # Compute impact scores
        results_df = compute_impact_scores(results_df)

        # Extract data
        row = results_df.iloc[0]

        gene = row.get("gene_name", "Unknown")
        region = row.get("region_class", "OTHER")
        impact_bed = row.get("Impact_Score_BED", np.nan)
        impact_bw = row.get("Impact_Score_BW", np.nan)
        magi_score = row.get("Global_z_sum_log", np.nan)
        llr = row.get("LLR", np.nan)
        kl_mean = row.get("MLM_KL_mean", np.nan)
        indel_size = row.get("indel_size", 0)
        transcript_set = row.get("transcript_set", set())

        variant_type = "Indel" if indel_size != 0 else "SNP"
        indel_label = f" ({int(indel_size):+d} bp)" if indel_size != 0 else ""

        if pd.isna(impact_bed):
            impact_level = "**unknown**"
        elif impact_bed > 0.1:
            impact_level = "**high**"
        elif impact_bed > 0.05:
            impact_level = "**moderate**"
        else:
            impact_level = "**low**"

        llr_note = " *(N/A for indels)*" if pd.isna(llr) else ""

        transcript_str = ""
        if transcript_set and len(transcript_set) > 0:
            transcript_ids = ", ".join(sorted(list(transcript_set))[:3])
            transcript_str = f"  \n**Transcripts:** `{transcript_ids}`"

        context_start = max(1, pos - CONTEXT_LEN // 2)
        context_end = pos + CONTEXT_LEN // 2

        # Inline genomic context tags (replaces separate Annotation Detail accordion)
        active_flags = [
            flag
            for flag in ANNOTATION_FLAGS
            if flag in row.index and row.get(flag, 0) == 1
        ]
        if not _is_human_species(species):
            context_tags = "Transcript-level human MANE annotation is not available for this species in the web app."
        elif active_flags:
            context_tags = " · ".join(f"`{f}`" for f in active_flags)
        else:
            context_tags = "Intergenic (no annotated features)"

        # === 1. Unified ranking → fingerprint, table, region tracks ===
        bed_names = _MODEL_CACHE.get("bed_names", [])
        bw_names_selected = _MODEL_CACHE.get("selected_bw_indices", [])
        bw_names_all = _MODEL_CACHE.get("bigwig_names", [])
        bw_names_filtered = (
            [bw_names_all[i] for i in bw_names_selected]
            if _is_human_species(species) and bw_names_selected and bw_names_all
            else None
        )

        ranked = rank_top_disrupted_tracks(
            row,
            bed_names,
            bw_names_filtered,
            metadata_df=TRACK_METADATA,
            metadata_dict=TRACK_META_DICT,
            top_k=None,
        )

        top_table_df = build_top_track_table(
            ranked,
            max_rows=50,
            min_rows_by_type={"BED": 5},
        )

        fingerprint_fig = make_fingerprint_plot(ranked, top_k=15)

        interpretation_md = build_signal_interpretation(row, ranked, variant_type)

        summary_signals = extract_top_summary_signals(
            row,
            ranked,
            min_abs_threshold=0.03,
        )
        top_signal_md = _build_top_signal_summary(
            summary_signals,
            include_bigwig=_is_human_species(species),
        )

        species_name = _species_label(species)
        bigwig_value = (
            _fmt(impact_bw)
            if _is_human_species(species)
            else "N/A (human-only context tracks)"
        )
        annotation_note = (
            ""
            if _is_human_species(species)
            else "\n**Non-human note:** BigWig context tracks and MANE transcript annotation are not available in this app for non-human species."
        )

        summary_md = f"""
### Variant Summary

**Variant:** `{chrom}:{pos} {ref_clean}>{alt_clean}` ({variant_type}{indel_label})  
**Species:** {species_name} (`{species}`)  
**Gene:** {gene if gene else "Intergenic"}  
**Region:** {region}{transcript_str}  
**Genomic context:** {context_tags}

{top_signal_md}

| Metric | Value |
|--------|-------|
| MAGI score (`Global_z_sum_log`) | {_fmt(magi_score)} |
| BED Impact Score | {_fmt(impact_bed)} |
| BigWig Impact Score | {bigwig_value} |
| LLR | {_fmt(llr)}{llr_note} |
| KL Divergence (mean) | {_fmt(kl_mean)} |

**Analysis window:** `{chrom}:{context_start:,}-{context_end:,}` ({context_end - context_start:,} bp)  
**BED-based summary tier:** {impact_level}. This tier reflects `Impact_Score_BED` only.{annotation_note}
"""

        # === 2. BED Table (all 21 elements) ===
        bed_table = []
        for name in bed_names:
            ref_col = f"REF_BED_{name}"
            delta_col = f"D_BED_{name}"
            if ref_col in row.index and delta_col in row.index:
                ref_val = row[ref_col]
                delta_val = row[delta_col]
                alt_val = (
                    ref_val + delta_val
                    if not pd.isna(ref_val) and not pd.isna(delta_val)
                    else np.nan
                )
                bed_table.append(
                    {
                        "Element": name,
                        "REF": f"{ref_val:.4f}" if not pd.isna(ref_val) else "N/A",
                        "ALT": f"{alt_val:.4f}" if not pd.isna(alt_val) else "N/A",
                        "Delta": f"{delta_val:+.5f}"
                        if not pd.isna(delta_val)
                        else "N/A",
                    }
                )
        bed_table_df = pd.DataFrame(bed_table)

        # === 3. Sequence model metrics ===
        kl_max = row.get("MLM_KL_max", np.nan)
        logprob_ref = row.get("MLM_logprob_ref", np.nan)
        logprob_alt = row.get("MLM_logprob_alt", np.nan)
        ref_5mer = row.get("REF_5mer", "N/A") or "N/A"
        alt_5mer = row.get("ALT_5mer", "N/A") or "N/A"

        mlm_md = f"""
These metrics summarize how the NTv3 model responds to the alternate sequence relative to the reference:

- **Variant type:** {variant_type}{indel_label}  
- **LLR (Log-Likelihood Ratio):** {_fmt(llr)}{"  *(SNPs only — N/A for indels)*" if pd.isna(llr) else ""}  
- **KL Divergence (mean):** {_fmt(kl_mean)}  
  Average KL(ALT \u2225\u2225 REF) across the evaluated span  
- **KL Divergence (max):** {_fmt(kl_max)}  
- **Log-prob REF:** {_fmt(logprob_ref, fmt=".4f")}  
- **Log-prob ALT:** {_fmt(logprob_alt, fmt=".4f")}  
- **Context (REF):** `{ref_5mer}`  
- **Context (ALT):** `{alt_5mer}`  
"""

        if abs(indel_size) > 0:
            emb_cosine = row.get("EMB_cosine_dist", np.nan)
            emb_l2 = row.get("EMB_l2_dist", np.nan)
            emb_max = row.get("EMB_max_pos_dist", np.nan)
            emb_mean = row.get("EMB_mean_pos_dist", np.nan)
            mlm_md += f"""
**Embedding Distance Metrics (Indel-specific):**
- Cosine Distance: {_fmt(emb_cosine)}
- L2 Distance: {_fmt(emb_l2)}
- Max per-position dist: {_fmt(emb_max)}
- Mean per-position dist: {_fmt(emb_mean)}
"""

        # === 4. CSV Export ===
        try:
            result_csv_path = str(BASE_DIR / "last_variant_result.csv")
            results_df.to_csv(result_csv_path, index=False)
        except Exception as e:
            print(f"⚠️  CSV export failed: {e}")
            result_csv_path = None

        # === 5. Region Tracks Plot (uses unified ranking) ===
        try:
            region_tracks_fig = generate_region_tracks_plot(
                ranked_tracks=ranked,
                max_ranked_tracks=10,
                visible_radius_bp=max(int(zoom_bp), 8),
            )
        except Exception as e:
            print(f"⚠️  Region tracks plot failed: {e}")
            region_tracks_fig = None

        return (
            summary_md,
            interpretation_md,
            top_table_df,
            fingerprint_fig,
            region_tracks_fig,
            result_csv_path,
            bed_table_df,
            mlm_md,
            ranked,  # stashed for slider re-renders
        )

    except Exception as e:
        import traceback

        error_msg = f"❌ **Error during prediction:**\n\n```\n{str(e)}\n```\n\n**Traceback:**\n```\n{traceback.format_exc()}\n```"
        return (error_msg, None, None, None, None, None, None, None, None)


# ============================================================================
# BATCH PREDICTION
# ============================================================================
def predict_batch_wrapper(species, file):
    """Wrapper with ZeroGPU decorator if available."""
    return predict_batch(species, file)


# Apply batch decorator if available
if HAS_SPACES and spaces is not None:
    predict_batch_wrapper = spaces.GPU(predict_batch_wrapper)


def predict_batch(species, file):
    """
    Predict functional impact for multiple variants from CSV.

    Args:
        file: Uploaded CSV file with columns [chrom, pos, ref, alt]

    Returns:
        Tuple of (preview_df, download_file_path)
    """
    try:
        if file is None:
            return (pd.DataFrame([{"Error": "No file uploaded"}]), None)

        # Read CSV — handle both Gradio 3.x (file object) and 4.x (path string)
        file_path = file if isinstance(file, str) else file.name
        input_df = pd.read_csv(file_path)

        # Validate columns
        required = {"chrom", "pos", "ref", "alt"}
        missing = required - set(input_df.columns)
        if missing:
            return (pd.DataFrame([{"Error": f"Missing columns: {missing}"}]), None)

        # Limit batch size
        if len(input_df) > 10:
            return (
                pd.DataFrame(
                    [
                        {
                            "Error": f"Batch size limited to 10 variants (received {len(input_df)})"
                        }
                    ]
                ),
                None,
            )

        species = str(species).strip()
        print(f"🔬 Processing batch of {len(input_df)} variants for {species}...")

        # Run inference
        results_df = predict_variants(input_df, device=DEVICE, species=species)

        # Annotate
        results_df = _apply_species_annotations(results_df, species)

        # Compute impact scores
        results_df = compute_impact_scores(results_df)

        # Save full results
        output_path = BASE_DIR / "batch_results.csv"
        results_df.to_csv(output_path, index=False)

        # Create preview (summary columns)
        preview_df = format_summary_table(results_df)

        print(f"✅ Batch processing complete: {len(results_df)} variants")

        return (preview_df, str(output_path))

    except Exception as e:
        import traceback

        error_msg = f"Error: {str(e)}\n\n{traceback.format_exc()}"
        return (pd.DataFrame([{"Error": error_msg}]), None)


# ============================================================================
# GRADIO INTERFACE
# ============================================================================
def build_interface():
    """Build the Gradio interface."""

    # Custom CSS
    custom_css = """
    .variant-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 2rem;
        border-radius: 10px;
        margin-bottom: 2rem;
    }
    .metric-card {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #667eea;
    }
    """

    # Build status line from cached model info (safe — loads at startup)
    _seq_src = _MODEL_CACHE.get("sequence_source", "ucsc")
    _bed_n = len(_MODEL_CACHE.get("bed_names", []) or [])
    _bw_n = len(_MODEL_CACHE.get("selected_bw_indices", []) or [])
    _species_choices = [
        (f"{label} ({species})", species)
        for species, label in SUPPORTED_UI_SPECIES.items()
    ]

    if _APP_READY:
        _status_line = (
            "🟢 **Model ready** — "
            f"Device: `{DEVICE}` • "
            f"Sequence source: `{_seq_src}` • "
            f"Tracks: {_bed_n} BED + {_bw_n} BigWig"
        )
    else:
        _status_line = "🔴 **Model failed to load** — Check configuration and logs"

    # Gradio 6.x: create Blocks without theme/css, then set as attributes
    app = gr.Blocks(title="MAGI Variant Interpreter")
    app.theme = gr.themes.Soft()
    app.css = custom_css

    with app:
        gr.Markdown(
            """
        # 🧬 MAGI Variant Interpreter
        
        **Variant interpretation demo built on the NTv3 foundation model**
        
        MAGI combines NTv3 sequence predictions with annotation, ranking, and rule-based summaries to help review variant-associated signals in a local genomic window.

        The app reports changes across:
        - **Regulatory elements** (promoters, enhancers, CTCF sites)
        - **Structural features** (coding sequence, splice sites, UTRs)
        - **Epigenetic marks** (histone modifications, chromatin accessibility)
        - **Context-dependent tracks** (for example CAGE and chromatin assays)
        
        **Attribution:** MAGI was developed by Dan Ofer, Stav Zok, and Michal Linial. 
        
        **This web app** extends MAGI with multi-species support: 15+ animals and 6 plants can be scored via Ensembl REST APIs for on-the-fly sequence fetching. Human remains the primary supported species with full BigWig and MANE transcript annotation.
        """,
            elem_classes="variant-header",
        )
        gr.Markdown(_status_line)

        with gr.Tabs():
            # ===================================================================
            # TAB 1: SINGLE VARIANT
            # ===================================================================
            with gr.Tab("🔬 Single Variant Analysis"):
                gr.Markdown("""
                ### Manual Variant Input
                Enter a single variant for detailed scoring and interpretation.
                Human defaults to **GRCh38/hg38**. For non-human species, coordinates should match the selected species assembly available through Ensembl.
                """)

                with gr.Row():
                    with gr.Column(scale=1):
                        species_input = gr.Dropdown(
                            choices=_species_choices,
                            label="Species",
                            value="human",
                        )
                        chrom_input = gr.Textbox(
                            label="Chromosome",
                            value="chr17",
                            max_lines=1,
                        )
                        gr.Markdown(
                            "*Human: `chr1`–`chr22`, `chrX`, `chrY`, `chrM`. "
                            "Non-human: bare names (e.g., `1`, `X`, `MT`) or with `chr` prefix. "
                            "Coordinates should match your species' current Ensembl assembly.*"
                        )
                        pos_input = gr.Number(
                            label="Position (1-based)",
                            value=7675088,
                            precision=0,
                        )
                    with gr.Column(scale=1):
                        ref_input = gr.Textbox(
                            label="Reference Allele",
                            value="C",
                            max_lines=1,
                        )
                        alt_input = gr.Textbox(
                            label="Alternate Allele",
                            value="T",
                            max_lines=1,
                        )

                predict_btn = gr.Button(
                    "🚀 Predict Impact", variant="primary", size="lg"
                )

                with gr.Accordion("📋 Example Variants", open=False):
                    gr.Examples(
                        examples=[
                            [
                                "human",
                                "chr17",
                                7675088,
                                "C",
                                "T",
                            ],  # TP53 R175H — Pathogenic missense, most common cancer driver mutation
                            [
                                "human",
                                "chr7",
                                117559593,
                                "ATCT",
                                "A",
                            ],  # CFTR F508del — Pathogenic 3-bp deletion, ~70% of cystic fibrosis alleles
                            [
                                "human",
                                "chr13",
                                32332771,
                                "AGAGA",
                                "AGA",
                            ],  # BRCA2 c.5946delT — Pathogenic frameshift, hereditary breast/ovarian cancer
                            [
                                "human",
                                "chr11",
                                5227002,
                                "T",
                                "A",
                            ],  # HBB E6V (rs334) — Pathogenic missense, causes sickle cell disease
                            [
                                "human",
                                "chr17",
                                43092418,
                                "T",
                                "C",
                            ],  # BRCA1 c.3113A>G (rs16941) — Benign synonymous variant
                        ],
                        inputs=[
                            species_input,
                            chrom_input,
                            pos_input,
                            ref_input,
                            alt_input,
                        ],
                        label="Click to load pre-filled examples",
                    )

                gr.Markdown("---")
                gr.Markdown("### Results")

                with gr.Row():
                    with gr.Column(scale=1):
                        summary_output = gr.Markdown()
                    with gr.Column(scale=1):
                        interpretation_output = gr.Markdown()

                with gr.Row():
                    with gr.Column(scale=2):
                        top_table_output = gr.DataFrame(
                            label="Top ranked tracks (BED + BigWig, ordered by |Δ|)",
                            wrap=True,
                            max_height=600,
                        )
                    with gr.Column(scale=3):
                        fingerprint_output = gr.Plot(label="Top signals bar chart")

                with gr.Accordion("📥 Download Results", open=False):
                    download_single_output = gr.File(label="Download CSV")
                    download_tracks_output = gr.File(
                        label="Download Track PNG (current zoom/filter)"
                    )

                ranked_state = gr.State(value=None)

                with gr.Accordion("🔬 Region Track View", open=True):
                    gr.Markdown(
                        "*Reference (dashed gray) vs alternate signal for the top disrupted tracks. "
                        "Red shading = gain of function, blue shading = loss of function. "
                        "The slider maximum updates to the available track span for the current prediction.*"
                    )
                    with gr.Row(equal_height=True):
                        zoom_slider = gr.Slider(
                            minimum=MIN_ZOOM_BP,
                            maximum=INITIAL_ZOOM_MAX_BP,
                            value=DEFAULT_ZOOM_BP,
                            step=8,
                            label="Zoom radius (bp around variant)",
                            info="Updates after each prediction to the exact available track radius",
                            scale=3,
                        )
                        track_filter_dd = gr.Dropdown(
                            choices=["All", "BED", "BigWig"],
                            value="All",
                            label="Track type filter",
                            scale=1,
                        )
                    region_tracks_output = gr.Plot(label="Region Track View")

                with gr.Accordion("📊 All BED Elements", open=False):
                    bed_table_output = gr.DataFrame(
                        wrap=True,
                    )

                with gr.Accordion("📋 Sequence Model Metrics", open=True):
                    mlm_output = gr.Markdown()

                # Connect button
                predict_btn.click(
                    fn=predict_single_variant_wrapper,
                    inputs=[
                        species_input,
                        chrom_input,
                        pos_input,
                        ref_input,
                        alt_input,
                        zoom_slider,
                    ],
                    outputs=[
                        summary_output,
                        interpretation_output,
                        top_table_output,
                        fingerprint_output,
                        region_tracks_output,
                        download_single_output,
                        download_tracks_output,
                        bed_table_output,
                        mlm_output,
                        zoom_slider,
                        ranked_state,
                    ],
                )

                # Zoom slider re-renders region tracks on release (not during drag)
                zoom_slider.release(
                    fn=_update_region_zoom,
                    inputs=[zoom_slider, ranked_state, track_filter_dd],
                    outputs=[region_tracks_output, download_tracks_output],
                )
                track_filter_dd.change(
                    fn=_update_region_zoom,
                    inputs=[zoom_slider, ranked_state, track_filter_dd],
                    outputs=[region_tracks_output, download_tracks_output],
                )

            # ===================================================================
            # TAB 2: BATCH UPLOAD
            # ===================================================================
            with gr.Tab("📤 Batch Analysis"):
                gr.Markdown("""
                ### Batch Variant Upload
                Upload a CSV file with multiple variants for batch processing.
                
                **Required columns:** `chrom`, `pos`, `ref`, `alt`  
                **Format:** 1-based coordinates matching the selected species assembly  
                **Limit:** Maximum 10 variants per batch  
                **Chromosome naming:** For non-human, bare names (e.g., `1`, `X`, `MT`) or `chr`-prefixed both work.
                
                **Example CSV format (human):**
                ```
                chrom,pos,ref,alt
                chr17,7675088,C,T
                chr7,117559593,ATCT,A
                chr13,32340300,G,A
                ```
                
                **Example CSV format (non-human):**
                ```
                chrom,pos,ref,alt
                2,50000000,A,T
                X,25000000,C,G
                ```
                """)

                file_input = gr.File(
                    label="Upload CSV File",
                    file_types=[".csv"],
                )

                batch_species_input = gr.Dropdown(
                    choices=_species_choices,
                    label="Species",
                    value="human",
                )

                batch_btn = gr.Button("🚀 Process Batch", variant="primary", size="lg")

                gr.Markdown("---")
                gr.Markdown("### Results Preview")

                preview_output = gr.DataFrame(
                    label="Summary Table (first 10 rows)",
                    wrap=True,
                )

                download_output = gr.File(
                    label="Download Full Results (CSV)",
                )

                # Connect button
                batch_btn.click(
                    fn=predict_batch_wrapper,
                    inputs=[batch_species_input, file_input],
                    outputs=[preview_output, download_output],
                )

            # ===================================================================
            # TAB 3: DOCUMENTATION
            # ===================================================================
            with gr.Tab("📖 Documentation"):
                gr.Markdown("""
                ## About MAGI
                
                MAGI is a variant interpretation workflow that uses a Genomic foundation model (NTv3) for sequence scoring and adds:
                - gene and region annotation from MANE Select transcripts
                - ranking of the strongest BED and BigWig changes
                - a compact rule-based text summary
                - a baseline-derived MAGI score (`Global_z_sum_log`)

                ### Model configuration
                - Current model: `InstaDeepAI/NTv3_650M_post`
                - Current sequence window: **32 kb** (`CONTEXT_LEN = 32768`)
                - Region Track View zoom is capped by the available NTv3 track-profile span for the current prediction
                - Sequence source:
                  - **Human:** local `hg38.fa` when available, otherwise UCSC REST fallback
                  - **Non-human animals:** Ensembl REST API
                  - **Plants:** Ensembl Plants REST API (with Ensembl REST fallback)
                
                ### Output groups
                
                **BED outputs**
                - Functional and structural elements supplied by the NTv3 model configuration
                - Typically include coding, splice, promoter, UTR, exon, intron, and related annotations
                
                **BigWig outputs**
                - A filtered subset of assay tracks such as histone marks, chromatin accessibility, and CAGE
                - Used as contextual signals rather than direct mechanistic proof
                - Currently reported for human only in this app
                
                **Sequence model metrics**
                - `LLR` for SNPs
                - `KL` divergence summaries
                - log-probability differences
                - embedding distance summaries for indels
                
                ### Gene Annotation
                Variants are automatically annotated with:
                - **Gene name** (from MANE Select RefSeq transcripts)
                - **Region class:** CODING, SPLICE, UTR_5, UTR_3, PROMOTER, INTRONIC, GENIC_OTHER, OTHER
                - **Annotation flags:** overlap with coding sequence, splice sites, promoters, UTRs, and related transcript features
                - Non-human species skip MANE transcript annotation and report sequence- and BED-based outputs only
                - All supported species (animals and plants) can be scored via BED elements and MLM sequence-model features
                
                ### Impact Scoring
                - **MAGI score (`Global_z_sum_log`)**: baseline-derived burden score across BED and BigWig deltas
                - **Impact_Score_BED**: mean absolute value of the top 3 BED deltas
                - **Impact_Score_BW**: mean absolute value of the top 10 BigWig deltas
                
                Larger values indicate stronger deviation from the reference prediction. The summary tier shown in the Variant Summary is currently derived from `Impact_Score_BED` only.

                Positive Δ indicates an increase in the predicted signal; negative Δ indicates a decrease.

                ### Interpretation panel
                The rule-based interpretation panel summarizes the strongest ranked BED, BigWig, and sequence-model signals already shown elsewhere in the app. It is a compact heuristic summary, not a calibrated pathogenicity assessment.
                
                ### Limitations
                1. **Predictions are computational** and require experimental follow-up.
                2. No phasing information is used.
                3. **Coordinates must match the species assembly available through Ensembl:**
                   - Human: GRCh38/hg38
                   - Other animals: latest Ensembl assembly for that species
                   - Plants: latest Ensembl Plants assembly for that species
                4. BigWig context tracks are currently human-only in this app.
                5. MANE transcript annotation is human-only; non-human predictions show BED-level and sequence-model features only.
                6. Baseline z-scores (`Global_z_sum_log`, `magi_baseline_stats.csv`) are derived from human variants and may be less meaningful for non-human species. Use with caution.
                
                ### Citation
                If you use MAGI in your research, please cite the MAGI manuscript. If you find this webserver useful, please mention it! If you rely on the underlying foundation model, we suggest cite NTv3.
                
                ```
                Ofer, D., Zok, S., & Linial, M. (2026). MAGI: Mechanistic Interpretation of 
                Genetic Variants Consequences via Genomic Foundation Models.
                ```
                
                ### External model credit
                **NTv3:** Dalla-Torre, H., et al. (2025). Nucleotide Transformer: building and evaluating robust foundation models for human genomics. Nature Methods .
                
                <!-- ### Contact & Source Code
                - **Model:** [InstaDeepAI/NTv3_650M_post](https://huggingface.co/InstaDeepAI/NTv3_650M_post)
                - **Paper:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.01.11.523679v2)
                - **GitHub:** [Source repository](https://github.com/instadeepai/nucleotide-transformer) -->
                
                ---
""")

    return app


# Build app at module level for Gradio/Spaces auto-detection
app = build_interface()

if __name__ == "__main__":
    app.launch()
