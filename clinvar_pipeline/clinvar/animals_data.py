"""Load OMIA animal LLM concordance evaluation tables for supplementary figures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from config import PACKAGE_ROOT, PATHS


def _species_lookup_candidates() -> list[Path]:
    keys = (
        "omia_all_species_inferred_full",
        "omia_all_species_inferred_relevant",
    )
    paths: list[Path] = []
    for key in keys:
        path_str = PATHS.get(key)
        if path_str:
            paths.append(Path(path_str))
    paths.append(PACKAGE_ROOT / "notebooks" / "parquet" / "omia_all_species_inferred_full.parquet")
    return paths


def _load_species_lookup_frame() -> pd.DataFrame:
    """Load table used for ``#VariationID`` → ``species_key`` (notebook cell)."""
    for path in _species_lookup_candidates():
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        if "#VariationID" in df.columns:
            return df.set_index("#VariationID", drop=False)
        if df.index.name == "#VariationID":
            return df.reset_index(drop=False).set_index("#VariationID", drop=False)
        if "species_key" in df.columns:
            return df.reset_index(drop=True)
    raise FileNotFoundError(
        "Missing OMIA species lookup parquet (omia_all_species_inferred_full.parquet with "
        "species_key). Copy it to data/parquet/ or notebooks/parquet/ — see "
        "animals_result_analysis.ipynb."
    )


def _species_from_deltas_by_row_order() -> tuple[pd.Series, pd.Series]:
    """Fallback when OMIA parquet is absent: eval rows match delta row order."""
    snp_d = pd.read_parquet(PATHS["animals_snp_deltas"], columns=["species_key"])
    indel_d = pd.read_parquet(PATHS["animals_indel_deltas"], columns=["species_key"])
    snp_eval = pd.read_parquet(PATHS["animals_snp_eval"], columns=["#VariationID"])
    indel_eval = pd.read_parquet(PATHS["animals_indel_eval"], columns=["#VariationID"])
    if len(snp_eval) != len(snp_d) or len(indel_eval) != len(indel_d):
        raise ValueError("Animal eval row counts do not match delta parquets for species merge")
    return snp_d["species_key"].reset_index(drop=True), indel_d["species_key"].reset_index(drop=True)


def _attach_species_keys(res: pd.DataFrame) -> pd.DataFrame:
    try:
        lookup = _load_species_lookup_frame()
        res = res.copy()
        res["species_key"] = _map_species_key(res, lookup)
        return res
    except FileNotFoundError:
        snp_species, indel_species = _species_from_deltas_by_row_order()
        res = res.copy()
        res.loc[res["variant_type"] == "SNP", "species_key"] = snp_species.values
        res.loc[res["variant_type"] == "INDEL", "species_key"] = indel_species.values
        return res


def _map_species_key(eval_df: pd.DataFrame, lookup: pd.DataFrame) -> pd.Series:
    """Match ``animals_result_analysis.ipynb``: ``#VariationID``.map(lookup['species_key'])."""
    if lookup.index.name == "#VariationID" or "#VariationID" in lookup.index.names:
        species = lookup["species_key"]
        return eval_df["#VariationID"].map(species)
    # Notebook uses positional index in the OMIA parquet (default RangeIndex).
    species_by_pos = lookup.reset_index(drop=True)["species_key"]
    return eval_df["#VariationID"].map(species_by_pos)


def load_animals_eval_merged() -> pd.DataFrame:
    """Return SNP + indel LLM evaluation rows (notebook: animals_result_analysis.ipynb)."""
    frames: list[pd.DataFrame] = []
    for variant_type, path_key in (
        ("SNP", "animals_snp_eval"),
        ("INDEL", "animals_indel_eval"),
    ):
        path = Path(PATHS[path_key])
        if not path.exists():
            raise FileNotFoundError(f"Missing animal evaluation parquet: {path}")
        df = pd.read_parquet(path).copy()
        df["variant_type"] = variant_type
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True, sort=False)
    merged["concordance"] = merged["concordance"].astype(str).str.strip().str.upper()
    return merged


def load_animals_eval_res() -> pd.DataFrame:
    """Combined eval + ``species_key``, matching ``animals_result_analysis.ipynb``."""
    res = load_animals_eval_merged()
    res = _attach_species_keys(res)
    if res["species_key"].isna().any():
        missing = int(res["species_key"].isna().sum())
        raise ValueError(f"{missing} evaluation rows missing species_key after species merge")
    res["is_correct"] = res["concordance"].eq("CONCORDANT")
    return res


def per_species_accuracy_long(res: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per species × variant_type accuracy table (notebook ``by`` frame)."""
    if res is None:
        res = load_animals_eval_res()

    by = (
        res.groupby(["species_key", "variant_type"])
        .agg(n=("is_correct", "size"), accuracy=("is_correct", "mean"))
        .reset_index()
    )
    by["accuracy_pct"] = by["accuracy"] * 100
    return by.sort_values(["species_key", "variant_type"]).reset_index(drop=True)


def per_species_accuracy_table(res: pd.DataFrame | None = None) -> pd.DataFrame:
    """Wide pivot for CSV sidecar (species × SNP/INDEL accuracy and counts)."""
    by = per_species_accuracy_long(res)
    acc = by.pivot(index="species_key", columns="variant_type", values="accuracy_pct").fillna(0)
    n_piv = by.pivot(index="species_key", columns="variant_type", values="n").fillna(0).astype(int)

    rows: list[dict[str, object]] = []
    for species in acc.index:
        rows.append(
            {
                "species": species,
                "n_snp": int(n_piv.loc[species].get("SNP", 0)),
                "n_indel": int(n_piv.loc[species].get("INDEL", 0)),
                "snp_accuracy_pct": float(acc.loc[species].get("SNP", 0)),
                "indel_accuracy_pct": float(acc.loc[species].get("INDEL", 0)),
            }
        )
    return pd.DataFrame(rows)


def concordance_counts_by_variant_type(df: pd.DataFrame | None = None) -> pd.DataFrame:
    if df is None:
        df = load_animals_eval_res()
    counts = (
        df.groupby(["variant_type", "concordance"], dropna=False)
        .size()
        .rename("count")
        .reset_index()
    )
    totals = counts.groupby("variant_type")["count"].transform("sum")
    counts["fraction"] = counts["count"] / totals
    return counts.sort_values(["variant_type", "concordance"]).reset_index(drop=True)
