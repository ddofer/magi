"""Build MANE_processed.csv and Promoter_processed.csv from MANE RefSeq GFF.

Ported from ``MANE_processing.ipynb``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

GFF_COLUMNS = [
    "Chromosome",
    "Source",
    "Feature",
    "Start",
    "End",
    "Score",
    "Strand",
    "Frame",
    "Attributes",
]

ATTR_FIELDS = [
    "ID",
    "Parent",
    "Dbxref",
    "Name",
    "description",
    "gbkey",
    "gene",
    "gene_biotype",
    "product",
    "tag",
    "transcript_id",
]

# Keys match MANE_processing.ipynb (GFF uses lnc_RNA; only mRNA matches in practice)
PROMOTER_DISTANCE = {
    "mRNA": 2000,
    "lncRNA": 2000,
    "snoRNA": 500,
    "snRNA": 500,
    "telomerase_RNA": 1500,
    "antisense_RNA": 1500,
}


def parse_attributes(attribute_string: str) -> dict:
    attributes: dict[str, str] = {}
    for field in ATTR_FIELDS:
        match = re.search(rf"{field}=([^;]+)", attribute_string)
        if match:
            attributes[field] = match.group(1)
    if "Dbxref" in attributes:
        for value in attributes["Dbxref"].split(","):
            if value.startswith("Ensembl:"):
                attributes["Ensembl"] = value.split(":", 1)[1]
                break
    return attributes


def _calculate_promoter(row: pd.Series) -> tuple[int, int]:
    dist = PROMOTER_DISTANCE[row["Feature"]]
    if row["Strand"] == "+":
        tss = int(row["Start"])
        return tss - dist, tss
    tss = int(row["End"])
    return tss, tss + dist


def build_mane_tables(gff_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse MANE GFF and return (MANE, Promoter) dataframes."""
    gff_path = Path(gff_path)
    mane = pd.read_csv(
        gff_path,
        sep="\t",
        compression="gzip" if gff_path.suffix == ".gz" else None,
        comment="#",
        header=None,
        names=GFF_COLUMNS,
    )
    parsed = mane["Attributes"].apply(parse_attributes)
    mane = pd.concat([mane, pd.DataFrame(parsed.tolist())], axis=1)

    promoter = mane[mane["Feature"].isin(PROMOTER_DISTANCE.keys())].copy()
    promoter[["Promoter_Start", "Promoter_End"]] = promoter.apply(
        _calculate_promoter, axis=1, result_type="expand"
    )
    promoter = promoter.reset_index(drop=True)
    return mane, promoter


def write_mane_metadata(
    gff_path: str | Path,
    mane_out: str | Path,
    promoter_out: str | Path,
) -> None:
    mane, promoter = build_mane_tables(gff_path)
    Path(mane_out).parent.mkdir(parents=True, exist_ok=True)
    mane.to_csv(mane_out, index=False)
    promoter.to_csv(promoter_out, index=False)
    print(f"  Saved → {mane_out}  ({len(mane):,} GFF rows)")
    print(f"  Saved → {promoter_out}  ({len(promoter):,} promoter rows)")
