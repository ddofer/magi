# Validation checks, October 2026

Checks of whether MAGI's explanations agree with lab assays and with established tools. They are not a ranking of tools. Each folder holds the summary (`results.json`) and the main per-variant tables. The scripts and the large MAGI delta files are not included.

| Folder | Check | Headline |
|---|---|---|
| `E1_mfass` | Splicing assay (MFASS, about 28,000 variants measured for exon skipping) | Among variants SpliceAI calls harmless and VEP does not flag, 13 of MAGI's 42 splice calls truly disrupt splicing (31%; 2% by chance). Six of the 13 sit in one exon. |
| `E2_mpra` | Promoter and enhancer assay (4,141 variants) | On MAGI's 403 calls, the direction (element lost or gained) is right 66% of the time; always guessing "loss" gets 71%. |
| `E3E4` | Animal splice variants (OMIA), and 62 ClinVar variants VEP calls missense where MAGI calls splicing | 27 of 49 animal splice variants get a splice rule, against 2 of 225 missense. 25 of the 62 ClinVar variants have lab RNA evidence of a splicing effect; SpliceAI flags all 25. |
| `E5_vus` | ClinVar variants of uncertain significance with no protein change (13,536) | 95% of MAGI's splice calls are also flagged by SpliceAI; `magi_only_vus.csv` lists the rest. |
| `E6_omia_rerun` | OMIA rerun with real indel alleles | SNV calls unchanged. 217 indels scored, 88 get a call; `skipped_indels.csv` lists the 44 that cannot be scored. |
| `E7_subtype` | Are MAGI's named splice mechanisms right? 118 ClinVar variants with a lab splicing outcome | The broad call holds: a specific splice call for 72, and none for the 18 with normal RNA. The specific name matches the lab outcome no more often than shuffled names (68% against 74%). The per-variant lab labels are held back until a human check of them is done. |

The judge results for the OMIA variants are in `../generated_tables/OMIA_CLAUDE_JUDGE.md`.

## Third-party annotations

`comparators/annotations.parquet` has VEP 116, AlphaGenome variant-impact scores and precomputed SpliceAI scores for the variants above (46,060 rows, hg38). Column `set` names the check.

- VEP: vepyr 0.10.0 with the Ensembl VEP 116 cache from Hugging Face (`biodatageeks/vepyr_116_GRCh38_ensembl`), without the known-variant files.
- AlphaGenome scores: `just-dna-seq/alphagenome_avi`. **These fall under the AlphaGenome Output Terms of Use (http://deepmind.google.com/science/alphagenome/output-terms): non-commercial use only, and not for training variant-effect models.**
- SpliceAI: `luoyiming1991/spliceai-mane-v1.5-d500-m1-snv` (SpliceAI 1.3.1, MANE Select, distance 500, masked). **CC BY-NC 4.0.**

So this file is for non-commercial research only, whatever the licence of the rest of the repository.
