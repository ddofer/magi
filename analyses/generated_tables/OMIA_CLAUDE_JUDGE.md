# OMIA judge results from Claude Sonnet

Two tables, one row per animal variant (dog, cat, chicken). Each was judged with the same prompts and system prompts as `clinvar_pipeline/scripts/10_animals_run_llm_parallel.py`, by Claude Sonnet instead of Gemini, on 1 November 2026. Gemini was not reachable from the cluster.

| File | Variants | Inputs | Concordant |
|---|---|---|---|
| `omia_snv_judge_claude_sonnet.csv` | 372 SNVs | the recorded run, the same inputs Gemini judged | 293 (79%) |
| `omia_indel_judge_claude_sonnet_real_alleles.csv` | 217 indels | a rerun with real alleles | 148 (68%) |

Per species, against the recorded Gemini run (Supplementary Fig. S1, read off the figure):

| Species | SNVs, Gemini | SNVs, Claude | Indels, Gemini, placeholder alleles | Indels, Claude, real alleles |
|---|---|---|---|---|
| Cat | 82% | 82% (103/125) | 17% (60 indels) | 70% (39/56) |
| Chicken | 70% | 59% (16/27) | 12% (8) | 70% (7/10) |
| Dog | 83% | 79% (174/220) | 31% (140) | 68% (102/151) |

On SNVs, Claude gives about the same concordance as Gemini, or slightly lower. So the rise in the indel bars comes from the real alleles, not from a more lenient judge. The recorded indel run passed `-` or `N` to the model in place of the inserted or deleted letters; `inference_omia.py` now builds VCF-style alleles from the species genome, and 44 OMIA indels still cannot be scored.

`rationale_names_mechanism` is false when the judge read the OMIA text as naming no mechanism (for example "accession only"). The system prompts count such a text plus weak signals as concordant. Where the text names a mechanism, 194 of 272 SNVs and 73 of 141 indels are concordant.
