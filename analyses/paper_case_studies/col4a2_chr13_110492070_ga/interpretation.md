### Rule-Based Signal Interpretation

**Primary hypothesis:** Localized structural disruption  
**Context anchor:** COL4A2 / CODING  
**Why this is suggested:** The strongest BED signal points to exon.

**Top evidence**
- BED `exon` shows a strong loss in exon (REF 0.875 → ALT 0.520, Δ=-0.355).
- BED `intron` shows a strong gain in intron (REF 0.551 → ALT 0.871, Δ=+0.320).
- Context track `brain (RNA-seq [1])` has a strong loss (REF 1.000 → ALT 0.590, Δ=-0.410).
- Context track `liver (RNA-seq [1])` has a strong loss (REF 0.965 → ALT 0.570, Δ=-0.395).
- Sequence-model evidence is supportive: LLR -1.849 penalizes the alternate sequence.
