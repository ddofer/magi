# Vendored dependencies

Standalone copies of scripts that previously lived outside `clinvar_pipeline/`.

| File | Used by |
|------|---------|
| `nt_mechanism_assignment.py` | `scripts/07_assign_mechanisms.py` (fig3e mechanism concordance) |

If `vendor/nt_mechanism_assignment.py` is missing, stage 07 falls back to `nt_mechanism_assignment.py` at the repository root (parent of `clinvar_pipeline/`).
