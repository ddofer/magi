PROVENANCE ONLY - these scripts are NOT runnable as-is.

They are the one-off helpers used to build fig5b (Ensembl SVG export -> row
extraction) and fig5c (Ensembl REST -> recomputed alignment) in July 2026.
Several of them hardcode a session scratch directory (`SD = "/tmp/claude-.../
scratchpad"`) that no longer exists, and they expect intermediate downloads that
were not committed.

They are kept so the method is auditable: what was fetched, from which endpoint,
and how the artwork was assembled. The authoritative description is
../fig5bc_PROVENANCE.md, which records the exact URLs, Ensembl release, POST
payloads and coordinates.

DO NOT try to regenerate fig5b/fig5c from these. The finished vector files
(../fig5b_MAT1A_transcripts.{pdf,svg}, ../fig5c_MAT1A_transcript_alignment.{pdf,svg})
are the deliverables, and at least part of the pipeline can no longer be
reproduced anyway: the Ensembl 115 archive - the release the published panel was
made on - now returns HTTP 403 to automated clients.

If these panels genuinely need rebuilding, start from fig5bc_PROVENANCE.md and
drive the current Ensembl release by hand, then re-check the transcript IDs and
aa lengths against the published panel.
