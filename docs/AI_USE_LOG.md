# AI and Tool Assistance Log

This log supports transparent attribution. The student should add entries whenever AI materially influences code, analysis, interpretation, or communication.

| Date | Tool | Assistance | Files or decisions affected | Student verification |
|---|---|---|---|---|
| 2026-07-25 | OpenAI Codex | Audited the v1 repository; identified external validation and IRIS/ISEF compliance as v2 priorities; implemented and ran the initial E-MTAB-1980 pipeline integration; added proportional-hazards and FDR gates; created internal readiness and acceptance checklists. | `analysis/14_external_survival_emtab1980.R`, pipeline configuration, validation checks, E-MTAB-1980 result tables, internal docs | Automated output validation passed. Pending: student must inspect the code, explain each modeling choice, independently reproduce the run, and record the interpretation in the student's own words. |
| 2026-07-25 | OpenAI Codex | Implemented selection stability, repeated held-out clinical-increment testing, Human Protein Atlas cell-source triangulation, GDC schema compatibility, GEO download resilience, source provenance, and run checksums; installed R dependencies; ran and debugged the complete pipeline. | `analysis/00_config.R`, `analysis/01_download_tcga.R`, `analysis/02_download_geo.R`, `analysis/05_*`, `analysis/12_validate_outputs.R`, `analysis/13_*` through `analysis/18_*`, `run_pipeline.ps1`, generated result tables and figures | Final uninterrupted run printed `Pipeline complete.` after output validation. Student still must inspect and explain code, independently check interpretations, and disclose assistance. |

## Boundaries

- AI output must not be copied into the IRIS/ISEF research plan, abstract, poster, citations, or video script.
- AI-generated code or analysis advice must be disclosed as a project resource and understood by the student.
- The student remains responsible for scientific decisions, execution, interpretation, citation checking, and every claim presented to judges.
