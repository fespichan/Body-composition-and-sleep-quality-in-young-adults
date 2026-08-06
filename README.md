# Body composition and sleep quality in young adults

Analysis code for the cross-sectional study **"Biological and Anthropometric Correlates of Sleep
Quality in a Young-Adult University Sample"** (Universidad Nacional del Callao, Peru).

The design is **explanatory**: associations are reported as odds ratios with 95% confidence
intervals. The study does **not** develop or evaluate a prediction model.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21814202.svg)](https://doi.org/10.5281/zenodo.21814202)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Citation

Cite the archived code through its **concept DOI**, which always resolves to the most recent
version:

> Espichán, F., Carbajal, L., & Siccha Macassi, A. L. (2026). *Analysis code for: Biological and
> Anthropometric Correlates of Sleep Quality in a Young-Adult University Sample* [Computer
> software]. Zenodo. https://doi.org/10.5281/zenodo.21814202

To cite the exact state of the code rather than the latest, use the version DOI —
[`10.5281/zenodo.21814202`](https://doi.org/10.5281/zenodo.21814202) for v1.2.

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).


---

## Requirements

- **R 4.3.3**
- Package versions are pinned in [`renv.lock`](renv.lock). `renv::restore()` installs them.

| Package | Version | Used for |
|---|---|---|
| `mice` | 3.17.0 | multiple imputation |
| `ranger` | 0.17.0 | the `rf` engine behind `mice` |
| `car` | 3.1.2 | GVIF |
| `FactoMineR` | 2.11 | multiple factor analysis |
| `factoextra` | 1.0.7 | MFA plots |
| `psych` + `GPArotation` | 2.4.3 / 2026.4.1 | reliability (α, ω) |
| `dplyr`, `ggplot2`, `patchwork`, `here`, `readxl` | | data handling, figures, paths |

`renv.lock` also pins `pROC`, `caret`, `randomForest`, `themis` and `recipes`. **None of these is
needed to run the pipeline** — they are required only by the `audit_*.R` diagnostics described
below, which reproduce analyses that were tested and then removed from the study.

`renv.lock` was written directly from the verified session versions rather than produced by
`renv::snapshot()`, so it pins versions but carries no package hashes.

> `mice` changed its `rf` engine from `randomForest` to `ranger` at some point in its history.
> Under a different engine the RNG stream differs and `set.seed(123456)` will not reproduce the
> imputations bit-for-bit. Use the pinned version.

---

## Repository layout

```
R/                       analysis scripts, numbered in run order (NN_name.R)
  lib_*.R                libraries, sourced by the numbered scripts
  audit_*.R              diagnostics — NOT part of the run order, kept for provenance
scripts/                 QA runner (reads the private source workbook)
data/raw/                de-identified input data
data/processed/          intermediates (created on first run; gitignored)
outputs/reports/         run logs, named after the script that writes them (gitignored)
outputs/figures/         generated figures (gitignored)
```

Paths resolve through `here::here()`. There is no orchestrator and no setup script **by design**:
each file runs on its own and is meant to be read line by line as the record of what was computed.

> ⚠ **Run with the working directory inside the project.** Open the `.Rproj`, or `setwd()` to the
> project root first. `here::here()` resolves from the working directory, not from the script's own
> location, so `Rscript /elsewhere/R/01_clinical_audit.R` launched from an unrelated directory will
> anchor somewhere else entirely.

---

## Run order

Each script produces one manuscript object. Run them in numerical order.

| # | Script | Produces | Reads | Writes |
|---|---|---|---|---|
| 00 | `00_deidentify.R` | the public source file | private raw CSV | `data/raw/antro_anon.csv` — **data owner only**, see below |
| 01 | `01_clinical_audit.R` | the keystone | `data/raw/antro_anon.csv` | `data/processed/sleep_data_clean_for_imputation.csv` (289 × 16) |
| 02 | `02_multiple_imputation.R` | the imputations | the keystone | `data/processed/mice_imputation_m100.rds` (m = 100, maxit = 50, seed 123456) |
| 03 | `03_table1_descriptive.R` | **Table 1** | the keystone | log |
| 04 | `04_reliability.R` | **Table S1** | private workbook | console |
| 05 | `05_table2_global_model.R` | **Table 2 + Table S2** | the imputations | log |
| 06 | `06_mfa.R` | **Figures 3 and 4** | the imputations | `outputs/figures/*.tif` (1200 dpi) |
| 07 | `07_component_grids.R` | **Table 3** | the imputations | log |
| 08 | `08_figure1_flow.R` | **Figure 1** (participant flow) | the keystone | `.tif` (1200 dpi) |
| 09 | `09_figure2_forest.R` | **Figure 2** (forest plot) | the imputations | `.tif` (1200 dpi) |
| 10 | `10_supp_figures.R` | **Figures S1–S5** | the imputations | `.tif` (600 dpi; S5 at 1200) |

Step 02 takes about 5 minutes; every other step runs in seconds. Total ≈ 7 minutes.

**Step 00 cannot be run from this repository, and does not need to be.** It is the only script that
touches participant-identifying data. Its output, `data/raw/antro_anon.csv`, is committed here, so
the pipeline starts at step 01.

### Self-checking steps

Steps 03, 05, 08 and 09 **hard-assert every published value they produce** and abort if anything
fails to regenerate:

- step 03 checks all 91 cells of Table 1;
- step 05 checks the five coefficients, Nagelkerke R², EPV and all four GVIF values;
- step 08 checks the six participant counts against the observed-component pattern;
- step 09 re-derives the five plotted odds ratios and compares them against Table 2 **before
  drawing**, so the figure and the table cannot drift apart.

If the data or the model change, these scripts stop rather than emit a silently wrong number.

### Not part of the run order

- `R/lib_psqi_scoring.R`, `R/lib_reliability.R` — libraries, sourced by step 04 and the validator.
- `R/audit_verify_sweep.R` — verification sweep. Six of its seven blocks are superseded or belong
  to a section removed from the manuscript.
- `R/audit_congenial_pilot.R` — retired diagnostic (block vs congenial predictor matrix). It reads
  `master_clean.csv`, which is **not** the analysis input.
- `scripts/validate_psqi_scoring.R` — validates the scoring layer against the canonical workbook.

The `audit_*` scripts are kept deliberately. They are the record that circularity, uncongeniality,
imputed-label inflation and VIF-as-imputation-artifact were each tested and refuted, rather than
assumed away.

---

## Analysis summary

- **Primary analytic sample n = 206** — participants with at least one observed PSQI component
  (112 with all seven, 92 with six of seven, 2 with five of seven). The missing component is
  Component 4 (habitual sleep efficiency) in all 94 partial cases; it is non-computable when
  reported bed and rise times are ambiguous. The 83 participants who completed no PSQI item are
  excluded: their outcome would rest entirely on imputation with no outcome-informative data of
  their own. The full cohort of **n = 289** is reported as a sensitivity analysis.
- **Model** — poor sleep (PSQI > 5) ~ sex + age + BMI category + hemoglobin. Fitted on each of 100
  imputed datasets and pooled by Rubin's rules. Unpenalized; no variable selection was performed on
  the outcome.
- **Imputation** — MICE with two independent blocks: physiological variables do **not** predict PSQI
  components. Components are imputed from the other components plus sex (`pmm`); physiological
  variables use random forest. This prevents circularity by construction — imputed PSQI values carry
  no anthropometric information, so regressing PSQI on anthropometry is not circular.
- **Body fat** remains in the imputation as an auxiliary variable but was dropped as a *model
  predictor*: it duplicated the adiposity cluster already represented by BMI category, and was the
  source of the only GVIF above 3.
- **Determinism verified** — two runs under seed 123456 produced identical values across all
  105,700 imputed cells.

---

## Data availability and participant privacy

`data/raw/antro_anon.csv` is the de-identified input. Participant names are removed and identifiers
are replaced by positional labels (`P001`, `P002`, …). The `ide → anon_id` mapping is **never
written to disk**: `anon_id` is derived from row position, so re-running `00_deidentify.R`
reproduces the same labels without any key existing anywhere.

`00_deidentify.R` carries a leak guard: no column may contain free text other than `sex` and
`anon_id`. It also asserts that its output is cell-for-cell identical to the manually prepared
anonymized file (4624 cells, 0 differences).

The identifying source files are held privately by the corresponding author and are not
distributed. Raw data beyond the committed de-identified file are available on reasonable request,
subject to the terms of the ethics approval.

---

## Reproducibility status

**Verified from a clean checkout on 4 August 2026.** A tree containing only the committed files —
no private data, no intermediates, no outputs — ran steps 01–10 end to end and reproduced the
keystone byte-for-byte, all 105,700 imputed cells identically, every comparable log line, and all
nine figures by checksum.

| Manuscript output | Regenerates from committed code? |
|---|---|
| Table 1 (descriptives by sex) | ✅ `03_table1_descriptive.R` — all 91 cells asserted |
| Table 2 (four-predictor model, n = 206) | ✅ `05_table2_global_model.R` — every value asserted |
| Table 3 (component-level grid) | ✅ `07_component_grids.R` |
| Table S2 (GVIF) | ✅ `05_table2_global_model.R` |
| Figure 1 (participant flow) | ✅ `08_figure1_flow.R` |
| Figure 2 (forest plot) | ✅ `09_figure2_forest.R` |
| Figures 3–4 (MFA) | ✅ `06_mfa.R` |
| Figures S1–S4 (imputation diagnostics) | ✅ `10_supp_figures.R` |
| Figure S5 (MFA, Component 7) | ✅ `10_supp_figures.R` — writes `FigureS8_*.tif`, an obsolete filename for what the supplementary numbers **S5** |
| **Table S1 (reliability)** | ⚠ `04_reliability.R` runs, but reads the **private** workbook, so it is the one published object a reader cannot regenerate from this repository |

### Two conventions worth knowing

**Quantile definition (Table 1).** The published quartiles correspond to `stats::quantile`
**type 2** (the SAS definition), not R's default type 7. Verified empirically: type 2 reproduces
all 63 continuous cells, type 7 only 58. Two cells fall exactly on a half unit and are rounded away
from zero rather than to even. Step 03 documents and applies both conventions.

**Three-decimal confidence bounds.** `09_figure2_forest.R` prints a bound at three decimals whenever
two would round it to `1.00` — hemoglobin is reported as `1.43 (0.997–2.06)`, because `1.00–2.06`
would read as excluding the null when it does not. The rule is asserted to fire.

### A note on the legacy scripts

The original 14 scripts are retained outside this repository for provenance and are **not**
production code. Eleven of them require an in-memory object produced by another script and nine
perform no disk read at all, so they only run inside a live session in which an earlier script has
already been executed. They also use bare relative filenames and an implicit working directory. The
pipeline in `R/` is a clean rebuild, not a copy.

---

## Ethics

The study was reviewed by the Comité de Ética en Investigación, Vicerrectorado de Investigación,
Universidad Nacional del Callao (**Constancia N° 006-2026-CEI-VRI-UNAC**, 7 July 2026). All
participants gave informed consent.

---

## License

[MIT](LICENSE) © 2026 Fabio Francisco Espichán J.
