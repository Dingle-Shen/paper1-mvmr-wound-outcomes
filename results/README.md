# Expected outputs

`results/` is populated by the scripts. It is mostly excluded from the repository (`.gitignore`) because the output tables contain FinnGen summary-statistic extracts whose redistribution is governed by FinnGen terms; the four *input* tables listed under "Intermediate inputs" below are shipped with the repository (explicitly excepted from the ignore rule), as they are needed to run the pipeline without first rebuilding the instrument list.

## Script to output mapping

| Script | Output(s) | Manuscript location |
|---|---|---|
| `run_mvmr_main.py` | `results/mvmr_docA_r9local_summary.csv` (15 rows) | Table 2, Table 3, Table 5 (Q) |
| `fix_singlevar_estimates.py` | `results/mvmr_docA_r9local_summary_corrected.csv` (single-exposure column recomputed) | Table 2 single-exposure column, Figure 2 |
| `run_pleiotropy_sensitivity.py` | `results/mvmr_pleiotropy_egger.tsv`, `results/mvmr_outlier_attribution.tsv`, `results/mvmr_singlevar_audit.csv` | Table 5, Section 3.5, Section 4.2 |
| `run_mvmr_r12_sensitivity.py` | `results/mvmr_docA_r12local_summary.csv`, `_withq.csv` | Table 4 |
| `steiger_direction_test.py` | `results/steiger_direction_bmi_crp.tsv` | Supplementary Table S2 |
| `run_mediation_forward.py` | `results/mediation_network_direction.tsv`, `results/mediation_results_5outcomes.tsv` | Supplementary Tables S1 and S3 |
| `sensitivity_mvmr_b.py` | `results/sensitivity_mediation.tsv` | Supplementary Table S4 |
| `run_mediation_reverse.py` | `results/mediation_reverse.tsv`, `results/mediation_reverse_diagnosis.tsv` | Supplementary Tables S10–S12 (reverse-mediation contrast, Section S14) |
| `make_figures_p1.py` | `figures/fig1_forest_mvmr.png`, `figures/fig2_direction_flip.png` (plus PDF) | Figures 1-2 |

All tabular outputs are tab-separated except `mvmr_singlevar_audit.csv`, which is comma-separated.

Note on numbering: the supplementary material uses **two separate series**. Supplementary Methods are numbered `S1`–`S4` (in `Supplementary_Methods.docx`), while supplementary tables are numbered `Table S1`–`Table S12` (in `Supplementary_Results.docx`). The two series are offset by four, so "S2" alone is ambiguous: as a table it is the Steiger result, as a method it is the Steiger procedure. This README always writes `Table Sn` for a table and `Section Sn` for a method or result section.

Section-to-table mapping in `Supplementary_Results.docx`:

| Result section | Table | Content |
|---|---|---|
| S5 | Table S1 | Pairwise network MR, six directions |
| S6 | Table S2 | Steiger directionality test |
| S7 | Table S3 | Mediation results |
| S8 | Table S4 | Mediation b-path sensitivity |
| S9 | Table S5 | Outlier attribution |
| S10 | Table S6 | R9 recompute consistency |
| S11 | Table S7 | Full main-analysis estimates |
| S12 | Table S8 | Instrument list (352 variants) |
| S13 | Table S9 | Leave-one-out summary |
| S14 | Tables S10–S12 | Reverse-mediation contrast and decomposition |

## Intermediate inputs

The following input tables are shipped with the repository (excepted from the `results/*` ignore rule; see the main `README.md`):

- `results/mvmr_docA_ivs_rebuilt.tsv` — pre-specified pooled instrument list (492 variants; rsID / chr / pos). Needed by `run_mvmr_main.py` and `run_mvmr_r12_sensitivity.py`.
- `results/mediation_exposure_raw.tsv` — 352 instruments x three exposures, with effect alleles.
- `results/mediation_ivs_5outcomes.tsv` — 352 instruments x five outcomes (1760 rows), with outcome effects and per-outcome exposure effects.
- `results/_outcome_alleles_cache.tsv` — outcome reference/alternate alleles for the 352 instruments (cache; regenerated from FinnGen files when absent).

`run_pleiotropy_sensitivity.py` reads `mediation_ivs_5outcomes.tsv` because the per-outcome instrument file used in the main pipeline omits the chronic skin ulcer endpoint.

## Reference numbers (R9 baseline, reproduced with this code)

Multivariable estimates (beta, BH-adjusted q < 0.05 for 10 of 15 combinations):

| Outcome | T2D | BMI | CRP |
|---|---|---|---|
| Skin infection | 0.070 | 0.328 | 0.113 |
| Cellulitis | 0.088 | 0.267 | not significant |
| Lower limb ulcer | 0.215 | 0.470 | not significant |
| Chronic skin ulcer | 0.105 | 0.687 | not significant |
| Pressure ulcer | 0.232 | not significant | not significant |

- Conditional F: 26 (T2D), 21 (BMI), 44 (CRP).
- MVMR-Egger intercepts: all non-significant (P = 0.22, 0.24, 0.85, 0.79, 0.46).
- Leave-one-out: maximum absolute change in any estimate after removing a single instrument <= 0.075.
- Global Q significant for skin infection (P = 5.6e-4) and lower limb ulcer (P = 0.020); removing the top five Q-contributing variants renders both non-significant while the main effects are unchanged.
- Single-exposure versus multivariable: directions agree in all 15 combinations; Pearson r = 0.9245.
- Steiger (BMI vs CRP): both directions pass, consistent with shared genetic architecture rather than a single causal direction.
- Mediation (BMI -> CRP -> outcome), multivariable-adjusted b: proportions 6.4% to 9.3% in the three nominally significant outcomes.
- Reverse mediation (CRP -> BMI -> outcome): reverse proportions 4.6 to 58.9 times the forward proportions, entirely attributable to exposure effect-strength asymmetry (a-ratio 0.7476 <= 1 in every outcome); no outcome meets the criterion for re-examining the assumed direction.

## Verification steps and how to check them

| Manuscript item | How it is reproduced |
|---|---|
| Supplementary Table S6 (R9 original versus recomputed estimates) | Re-run `code/run_mvmr_main.py`. The pipeline is deterministic given the same input files, so the 15 multivariable estimates and their P values reproduce the published values exactly. Compare against the reference table below. |
| Section 3.3 / Section 4.2 single-exposure estimator correction | `code/run_pleiotropy_sensitivity.py` writes `mvmr_singlevar_audit.csv`, which contains both the pre-correction column (outcome-side weights only) and the standard inverse-variance column. The direction-reversal flags (`翻转_现有` versus `翻转_标准`) show 9 of 15 reversals before the correction and 0 of 15 after. |
| Section 3.6 R9 versus R12 direction consistency | `code/run_mvmr_r12_sensitivity.py` followed by comparison of the two summary tables. All 15 estimates keep their sign; the number of BH-significant combinations rises from 10 to 13. |
| Section 2.8 / Table 3 minimum detectable effects | `code/run_mde_power.py` writes `results/mde_power.tsv` and prints the OR matrix of Table 3. |
