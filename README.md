# Multivariable Mendelian randomization: metabolic exposures and chronic skin wound outcomes

Analysis code for the paper: *Independent causal effects of metabolic factors on chronic skin wound outcomes: a multivariable Mendelian randomization study* (paper 1 of the project).

All scripts reproduce the numbers reported in the manuscript tables and figures from publicly available summary statistics. Scripts were reorganised from the original project workspace; statistical logic is unchanged except for (a) hard-coded absolute paths replaced with repository-relative paths, and (b) FinnGen data directories made configurable through environment variables.

## Repository structure

```
code/
  mvmr.py                          Core MVMR routines: multivariable IVW, conditional F.
  ogw.py                           OpenGWAS client wrapper (exposure summary statistics).
  scipy_shim.py                    Minimal scipy compatibility layer (norm.sf/cdf/ppf, chi2.sf).
                                   Imported by scripts so the analysis runs without scipy installed.
  run_mvmr_main.py                 Main analysis: MVMR of T2D / BMI / CRP on the five FinnGen R9
                                   wound outcomes (352 instruments), conditional F, Cochran Q,
                                   BH-FDR over the 15-test family.
                                   -> Table 2, Table 5 (Q), Figure 1.
  run_mvmr_r12_sensitivity.py      Same pipeline on FinnGen R12 (version robustness).
                                   -> Table 4.
  fix_singlevar_estimates.py       Recomputes single-exposure estimates with the standard IVW
                                   estimator (weights w_j = w * beta^2) and writes an audited copy
                                   of the summary table.
                                   -> Table 2 single-exposure column, Figure 2.
  run_pleiotropy_sensitivity.py    MVMR-Egger intercept, Cochran Q with top-5 outlier attribution,
                                   single-exposure estimator audit.
                                   -> Table 5, Section 3.5.
  run_mde_power.py                 Minimum detectable effect per outcome at 80% power
                                   (OR = exp(2.801585 * SE)).
                                   -> Table 3.
  steiger_direction_test.py        Steiger directionality test for BMI and CRP (approximate and
                                   exact R^2 estimators).
                                   -> Supplementary Table S2.
  run_mediation_forward.py         Pairwise network MR (six directions) and two-step mediation
                                   BMI -> CRP -> outcome.
                                   -> Supplementary Tables S1 and S3.
  run_mediation_reverse.py         Reverse mediation CRP -> BMI -> outcome, read as a direction
                                   sensitivity analysis with a pre-specified three-component
                                   decomposition.
                                   -> Supplementary Table (mediation sensitivity).
  sensitivity_mvmr_b.py            Mediation b-path under two estimator conventions
                                   (multivariable-adjusted versus single-exposure).
                                   -> Supplementary Table S4.
  make_figures_p1.py               Figures 1-2 (English, Arial).
data/README.md                     How to obtain exposure and outcome data.
results/README.md                  Expected outputs; script-to-table mapping; FinnGen terms.
```

`results/` additionally ships four input tables that the scripts read (they are kept in `results/` because the scripts resolve their inputs there; they are excluded from the `results/*` ignore rule):

| File | Contents | Needed by |
|---|---|---|
| `results/mvmr_docA_ivs_rebuilt.tsv` | Pre-specified pooled instrument list (492 variants, rsID / chr / pos) | `run_mvmr_main.py`, `run_mvmr_r12_sensitivity.py` |
| `results/mediation_exposure_raw.tsv` | Per-variant effects of the three exposures (beta / se / ea / nea) on the 352 shared instruments, from OpenGWAS | `steiger_direction_test.py`, `run_mediation_forward.py`, `sensitivity_mvmr_b.py`, `run_mediation_reverse.py` |
| `results/mediation_ivs_5outcomes.tsv` | Harmonised exposure and outcome effects for 352 instruments x 5 outcomes | `run_mediation_forward.py`, `sensitivity_mvmr_b.py`, `run_pleiotropy_sensitivity.py`, `fix_singlevar_estimates.py`, `run_mediation_reverse.py` |
| `results/_outcome_alleles_cache.tsv` | Outcome reference/alternate alleles for the 352 instruments (ref / alt per outcome) | `run_mediation_forward.py`, `run_mediation_reverse.py`, `sensitivity_mvmr_b.py` |

These tables are derived from the public summary statistics named in `data/README.md`; their redistribution is covered under the OpenGWAS API terms and the FinnGen CC BY 4.0 licence.

## Reproduction order

Run all commands from the repository root. Scripts locate the root as the parent of `code/`.

1. `python code/run_mvmr_main.py` — regresses the five outcomes on the three exposures. Requires the FinnGen R9 files under `results/finngen_full/` (see `data/README.md`). Writes the main summary table.
2. `python code/fix_singlevar_estimates.py` — replaces the single-exposure column with standard IVW estimates and writes `mvmr_docA_r9local_summary_corrected.csv`.
3. `python code/run_pleiotropy_sensitivity.py` — MVMR-Egger intercepts, heterogeneity Q with outlier attribution, and the single-exposure audit (before/after the estimator correction).
4. `python code/run_mde_power.py` — minimum detectable effects and the Table 3 OR matrix.
5. `python code/run_mvmr_r12_sensitivity.py` — reruns the pipeline on FinnGen R12. Requires `results/finngen_R12/` or `MR_FINNGEN_R12_DIR`.
6. `python code/steiger_direction_test.py` — Steiger directionality test.
7. `python code/run_mediation_forward.py` — pairwise network MR and forward mediation.
8. `python code/sensitivity_mvmr_b.py` — mediation b-path sensitivity.
9. `python code/run_mediation_reverse.py` — reverse mediation direction-sensitivity analysis.
10. `python code/make_figures_p1.py` — writes the figures.

For a reviewer who wants to check one specific number, `results/README.md` maps every manuscript table and supplementary table to the script that produces it, and gives reference values to compare against.

`scipy` is optional: `scipy_shim.py` provides the small subset of `scipy.stats` used here and is imported at the top of each script that needs it. If scipy is installed it is used as normal.

## Analysis choices and notes (as actually implemented)

- **Instruments**: genome-wide significant variants (P < 5e-8) per exposure with F > 10, pooled across the three exposures into one set, then LD-pruned between exposures (10 Mb window, r^2 < 0.001, 1000 Genomes European panel). 492 variants are pooled and 352 remain after harmonisation with the outcome data.
- **Single-exposure contrast**: computed on the same 352 instruments, with standard inverse-variance weights `w_j = w * beta_j^2`. An earlier version of the code used outcome-side weights only, which put the point estimate and its standard error on inconsistent weightings and produced spurious direction reversals in 9 of 15 combinations. The corrected estimator yields direction agreement in all 15 combinations (Pearson r = 0.9245 between the two estimate sets).
- **Conditional F**: residual-based formulation following Sanderson et al. (2019); values 26 / 21 / 44 for T2D / BMI / CRP.
- **FDR**: one family of 15 tests (3 exposures x 5 outcomes), Benjamini-Hochberg step-up with cumulative minimum; threshold q < 0.05.
- **MVMR-Egger**: intercept estimated in a weighted regression that includes the three exposures plus an intercept, weights `1/se_out^2`. Q statistic degrees of freedom = 349.
- **Outlier attribution**: the five variants contributing most to Q are removed and the model refitted.
- **Mediation**: two-step product-of-coefficients. `a` from exposure instruments, `b` from either multivariable-adjusted or single-exposure instruments, `c` as the exposure total effect. Delta-method standard error for `ab`; proportion = `ab / c`. The BMI -> CRP coefficient `a = 0.3378` (SE 0.0087) used by `sensitivity_mvmr_b.py` comes from the pairwise network MR estimate of `run_mediation_forward.py` (written to `results/mediation_network_direction.tsv`).
- **Reverse mediation**: reported as a direction *sensitivity* analysis, not as direction confirmation. The proportion ratio decomposes exactly as `(a_rev/a_fwd) * (b_rev/b_fwd) * (c_fwd/c_rev)`; because the `a` ratio is 0.75 (<= 1) in every outcome, the higher reverse proportion is attributable to asymmetry in exposure effect strength rather than to directional evidence. See `run_mediation_reverse.py` for the implementation and the three-tier interpretation rule.

## Environment

- Python 3.10+ with `numpy` and `pandas`; `scipy` and `matplotlib` recommended (see `requirements.txt`).
- No R is required for the analyses in this repository.

## License

Code: MIT (see `LICENSE`). Data belong to their respective providers (OpenGWAS, FinnGen) and are subject to their terms. The repository ships only the small derived tables listed above (instrument list and harmonised effects on 352 variants), which are redistributable under the OpenGWAS API terms and the FinnGen CC BY 4.0 licence; no full summary-statistics rows are redistributed.

## Manuscript code-availability statement (suggested)

> Analysis code is available at [repository URL / DOI]. Public summary statistics from OpenGWAS and FinnGen (R9 and R12) were used; see the data availability statement for access procedures.
