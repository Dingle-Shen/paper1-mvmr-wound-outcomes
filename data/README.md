# Data acquisition

The repository ships only small derived tables (the pre-specified instrument list and harmonised effects on 352 variants, see `README.md`); the full summary statistics required to rerun the pipeline from raw data must be obtained as below.

## Exposure summary statistics (OpenGWAS)

Fetched by `code/ogw.py` through the OpenGWAS API. An OpenGWAS token is required; set it in the environment as expected by your `ogw` configuration.

| Exposure | OpenGWAS ID | Sample size |
|---|---|---|
| Type 2 diabetes | `ebi-a-GCST006867` | 62,892 cases / 596,424 controls |
| Body mass index | `ebi-a-GCST006368` | 315,347 |
| C-reactive protein | `ebi-a-GCST90029070` | 575,531 |

## FinnGen R9 (baseline, required)

Five endpoints, downloaded from the FinnGen data portal (https://www.finngen.fi/en/access_results) after registration and approval of the summary-statistics access request:

| Endpoint | File name (place in `results/finngen_full/`) |
|---|---|
| Skin infection | `finngen_R9_L12_INFECT_SKIN.gz` |
| Cellulitis | `finngen_R9_L12_CELLULITIS.gz` |
| Lower limb ulcer | `finngen_R9_L12_ULCERLOWLIMB.gz` |
| Chronic skin ulcer | `finngen_R9_L12_CHRONICULCEROFSKIN.gz` |
| Pressure ulcer | `finngen_R9_L12_DECUBITANSULCERANDPRESSURE.sumstats.txt.gz` |

Set `MR_FINNGEN_R9_DIR` to point at the directory instead of placing the files under `results/`.

## FinnGen R12 (version-robustness rerun, required for `run_mvmr_r12_sensitivity.py`)

Same five endpoints from the R12 release. Set `MR_FINNGEN_R12_DIR`, or place the files under `results/finngen_R12/`.

## Data-use terms

- OpenGWAS: terms stated at the OpenGWAS website apply.
- FinnGen: summary statistics are made available under the terms stated on the FinnGen website at the time of download, and are licensed CC BY 4.0. Registering for access and citing the FinnGen papers is required. The repository redistributes only the small derived extracts shipped under `results/` (effects on the 352 pre-specified instruments), which is consistent with the CC BY 4.0 licence with attribution; it does not redistribute any full summary-statistics rows. Script outputs contain small FinnGen-derived extracts, so check the FinnGen terms before publishing or re-sharing those output files.

## Sample overlap between exposure and outcome data

The manuscript states that exposures and outcomes come from non-overlapping samples. The verification is documented here because it is a data-provenance check rather than a computation.

| Exposure | Summary statistics | Contributing cohorts | Contains FinnGen? |
|---|---|---|---|
| Type 2 diabetes | `ebi-a-GCST006867` (Xue et al. 2018, *Nat Commun* 9:2941) | DIAGRAM (12,171 + 22,669 cases) + GERA (6,905) + UK Biobank (21,147); 62,892 cases / 596,424 controls | No |
| Body mass index | `ebi-a-GCST006368` (Hoffmann et al. 2018, *Genetics* 210:499-515) | GERA (100,418 adults) with replication in GIANT and UK Biobank; 315,347 European-ancestry individuals | No |
| C-reactive protein | `ebi-a-GCST90029070` (Said et al. 2022, *Nat Commun* 13:2198) | UK Biobank (427,367) + CHARGE consortium; 575,531 European-ancestry individuals | No |

FinnGen is built from nine Finnish biobanks and national health registries and is an independent Finnish resource; none of the three exposure GWAS includes Finnish cohorts. The T2D GWAS predates the first public FinnGen release.

Two residual points are reported in the manuscript rather than resolved:

1. A cohort-level list for the CHARGE CRP meta-analysis was not verified item by item. If a small Finnish cohort were included, the affected sample size would be far smaller than the UK Biobank component and would not change the conclusion; this is stated as a residual uncertainty.
2. The three exposures partially overlap with each other (all three draw on UK Biobank). Exposure-to-exposure overlap does not violate the two-sample MR independence assumption, which concerns the separation of instrument and outcome samples; it is reported under STROBE-MR items 4b and 10c.
