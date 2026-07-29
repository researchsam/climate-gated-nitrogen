# Climate-Gated Nitrogen Export — corrected code and data

Reproducibility archive for:

> Soetan, S. M., & Kaleita, A. *Climate-Gated Nitrogen Export: An Equal-Effort
> Analysis of Nitrogen-Surplus Reduction Timing in Iowa.*

This repository supports an observational analysis of how annual precipitation
conditions the association between agricultural nitrogen surplus and estimated
riverine nitrogen load in Iowa. It combines leakage-resistant Random Forest
validation, two-way fixed-effects regression, an equal-effort timing
counterfactual, and multiple sensitivity analyses.

The analytical sample contains 1,750 county-years from 98 Iowa counties during
2001–2019. One county is absent from the complete-case model frame because it
does not have all variables required for the analysis; this is not a claim that
Iowa has only 98 counties.

## Repository layout

```text
code/
  01_data_construction.ipynb
      Builds the county-level load products from the public source data.
  02_analysis_corrected.py
      Canonical corrected analysis. Produces all four manuscript figures,
      headline statistics, robustness results, and machine-readable tables.
  03_confirmatory_analyses.py
      Additional load-coverage and flow-weighted-concentration checks.

data/
  county_year_model_frame.csv
      Complete-case analytical panel used by 02_analysis_corrected.py.
  county_month_load_rebuilt.csv
      County-month load and monitoring-coverage data.
  annual_load_rebuilt.csv
      County-year loads before merging the model covariates.
  county_huc8_audit.csv
      County-by-HUC8 monitoring-coverage audit.
  statewide_load_diagnostic_corrected.csv
      Statewide annual load diagnostics.
  Tile_data.csv
      USDA Census tile-drained acreage.

results/
  final_results.json
      Machine-readable headline and sensitivity results.
  Fig1_RF_corrected.*
  Fig2_marginal_effect_corrected.*
  Fig3_timing_corrected.*
  Fig4_robustness_corrected.*
  *.csv
      Cross-validation, counterfactual, influence, coverage, and importance
      results.

docs/
  data_dictionary.md
      Column definitions, units, and analytical-sample notes.
```

## Quick start

Create an environment and install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.\.venv\Scripts\Activate.ps1
```

Run the corrected analysis from the repository root:

```bash
python code/02_analysis_corrected.py
```

By default, the script reads `data/county_year_model_frame.csv`, uses 300
coefficient-refitting bootstrap replicates, and writes outputs to `results/`.
Paths and the bootstrap count can be changed explicitly:

```bash
python code/02_analysis_corrected.py \
  --data data/county_year_model_frame.csv \
  --output results \
  --bootstrap-replicates 300
```

## Headline results for verification

| Quantity | Corrected value |
|---|---:|
| Analytical sample | 1,750 county-years; 98 counties; 2001–2019 |
| RF temporal hold-out NSE / KGE / RMSE | 0.405 / 0.651 / 9.910 |
| RF grouped-CV mean NSE / KGE / RMSE | 0.589 / 0.688 / 8.273 |
| FE NS × precipitation coefficient | 1.853 × 10⁻⁴ |
| FE interaction p-value | 1.02 × 10⁻⁹ |
| Marginal-association zero crossing | 843.9 mm yr⁻¹ |
| Within-R² | 0.433 |
| Primary equal-effort timing dividend | 980 Mg N yr⁻¹ |
| Bootstrap median and 95% interval | 994 Mg N yr⁻¹; 406–1,602 |
| Unconstrained sensitivity dividend | 1,401 Mg N yr⁻¹ |
| 20% wet-year misclassification dividend | 557 Mg N yr⁻¹ |
| Later-period hold-out dividend | 1,026 Mg N yr⁻¹ |
| Random-year placebo mean | −14 Mg N yr⁻¹ |

## Counterfactual interpretation

Both schedules impose the same cumulative 15% reduction in positive nitrogen
surplus within each county. The primary analysis is deliberately conservative:
negative fitted marginal responses in dry years are set to zero. The 980 Mg N
yr⁻¹ estimate is therefore the manuscript’s primary timing dividend. The
1,401 Mg N yr⁻¹ estimate retains model-implied negative dry-year responses and
is reported only as a sensitivity analysis.

These are model-based equal-effort comparisons, not forecasts of realized
policy outcomes or estimates from a randomized intervention.

## Data provenance

The source observations are publicly available:

- Streamflow: USGS National Water Information System.
- Nitrate: Water Quality Portal.
- Precipitation: Iowa Environmental Mesonet IACLIMATE.
- Tile drainage and crop acreage: USDA NASS QuickStats.
- Nitrogen surplus: Iowa Food–Energy–Water nitrogen-balance product.

Derived data products are released under CC BY 4.0. Code is released under the
MIT License.

## Citation and archive

Use the version-specific Zenodo DOI associated with the GitHub release used in
your analysis. The concept DOI for all versions is:

<https://doi.org/10.5281/zenodo.20721859>

