# Climate-Gated Nitrogen Export — code and data

Reproducibility archive for:

> Soetan, S. M., & Kaleita, A. *Climate-Gated Nitrogen Export: Hydroclimatic
> Timing Improves the Efficiency of Nitrogen-Surplus Reductions in Iowa.*

This repository regenerates every figure, table, and statistic in the paper
from publicly available inputs. The study estimates how precipitation and
subsurface tile drainage gate the conversion of agricultural nitrogen surplus
into riverine nitrogen load across 99 Iowa counties (2001–2019), and tests an
equal-effort counterfactual that reallocates a fixed surplus reduction toward
wet years.

## Repository layout

```
code/
  01_data_construction.ipynb   Build the county-year panel from USGS flow,
                               water-quality nitrate, precipitation, tile
                               drainage, and nitrogen-surplus inputs.
  02_analysis.ipynb            Random Forest + two-way fixed-effects model,
                               equal-effort timing counterfactual, robustness
                               (placebo, bootstrap, hold-out, leverage, collinearity).
  03_confirmatory_analyses.py  Load-coverage check and flow-weighted-concentration
                               test (manuscript Section 3.11).
data/
  county_year_panel.csv        Analytical panel (one row per county-year). PRIMARY FILE.
  county_month_load_rebuilt.csv County-month loads with piece-level coverage.
  annual_load_rebuilt.csv      County-year loads prior to merging covariates.
  county_huc8_audit.csv        County x HUC8 intersection coverage audit.
  statewide_load_diagnostic_corrected.csv  Statewide annual load diagnostic.
  Tile_data.csv                USDA Census tile-drained acreage (raw).
docs/
  data_dictionary.md           Column definitions and units.
  CITATION.cff                 How to cite.
LICENSE-CODE.txt               MIT (all code).
LICENSE-DATA.txt               CC-BY-4.0 (derived data products).
requirements.txt               Python dependencies.
.zenodo.json                   Deposit metadata.
```

## Quick start

```bash
python -m venv venv && source venv/bin/activate     # optional
pip install -r requirements.txt

# Reproduce the headline analysis and figures:
jupyter nbconvert --to notebook --execute code/02_analysis.ipynb

# Reproduce the Section 3.11 confirmatory tests (uses data/county_year_panel.csv):
cd code && python 03_confirmatory_analyses.py
```

`02_analysis.ipynb` reads `data/county_year_panel.csv`. Update the `DATA_DIR`
path near the top of the notebook to point at the `data/` folder in this repo.

## Headline results (for verification)

| Quantity | Value |
|---|---|
| RF temporal hold-out NSE / KGE / RMSE | 0.50 / 0.66 / 9.5 |
| RF grouped-CV NSE / KGE / RMSE | 0.64 / 0.71 / 7.8 |
| FE NS×precipitation interaction (β) | 1.97 × 10⁻⁴ (p ≈ 6×10⁻¹²) |
| Marginal-effect zero-crossing | 814 mm yr⁻¹ |
| Within-R² | 0.45 |
| Equal-effort timing dividend | ≈ 1,701 Mg N yr⁻¹ (95% CI 1,274–2,210) |
| Placebo (random-year) dividend | ≈ 2 Mg N yr⁻¹ |
| Coverage robustness (β range, ≥8–≥12 mo) | 1.83–1.95 × 10⁻⁴ |
| Flow-weighted-concentration interaction (β) | 1.85 × 10⁻⁵ (p ≈ 1.6×10⁻⁶) |

## Data provenance

All raw inputs are public:

- **Streamflow & nitrate:** USGS National Water Information System
  (https://waterdata.usgs.gov) and the Water Quality Portal
  (https://www.waterqualitydata.us). Public domain.
- **Precipitation:** Iowa Environmental Mesonet, IACLIMATE
  (https://mesonet.agron.iastate.edu).
- **Tile drainage & crop acreage:** USDA NASS QuickStats
  (https://quickstats.nass.usda.gov). Public domain.
- **Nitrogen surplus:** Iowa Food–Energy–Water (IFEWS) nitrogen-balance product
  (Raul et al. 2022; Tuthill & Kaleita 2024).

The derived products in `data/` are released under CC-BY-4.0; the code under MIT.
