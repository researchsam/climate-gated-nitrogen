# Data dictionary

## `county_year_panel.csv` (primary analytical file; one row per county-year)

| Column | Units | Description |
|---|---|---|
| `CountyName` | — | Iowa county name (99 counties). |
| `Year` | — | Calendar year, 2001–2019 (2000 retained only to form the first lag). |
| `loading_kgha` | kg N ha⁻¹ yr⁻¹ | Annual riverine nitrogen load, county-area-normalized; sum of monthly runoff × nitrate concentration × 0.01 over sampled months, aggregated across county–HUC8 intersection pieces. **Dependent variable.** |
| `months_available` | count | Number of calendar months with a valid nitrate sample contributing to the annual load (coverage indicator; used in the Section 3.11 coverage check). |
| `mean_concentration_mgL` | mg L⁻¹ | Mean monthly nitrate concentration for the county-year. |
| `flow_mm` | mm yr⁻¹ | Annual runoff depth. Flow-weighted mean concentration = `loading_kgha` / (0.01 × `flow_mm`). |
| `area_ha` | ha | County area. |
| `NS` | kg N ha⁻¹ yr⁻¹ | Nitrogen surplus = inputs (fertilizer + manure + fixation) − harvested-grain N removal. Negative where removal exceeds inputs. |
| `precipmm` | mm yr⁻¹ | Annual precipitation (IEM IACLIMATE, county-aggregated). |
| `Value` | acres | USDA Census tile-drained cropland acreage (drainage proxy); intercensal years filled forward/backward. |
| `NS_lag1`, `NS_lag2` | kg N ha⁻¹ | First/second-order lag of `NS`. |
| `precip_lag1` | mm | First-order lag of `precipmm`. |
| `loading_lag1` | kg N ha⁻¹ yr⁻¹ | First-order lag of `loading_kgha` (antecedent load). |
| `NS_precip` | — | Interaction term `NS` × `precipmm`. |
| `precip_z` | — | County-standardized annual precipitation (z-score within county). |
| `Value_std` | — | Standardized tile acreage. |

## `county_month_load_rebuilt.csv` (county-month, with coverage)

| Column | Units | Description |
|---|---|---|
| `CountyName`, `Year`, `Month` | — | County-month key. |
| `loading_kgha` | kg N ha⁻¹ | Monthly load (runoff × concentration × 0.01) for the covered area. |
| `load_covered_area_ha` | ha | County area with both flow and concentration data that month. |
| `n_valid_pieces`, `n_total_pieces` | count | County–HUC8 intersection pieces with valid data vs total. |
| `county_area_ha_from_pieces` | ha | County area reconstructed from intersection pieces. |
| `coverage_frac` | 0–1 | Fraction of county area with valid data that month. |

## `county_huc8_audit.csv`
Per-county count of HUC8 intersections and how many have flow, concentration, or both — documents spatial data coverage.

## `statewide_load_diagnostic_corrected.csv`
Statewide annual load (Mg), mean/median county yield, and county count per year — sanity-check series.

## `annual_load_rebuilt.csv`
County-year loads before merging surplus/precip/tile covariates.

## `Tile_data.csv`
Raw USDA Census tile-drained acreage by county and census year (pre-fill).

### Unit note
1 mm of runoff over 1 ha at 1 mg L⁻¹ = 0.01 kg N (the factor in `loading_kgha`).
