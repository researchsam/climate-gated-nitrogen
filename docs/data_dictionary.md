# Data dictionary

## `county_year_model_frame.csv`

Complete-case analytical panel used in the corrected manuscript. It contains
1,750 county-years from 98 Iowa counties during 2001–2019. Iowa has 99
counties; the sample count is 98 because one county lacks a complete set of
model variables after lag construction and complete-case filtering.

| Column | Units | Description |
|---|---|---|
| `CountyName` | — | Iowa county name. |
| `Year` | year | Calendar year. |
| `loading_kgha` | kg N ha⁻¹ yr⁻¹ | Annual estimated riverine nitrogen load normalized by monitored intersection area. Primary dependent variable. |
| `loading_kgha_fullarea` | kg N ha⁻¹ yr⁻¹ | Sensitivity load normalized by full county area. |
| `loading_kgha_lag1`–`loading_kgha_lag3` | kg N ha⁻¹ yr⁻¹ | First- through third-order load lags. |
| `NS` | kg N ha⁻¹ yr⁻¹ | Agricultural nitrogen surplus: inputs minus harvested-grain nitrogen removal. |
| `NS_lag1`, `NS_lag2` | kg N ha⁻¹ yr⁻¹ | First- and second-order nitrogen-surplus lags. |
| `precipmm` | mm yr⁻¹ | Annual county precipitation. |
| `precipmm_lag1`–`precipmm_lag3` | mm yr⁻¹ | First- through third-order precipitation lags. |
| `tile_ha` | ha | Interpolated county tile-drained acreage. |
| `county_area_ha` | ha | County area used for equal-effort mass conversion. |
| `corng_ha` | ha | Corn-grain harvested area. |
| `soy_ha` | ha | Soybean harvested area. |
| `flow_mm` | mm yr⁻¹ | Annual runoff depth used in the load calculation. |
| `mean_concentration_mgL` | mg N L⁻¹ | Annual mean nitrate concentration. |
| `mean_monitored_area_ha` | ha | Mean monitored area contributing to the annual estimate. |
| `mean_monitored_area_frac` | fraction | Mean monitored-area coverage divided by county area. |
| `months_available` | count | Number of months contributing to the annual estimate. |

The analysis script constructs `NS_precip = NS × precipmm`, county-standardized
precipitation, row-crop fraction, and other modeling variables at runtime.

### Load unit conversion

One millimeter of runoff over one hectare at a concentration of one milligram
per liter equals 0.01 kilograms of nitrogen:

`loading_kgha = flow_mm × concentration_mgL × 0.01`

