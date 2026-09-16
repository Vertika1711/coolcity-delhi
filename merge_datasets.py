import pandas as pd

# ============================================================
# CoolCity — Dataset Merge Script
# Combines: ERA5-Land weather + Landsat satellite covariates +
# supplementary ERA5 variables (pressure, radiation components) +
# corrected direct-radiation values (resampled to fix null gaps)
# Produces: coolcity_full_dataset.csv, coolcity_grid_reference.csv
# ============================================================

# ---- Step 1: Load base ERA5-Land weather dataset ----
era5 = pd.read_csv("raw_data/coolcity_era5_full.csv")
print(f"ERA5 base loaded: {era5.shape}")

# Extract grid geometry once into a separate reference file
# (identical across all ~728 date-rows per cell -- de-duplicated,
# not discarded) before dropping it from the main dataset.
if '.geo' in era5.columns:
    grid_reference = era5[['cell_id', '.geo']].drop_duplicates(subset='cell_id').reset_index(drop=True)
    grid_reference.to_csv("coolcity_grid_reference.csv", index=False)
    print(f"Saved coolcity_grid_reference.csv ({grid_reference.shape[0]} unique cells)")
    era5 = era5.drop(columns=['.geo'])

# system:index is GEE's internal per-export row numbering, not a
# stable identifier -- cell_id is the real one, safe to drop.
if 'system:index' in era5.columns:
    era5 = era5.drop(columns=['system:index'])

print(f"ERA5 cleaned shape: {era5.shape}\n")


# ---- Step 2: Load and combine the 8 yearly satellite files ----
satellite_dfs = [pd.read_csv(f"raw_data/coolcity_landsat_phase1b_{y}.csv") for y in range(2019, 2027)]
satellite = pd.concat(satellite_dfs, ignore_index=True)
print(f"Satellite combined shape: {satellite.shape}")
print(f"Satellite null counts:\n{satellite.isnull().sum().to_string()}\n")


# ---- Step 3: Load and combine the 8 yearly supplementary ERA5 files ----
# (surface pressure + radiation components needed for MRT/Brimicombe)
supplementary_dfs = [pd.read_csv(f"raw_data/coolcity_supplementary_{y}.csv") for y in range(2019, 2027)]
supplementary = pd.concat(supplementary_dfs, ignore_index=True)

# surface_solar_radiation_downwards_sum already exists in the base
# ERA5 file -- drop the duplicate here rather than after merging.
supplementary = supplementary.drop(columns=['surface_solar_radiation_downwards_sum'])

# The direct-radiation columns in this file were pulled with a
# daylight-only mask that caused a structural null gap (~43% of
# rows, ~190/445 cells failing at full ERA5's coarse ~28km
# resolution). Corrected versions are loaded separately in Step 4
# -- drop the broken columns here so they're never merged in.
supplementary = supplementary.drop(columns=[
    'total_sky_direct_solar_radiation_at_surface',
    'mean_surface_direct_short_wave_radiation_flux'
])

print(f"Supplementary combined shape (deduped): {supplementary.shape}\n")


# ---- Step 4: Load and combine the 8 yearly CORRECTED direct-radiation files ----
# (resampled to EPSG:32643 @ 2000m via bilinear interpolation before
# extraction, which eliminated the coarse-resolution null gap entirely)
directrad_dfs = [pd.read_csv(f"raw_data/coolcity_directrad_resampled_{y}.csv") for y in range(2019, 2027)]
directrad = pd.concat(directrad_dfs, ignore_index=True)
print(f"Direct radiation (resampled) combined shape: {directrad.shape}")
print(f"Direct radiation null counts:\n{directrad.isnull().sum().to_string()}\n")


# ---- Step 5: Merge everything together on cell_id + date ----
final_dataset = era5.merge(satellite, on=['cell_id', 'date'], how='left')
final_dataset = final_dataset.merge(supplementary, on=['cell_id', 'date'], how='left')
final_dataset = final_dataset.merge(directrad, on=['cell_id', 'date'], how='left')

print(f"Final merged shape: {final_dataset.shape}")
print(f"Final columns: {list(final_dataset.columns)}\n")
print(f"Final null counts:\n{final_dataset.isnull().sum().to_string()}\n")


# ---- Step 6: Save ----
final_dataset.to_csv("coolcity_full_dataset.csv", index=False)
print("Saved coolcity_full_dataset.csv")
print("Saved coolcity_grid_reference.csv")

# ---- Step 7: Final sanity check ----
final_check = pd.read_csv("coolcity_full_dataset.csv")
print(f"\nFinal dataset: {final_check.shape[0]} rows, {final_check.shape[1]} columns")
print(f"Expected: 323960 rows, 20 columns")
print(f"Match: {final_check.shape == (323960, 20)}")