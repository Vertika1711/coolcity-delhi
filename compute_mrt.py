import pandas as pd
import numpy as np
from earthkit.meteo import solar
import thermofeel
import datetime

df = pd.read_csv("coolcity_full_dataset.csv")

all_mrt = []
unique_dates = sorted(df['date'].unique())
print(f"Processing {len(unique_dates)} unique dates...")

for i, date_str in enumerate(unique_dates):
    subset = df[df['date'] == date_str]

    begin = datetime.datetime.strptime(date_str, "%Y-%m-%d")
    end = begin + datetime.timedelta(days=1)

    cossza = solar.cos_solar_zenith_angle_integrated(
        begin_date=begin,
        end_date=end,
        latitudes=subset['centroid_lat'].values,
        longitudes=subset['centroid_lon'].values
    )

    ssrd = subset['surface_solar_radiation_downwards_sum'].values / 86400
    ssr = subset['surface_net_solar_radiation_sum'].values / 86400
    strd = subset['surface_thermal_radiation_downwards_sum'].values / 86400
    strr = subset['surface_net_thermal_radiation_sum'].values / 86400
    fdir = subset['total_sky_direct_solar_radiation_at_surface'].values / 3600
    dsrp = subset['mean_surface_direct_short_wave_radiation_flux'].values

    mrt = thermofeel.calculate_mean_radiant_temperature(
        ssrd=ssrd, ssr=ssr, dsrp=dsrp, strd=strd, fdir=fdir, strr=strr, cossza=cossza
    )

    result = pd.DataFrame({
        'cell_id': subset['cell_id'].values,
        'date': date_str,
        'MRT': mrt
    })
    all_mrt.append(result)

    if (i + 1) % 100 == 0:
        print(f"  Processed {i + 1}/{len(unique_dates)} dates...")

mrt_full = pd.concat(all_mrt, ignore_index=True)
print(f"\nMRT computed for all dates. Shape: {mrt_full.shape}")
print(f"Expected: {df.shape[0]} rows")
print(f"MRT overall range: {mrt_full['MRT'].min():.2f}K to {mrt_full['MRT'].max():.2f}K")
print(f"Null count: {mrt_full['MRT'].isnull().sum()}")

# Merge MRT into the full dataset
final = df.merge(mrt_full, on=['cell_id', 'date'], how='left')
final.to_csv("coolcity_full_dataset.csv", index=False)
print(f"\nSaved coolcity_full_dataset.csv with MRT column. Final shape: {final.shape}")