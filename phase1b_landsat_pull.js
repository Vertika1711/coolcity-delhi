// ============================================================
// PHASE 1B — Satellite Land-Surface Covariates
// NDVI, NDBI, Albedo, LST — Delhi 2km Grid, April-June 2019-2026
// CoolCity Project

// REVISED VERSION: creates one export task PER YEAR (8 tasks total)
// instead of one giant 8-year task, to avoid the
// "Too many concurrent aggregations" error.
// ============================================================


// ------------------------------------------------------------
// STEP 1: Reuse exact grid-creation code, unchanged.
// This guarantees cell_id values match exactly so the
// two CSVs can be joined later on cell_id + date.
// ------------------------------------------------------------

var admin2 = ee.FeatureCollection("FAO/GAUL/2015/level2");

// Filter to Delhi
var delhi = admin2.filter(ee.Filter.eq('ADM1_NAME', 'Delhi'));

// Get the bounding box of Delhi
var delhiBounds = delhi.geometry().bounds();

// Create the same 2 km grid as Phase 1a
var proj = ee.Projection('EPSG:32643').atScale(2000);

// Generate grid cells covering Delhi's bounding box
var grid = delhiBounds.coveringGrid(proj);

// Keep only cells that intersect Delhi
var delhiGrid = grid.filterBounds(delhi);

// Assign the same cell_id as phase1a script
var delhiGridWithID = delhiGrid.map(function(feature) {
  return feature.set('cell_id', feature.get('system:index'));
});

// Check number of grid cells
print('Number of grid cells:', delhiGridWithID.size());

// Visualize
Map.centerObject(delhi, 9);
Map.addLayer(delhiGridWithID, {color: 'blue'}, 'Delhi 2km Grid');


// ------------------------------------------------------------
// STEP 2 (REVISED): Cloud masking for Landsat Collection 2, Level 2.
//
// QA_PIXEL is a quality-assurance band, one bit per condition.
// Our original version only checked bit 3 (cloud) and bit 4
// (cloud shadow). That let thin/edge cloud pixels through —
// they show up later as impossibly cold LST values (thermal
// sensors read cloud tops as very cold), which is exactly what
// we saw in the exported CSVs.
//
// bit 1 = dilated cloud  (a buffer zone around detected clouds —
//         catches cloud edges the core cloud bit misses)
// bit 2 = cirrus          (thin high-altitude cloud)
// bit 3 = cloud
// bit 4 = cloud shadow
//
// We now keep a pixel only if ALL FOUR bits are 0.
// ------------------------------------------------------------

function maskLandsatClouds(image) {

  var qa = image.select('QA_PIXEL');

  var dilatedCloudBit = 1 << 1;
  var cirrusBit = 1 << 2;
  var cloudBit = 1 << 3;
  var cloudShadowBit = 1 << 4;

  var mask = qa.bitwiseAnd(dilatedCloudBit).eq(0)
      .and(qa.bitwiseAnd(cirrusBit).eq(0))
      .and(qa.bitwiseAnd(cloudBit).eq(0))
      .and(qa.bitwiseAnd(cloudShadowBit).eq(0));

  return image.updateMask(mask);
}


// ------------------------------------------------------------
// STEP 3: Convert Landsat Collection 2 scaled integers into
// physical units.
//
// Optical bands:
// reflectance = DN * 0.0000275 - 0.2
//
// Thermal band:
// temperature (Kelvin) = DN * 0.00341802 + 149.0
// ------------------------------------------------------------

function applyScaleFactors(image) {

  var opticalBands = image.select('SR_B.')
    .multiply(0.0000275)
    .add(-0.2);

  var thermalBand = image.select('ST_B10')
    .multiply(0.00341802)
    .add(149.0);

  return image
    .addBands(opticalBands, null, true)
    .addBands(thermalBand, null, true);
}


// ------------------------------------------------------------
// STEP 4: Calculate NDVI, NDBI, Albedo and LST.
//
// Landsat 8/9 bands:
// B2 = Blue
// B4 = Red
// B5 = NIR
// B6 = SWIR1
// B7 = SWIR2
// B10 = Thermal
// ------------------------------------------------------------

function addIndices(image) {

  // NDVI = (NIR - Red) / (NIR + Red)
  var ndvi = image
    .normalizedDifference(['SR_B5', 'SR_B4'])
    .rename('NDVI');


  // NDBI = (SWIR1 - NIR) / (SWIR1 + NIR)
  var ndbi = image
    .normalizedDifference(['SR_B6', 'SR_B5'])
    .rename('NDBI');


  // Broadband albedo approximation
  // Liang (2001) empirical formula for Landsat OLI
  var albedo = image.expression(
    '0.356*blue + 0.130*red + 0.373*nir + ' +
    '0.085*swir1 + 0.072*swir2 - 0.0018', {

      'blue': image.select('SR_B2'),
      'red': image.select('SR_B4'),
      'nir': image.select('SR_B5'),
      'swir1': image.select('SR_B6'),
      'swir2': image.select('SR_B7')

    }).rename('ALBEDO');


  // LST is already in Kelvin after scaling
  var lst = image
    .select('ST_B10')
    .rename('LST');


  return image.addBands([
    ndvi,
    ndbi,
    albedo,
    lst
  ]);
}


// ------------------------------------------------------------
// STEP 5: Load Landsat 8 and Landsat 9.
//
// We:
// 1. Filter to Delhi
// 2. Remove clouds/cloud shadows
// 3. Convert values to physical units
// 4. Calculate NDVI, NDBI, Albedo and LST
// ------------------------------------------------------------

var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(delhiGridWithID)
  .map(maskLandsatClouds)
  .map(applyScaleFactors)
  .map(addIndices);


var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(delhiGridWithID)
  .map(maskLandsatClouds)
  .map(applyScaleFactors)
  .map(addIndices);


// Merge Landsat 8 and Landsat 9
var landsatMerged = l8
  .merge(l9)
  .select([
    'NDVI',
    'NDBI',
    'ALBEDO',
    'LST'
  ]);


// ------------------------------------------------------------
// STEP 6: Create a composite for one time window.
//
// Satellite images aren't available every day.
// Therefore, for each 16-day window:
//
// 1. Take all available Landsat images
// 2. Take the median value for each pixel
// 3. Average those pixels inside each 2 km grid cell
//
// Result = one row per grid cell.
// ------------------------------------------------------------

function getCompositeForWindow(windowStart, windowEnd) {

  var windowImgs = landsatMerged
    .filterDate(windowStart, windowEnd);

  var composite = windowImgs.median();

  var reduced = composite.reduceRegions({

    collection: delhiGridWithID,

    reducer: ee.Reducer.mean(),

    scale: 30

  });

  return reduced;
}


// ------------------------------------------------------------
// STEP 7: Give every calendar date inside a 16-day window
// the values from that window's satellite composite.
//
// Example:
//
// April 1 → composite value
// April 2 → same composite value
// ...
// April 16 → same composite value
//
// Then a new satellite window begins.
// ------------------------------------------------------------

function expandToDates(compositeFC, dateStrList) {

  return ee.FeatureCollection(

    dateStrList.map(function(dateStr) {

      return compositeFC.map(function(f) {

        return f.set('date', dateStr);

      });

    })

  ).flatten();
}


// ------------------------------------------------------------
// STEP 8 (REVISED): Build one year's worth of rows.
//
// This is exactly the same window/date logic as before,
// just wrapped in a function so we can call it once per year
// and export each year separately, instead of building one
// giant 8-year FeatureCollection.
// ------------------------------------------------------------

var WINDOW_DAYS = 16;

function buildYearFC(year) {

  year = ee.Number(year);

  // April 1
  var seasonStart = ee.Date.fromYMD(year, 4, 1);

  // July 1 (exclusive), so June 30 is included
  var seasonEnd = ee.Date.fromYMD(year, 6, 30).advance(1, 'day');

  var totalDays = seasonEnd.difference(seasonStart, 'day');

  // Generate 16-day window offsets
  var windowOffsets = ee.List.sequence(0, totalDays, WINDOW_DAYS);

  var yearResults = windowOffsets.map(function(offset) {

    var windowStart = seasonStart.advance(offset, 'day');

    // Candidate end of the 16-day window
    var candidateEndMillis = windowStart.advance(WINDOW_DAYS, 'day').millis();

    // Do not go beyond June 30
    var windowEndMillis = ee.Number(candidateEndMillis).min(seasonEnd.millis());

    var windowEnd = ee.Date(windowEndMillis);

    // Number of calendar days in this window
    var numDaysInWindow = windowEnd.difference(windowStart, 'day');

    // Create dates for this window
    var dateStrList = ee.List.sequence(0, numDaysInWindow.subtract(1))
      .map(function(d) {
        return windowStart.advance(d, 'day').format('YYYY-MM-dd');
      });

    // Create satellite composite
    var composite = getCompositeForWindow(windowStart, windowEnd);

    // Give the composite value to every date inside the window
    return expandToDates(composite, dateStrList);

  });

  return ee.FeatureCollection(yearResults).flatten();
}


// ------------------------------------------------------------
// STEP 9 (REVISED): Create ONE export task per year.
//
// Instead of one 8-year export, this registers 8 separate
// tasks in the Tasks tab — one per year. Each one is a much
// smaller computation, so it stays under GEE's aggregation
// limits. You'll run each task individually.
// ------------------------------------------------------------

var allYears = [2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026];

allYears.forEach(function(year) {

  var yearFC = buildYearFC(year);

  Export.table.toDrive({
    collection: yearFC,
    description: 'CoolCity_Landsat_Phase1b_' + year,
    folder: 'CoolCity_Exports',
    fileNamePrefix: 'coolcity_landsat_phase1b_' + year,
    fileFormat: 'CSV',
    selectors: ['cell_id', 'date', 'NDVI', 'NDBI', 'ALBEDO', 'LST']
  });

});

print('8 export tasks created — go to the Tasks tab and click Run on each one.');


// ------------------------------------------------------------
// OPTIONAL EXTENSION — SENTINEL-2
//
// Sentinel-2 has 10m resolution and could provide higher
// resolution NDVI/NDBI/Albedo.
//
// However, Sentinel-2 does NOT have a thermal band,
// so it cannot provide LST.
//
// Therefore, the current Phase 1b implementation uses
// Landsat 8 + Landsat 9 consistently for all four variables.
// ------------------------------------------------------------
