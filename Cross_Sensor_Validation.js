// =====================================================================
// Pipeline 4 + Part B: masks for all years AND cross sensor validation
// Change ONLY 'year' at the top. Run 7 times for the masks.
// Part B at the bottom runs every time, note its output once.
// =====================================================================

var year = 2019;   // <-- শুধু এটি বদলান

var gopalganj = ee.FeatureCollection("FAO/GAUL/2025/level2")
    .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
    .filter(ee.Filter.eq('GAUL1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('GAUL2_NAME', 'Gopalganj'));

var aoi = gopalganj.geometry();
Map.centerObject(gopalganj, 10);
Map.setOptions('SATELLITE');

function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(1).and(scl.neq(3)).and(scl.neq(8))
      .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return image.select(['B2','B3','B4','B8','B11'])
      .divide(10000).updateMask(mask);
}

var SEASONS = ['DryWinter','PreMonsoon','Monsoon','PostMonsoon'];

function seasonRange(yr, s) {
  var y = String(yr), p = String(yr - 1);
  if (s === 'DryWinter')   return [p+'-12-01', y+'-02-28'];
  if (s === 'PreMonsoon') return [y+'-03-01', y+'-05-31'];
  if (s === 'Monsoon')    return [y+'-06-01', y+'-09-30'];
  return [y+'-10-01', y+'-11-30'];
}

var persistent = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
    .select('occurrence').gte(75).rename('persistent');

var distToPersist = persistent.selfMask().clip(aoi.buffer(50000))
    .fastDistanceTransform(512).sqrt().multiply(30).rename('dist');

// ==================== মাস্ক ও RGB ====================
var wList = [];
SEASONS.forEach(function (s) {
  var r = seasonRange(year, s);
  var c = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(aoi).filterDate(r[0], r[1])
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
      .map(maskS2clouds).median().clip(gopalganj);
  wList.push(c.normalizedDifference(['B3','B11']).gt(0)
      .or(c.normalizedDifference(['B3','B8']).gt(0)).rename('w'));
});

var seasonal = ee.ImageCollection(wList).max().rename('w');
var core = persistent.and(seasonal).rename('core');
var margin = seasonal.and(distToPersist.gt(1000))
    .and(persistent.unmask(0).eq(0)).rename('margin');

var rgb = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(aoi).filterDate(year+'-01-01', year+'-12-31')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .map(maskS2clouds).median().clip(gopalganj)
    .select(['B4','B3','B2']).multiply(3).clamp(0,1)
    .rename(['R','G','B']).toFloat();

Map.addLayer(rgb, {bands:['R','G','B']}, 'RGB ' + year, true);
Map.addLayer(core, {palette:['0000FF']}, 'Core ' + year, false);
Map.addLayer(margin, {palette:['00BFFF']}, 'Margin ' + year, false);
Map.addLayer(persistent, {palette:['000080']}, 'Persistent water', false);

print(year + ' | core px = ' + core.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 30,
  maxPixels: 1e9, tileScale: 8}).get('core').getInfo() +
  ' | margin px = ' + margin.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 30,
  maxPixels: 1e9, tileScale: 8}).get('margin').getInfo());

Export.image.toDrive({
  image: core.toFloat(), description: 'Mask_Core_' + year,
  folder: 'Gopalganj_Turbidity_Data',
  region: aoi, scale: 30, maxPixels: 1e13
});
Export.image.toDrive({
  image: margin.toFloat(), description: 'Mask_Margin_' + year,
  folder: 'Gopalganj_Turbidity_Data',
  region: aoi, scale: 30, maxPixels: 1e13
});
Export.image.toDrive({
  image: seasonal.toFloat(), description: 'Mask_SeasonalWater_' + year,
  folder: 'Gopalganj_Turbidity_Data',
  region: aoi, scale: 30, maxPixels: 1e13
});
Export.image.toDrive({
  image: rgb, description: 'RGB_' + year,
  folder: 'Gopalganj_Turbidity_Data',
  region: aoi, scale: 30, maxPixels: 1e13
});

// বছর বদলালেও একই থাকে, একবারই Run করতে হবে
Export.image.toDrive({
  image: distToPersist.toFloat(), description: 'Raster_DistToPermanentWater',
  folder: 'Gopalganj_Turbidity_Data',
  region: aoi, scale: 30, maxPixels: 1e13
});
Export.table.toDrive({
  collection: gopalganj, description: 'Boundary_Gopalganj',
  folder: 'Gopalganj_Turbidity_Data', fileFormat: 'GeoJSON'
});

// ==================== PART B: CROSS SENSOR ====================
print('');
print('===== PART B: Sentinel 2 vs Landsat 8/9 =====');

function lsMask(im) {
  var qa = im.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(2).eq(0)
      .and(qa.bitwiseAnd(4).eq(0))
      .and(qa.bitwiseAnd(8).eq(0))
      .and(qa.bitwiseAnd(16).eq(0))
      .and(qa.bitwiseAnd(32).eq(0))
      .and(qa.bitwiseAnd(128).eq(0));
  return im.select(['SR_B3','SR_B4'])
      .multiply(0.0000275).add(0.2).updateMask(clear);
}

var LSraw = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    .merge(ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"))
    .filterBounds(aoi)
    .filterDate('2019-01-01','2026-01-01')
    .filter(ee.Filter.lt('CLOUD_COVER', 20));

print('Landsat 8+9 scenes 2019-2025:', LSraw.size().getInfo());

function s2Annual(yr) {
  var nd = [], wt = [];
  for (var i = 0; i < 4; i++) {
    var r = seasonRange(yr, SEASONS[i]);
    var c = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi).filterDate(r[0], r[1])
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
        .map(maskS2clouds).median().clip(gopalganj);
    nd.push(c.normalizedDifference(['B4','B3']).rename('ndti'));
    wt.push(c.normalizedDifference(['B3','B11']).gt(0)
        .or(c.normalizedDifference(['B3','B8']).gt(0)).rename('w'));
  }
  return {
    ndti: ee.ImageCollection(nd).mean().rename('NDTI'),
    water: ee.ImageCollection(wt).max().rename('water')
  };
}

// কোন কী তে cross-correlation আছে তা খুঁজে বের করা
function findCrossKey(d) {
  if (!d) { return null; }
  var ks = Object.keys(d);
  for (var i = 0; i < ks.length; i++) {
    var p = ks[i].split('_');
    if (p.length === 2 && p[0] !== p[1]) { return ks[i]; }
  }
  return null;
}

for (var yr = 2019; yr <= 2025; yr++) {
  try {
    var s2 = s2Annual(yr);
    var wm = s2.water;

    var lsYear = LSraw.filter(ee.Filter.calendarRange(yr, yr, 'year'));
    var lsN = lsYear.size().getInfo();
    if (lsN === 0) { print('YEAR ' + yr + ': no Landsat, skipped'); continue; }

    var s2v = s2.ndti.updateMask(wm).rename('s');
    var lsv = lsYear.map(lsMask)
        .map(function(im) {
          return im.normalizedDifference(['SR_B4','SR_B3']).rename('v');
        })
        .mean().rename('v').updateMask(wm);

    var two = s2v.addBands(lsv);

    var meanS = two.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).get('s').getInfo();
    var meanL = two.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).get('v').getInfo();
    var sdS = two.reduceRegion({
      reducer: ee.Reducer.stdDev(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).get('s').getInfo();
    var sdL = two.reduceRegion({
      reducer: ee.Reducer.stdDev(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).get('v').getInfo();
    var cnt = two.reduceRegion({
      reducer: ee.Reducer.count(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).get('s').getInfo();

    var corrRaw = two.reduceRegion({
      reducer: ee.Reducer.pearsonsCorrelation(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).getInfo();

    var ck = findCrossKey(corrRaw);
    var r = (ck !== null) ? corrRaw[ck] : null;

    var fit = two.reduceRegion({
      reducer: ee.Reducer.linearFit(), geometry: aoi, scale: 60,
      maxPixels: 1e9, tileScale: 8}).getInfo();

    print('YEAR ' + yr + ' | lsScenes=' + lsN +
          ' | meanS2=' + meanS + ' | meanLS=' + meanL +
          ' | sdS2=' + sdS + ' | sdLS=' + sdL +
          ' | n=' + cnt + ' | r=' + r +
          ' | scale=' + (fit ? fit.scale : null) +
          ' | offset=' + (fit ? fit.offset : null));

    if (yr === 2019) {
      print('   CORR KEYS: ' + JSON.stringify(corrRaw));
      print('   FIT KEYS : ' + JSON.stringify(fit));
    }

  } catch (err) {
    print('YEAR ' + yr + ' FAILED: ' + err);
  }
}

print('=== ALL DONE ===');