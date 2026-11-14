// =====================================================================
// Pipeline 4: masks for all years, one variable controls everything
// =====================================================================

var year = 2025;   // <-- এখানে বছর বদলান

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

var persistent = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
    .select('occurrence').gte(75).rename('persistent');

var seasonKeys = ['DryWinter','PreMonsoon','Monsoon','PostMonsoon'];

function seasonRange(yr, s) {
  var y = String(yr), p = String(yr - 1);
  if (s === 'DryWinter')   return [p+'-12-01', y+'-02-28'];
  if (s === 'PreMonsoon') return [y+'-03-01', y+'-05-31'];
  if (s === 'Monsoon')    return [y+'-06-01', y+'-09-30'];
  return [y+'-10-01', y+'-11-30'];
}

// বছরের চার ঋতুর ঋতুগত পানি
var wList = [];
seasonKeys.forEach(function (s) {
  var r = seasonRange(year, s);
  var c = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(aoi).filterDate(r[0], r[1])
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
      .map(maskS2clouds).median().clip(gopalganj);
  wList.push(c.normalizedDifference(['B3','B11']).gt(0)
      .or(c.normalizedDifference(['B3','B8']).gt(0)).rename('w'));
});

var seasonal = ee.ImageCollection(wList).max().rename('w');

var distToPersist = persistent.selfMask().clip(aoi.buffer(50000))
    .fastDistanceTransform(512).sqrt().multiply(30).rename('dist');

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
Map.addLayer(distToPersist, {min:0,max:8000,palette:['00ff00','ffff00','ff0000']},
    'Dist to permanent water', false);

// ---- যাচাই ----
var cc = core.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 30,
  maxPixels: 1e9, tileScale: 8}).get('core').getInfo();
var mc = margin.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 30,
  maxPixels: 1e9, tileScale: 8}).get('margin').getInfo();
print(year + ' | core pixels = ' + cc + ' | margin pixels = ' + mc);

// ---- এক্সপোর্ট, নামে বছর সহ ----
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

print('4 export tasks created for ' + year + '. Run them from the Tasks tab.');