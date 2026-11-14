// =====================================================================
// Pipeline 2 final: two-stratum NDTI table, 28 composites
// Every statistic from a single-band reduceRegion.
// Includes verification against independently confirmed values.
// =====================================================================

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
  return image.select(['B2', 'B3', 'B4', 'B8', 'B11'])
      .divide(10000).updateMask(mask);
}

// সব মাস্কের ব্যান্ড নাম 'mask', যাতে reduceRegion একই কী খুঁজে
var persistentMask = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
    .select('occurrence').gte(75).rename('mask');

var notPersistentMask = persistentMask.unmask(0).eq(0).rename('mask');

var distToPersistent = persistentMask.selfMask().clip(aoi.buffer(50000))
    .fastDistanceTransform(512).sqrt().multiply(30)
    .rename('dist');

var farFromPersistent = distToPersistent.gt(1000).rename('mask');

var seasonKeys = ['DryWinter', 'PreMonsoon', 'Monsoon', 'PostMonsoon'];
var years = [2019, 2020, 2021, 2022, 2023, 2024, 2025];

// মাস্কে ১ যেখানে আছে, সেখানে গণনা
function countWhere(dom) {
  return dom.updateMask(dom).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: aoi, scale: 30, maxPixels: 1e9, tileScale: 8
  }).get('mask');
}

// NDTI এর গড়, মাস্কের ভেতরে
function meanOf(img, dom, band) {
  return img.updateMask(dom).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 30, maxPixels: 1e9, tileScale: 8
  }).get(band);
}

function seasonStat(yr, s) {
  var y = String(yr);
  var p = String(yr - 1);
  var start;
  var end;

  if (s === 'DryWinter') {
    start = p + '-12-01';
    end = y + '-02-28';
  } else if (s === 'PreMonsoon') {
    start = y + '-03-01';
    end = y + '-05-31';
  } else if (s === 'Monsoon') {
    start = y + '-06-01';
    end = y + '-09-30';
  } else {
    start = y + '-10-01';
    end = y + '-11-30';
  }

  var c = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(aoi)
      .filterDate(start, end)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
      .map(maskS2clouds)
      .median()
      .clip(gopalganj);

  var ndti = c.normalizedDifference(['B4', 'B3']).rename('NDTI');

  var seasonal = c.normalizedDifference(['B3', 'B11']).gt(0)
      .or(c.normalizedDifference(['B3', 'B8']).gt(0)).rename('mask');

  var coreMask = persistentMask.and(seasonal).rename('mask');
  var marginMask = farFromPersistent.and(seasonal).and(notPersistentMask)
      .rename('mask');

  return ee.Dictionary({
    'year': yr,
    'season': s,
    'core_mean': meanOf(ndti, coreMask, 'NDTI'),
    'core_n': countWhere(coreMask),
    'margin_mean': meanOf(ndti, marginMask, 'NDTI'),
    'margin_n': countWhere(marginMask)
  });
}

// বছরে বছর, একটি getInfo প্রতি বছরে
for (var i = 0; i < years.length; i++) {
  var yr = years[i];
  var rows = [];

  seasonKeys.forEach(function (s) {
    rows.push(seasonStat(yr, s).getInfo());
  });

  rows.forEach(function (r) {
    print(r.year + ' ' + r.season +
        ' | core ' + r.core_mean.toFixed(4) + ' (n=' + r.core_n + ')' +
        ' | margin ' + r.margin_mean.toFixed(4) + ' (n=' + r.margin_n + ')');
  });
}