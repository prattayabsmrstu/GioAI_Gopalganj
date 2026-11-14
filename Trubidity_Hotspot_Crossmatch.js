// ============================================================
// VERIFY ব্লক: সম্পূর্ণ স্বয়ংসম্পূর্ণ
// ============================================================

var VAoi = ee.FeatureCollection("FAO/GAUL/2025/level2")
    .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
    .filter(ee.Filter.eq('GAUL1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('GAUL2_NAME', 'Gopalganj'))
    .geometry();

var VPersistent = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
    .select('occurrence').gte(75).rename('persistent');

function VMask(image) {
  var scl = image.select('SCL');
  var ok = scl.neq(1).and(scl.neq(3)).and(scl.neq(8))
      .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return image.select(['B2', 'B3', 'B4', 'B8', 'B11'])
      .divide(10000).updateMask(ok);
}

// ---- ধাপ ১: ডায়াগনস্টিক ----
var Vraw = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(VAoi)
    .filterDate('2024-06-01', '2024-09-30')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60));

print('VERIFY | raw scene count =', Vraw.size());
print('VERIFY | raw band names =', Vraw.first().bandNames());

var Vmasked = Vraw.map(VMask);
print('VERIFY | masked band names =', Vmasked.first().bandNames());

var Vcomp = Vmasked.median().clip(VAoi);
print('VERIFY | composite band names =', Vcomp.bandNames());

// ---- ধাপ ২: জেলার আয়তন ও স্থায়ী জলভাগ ----
print('VERIFY | district area km2 =', VAoi.area(1).divide(1e6));

var VpCount = VPersistent.reduceRegion({
  reducer: ee.Reducer.count(), geometry: VAoi, scale: 30,
  maxPixels: 1e9, tileScale: 8
}).get('persistent');
print('VERIFY | persistent water pixels at 30m =', VpCount);
print('VERIFY | persistent water area km2 =',
      ee.Number(VpCount).multiply(900).divide(1e6));

// ---- ধাপ ৩: NDTI তিন ভাবে ----
var Vndti = Vcomp.normalizedDifference(['B4', 'B3']).rename('NDTI');
var Vmndwi = Vcomp.normalizedDifference(['B3', 'B11']);
var Vndwi = Vcomp.normalizedDifference(['B3', 'B8']);
var Vwater = Vmndwi.gt(0).or(Vndwi.gt(0)).rename('water');

function VStats(img, dom, label) {
  var r = img.updateMask(dom).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: VAoi, scale: 30,
    maxPixels: 1e9, tileScale: 8
  }).get('NDTI');
  var c = img.updateMask(dom).reduceRegion({
    reducer: ee.Reducer.count(), geometry: VAoi, scale: 30,
    maxPixels: 1e9, tileScale: 8
  }).get('NDTI');
  print(label + ' | mean NDTI =', r, '| pixel count =', c);
}

VStats(Vndti, VPersistent, 'A persistent only');
VStats(Vndti, VPersistent.and(Vwater), 'B persistent + seasonal');
VStats(Vndti, Vwater, 'C all seasonal water');
VStats(Vndti, Vwater.and(VPersistent.not()), 'D seasonal-only, outside persistent');

// ---- ধাপ ৪: ২৮-ব্যান্ড স্ট্যাকের গণনা যাচাই ----
var Vlist = [];
['DryWinter', 'PreMonsoon', 'Monsoon', 'PostMonsoon'].forEach(function (s) {
  var start, end;
  if (s === 'DryWinter')       { start = '2023-12-01'; end = '2024-02-28'; }
  else if (s === 'PreMonsoon') { start = '2024-03-01'; end = '2024-05-31'; }
  else if (s === 'Monsoon')    { start = '2024-06-01'; end = '2024-09-30'; }
  else                          { start = '2024-10-01'; end = '2024-11-30'; }

  var c = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(VAoi).filterDate(start, end)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
      .map(VMask).median().clip(VAoi);

  var n = c.normalizedDifference(['B4', 'B3']).rename('NDTI');
  var m = c.normalizedDifference(['B3', 'B11']).gt(0)
      .or(c.normalizedDifference(['B3', 'B8']).gt(0));

  Vlist.push(n.updateMask(VPersistent).updateMask(m).rename(s));
});

var Vstack = ee.Image.cat(Vlist);
print('VERIFY | stacked band count =', Vstack.bandNames().size());

var VstackCount = Vstack.mask().reduceRegion({
  reducer: ee.Reducer.count(), geometry: VAoi, scale: 30,
  maxPixels: 1e9, tileScale: 8
});
print('VERIFY | stacked mask count (all bands) =',
      VstackCount.get('DryWinter'));

var VstackMeans = Vstack.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: VAoi, scale: 30,
  maxPixels: 1e9, tileScale: 8
});
print('VERIFY | stacked means =', VstackMeans);

// ---- ধাপ ৫: ঋতুভিত্তিক গড়, stack থেকে ----
print('VERIFY | seasonal means from 28 band stack:');
['DryWinter', 'PreMonsoon', 'Monsoon', 'PostMonsoon'].forEach(function (s) {
  var m = Vstack.select(s).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: VAoi, scale: 30,
    maxPixels: 1e9, tileScale: 8
  }).get(s);
  var c = Vstack.select(s).reduceRegion({
    reducer: ee.Reducer.count(), geometry: VAoi, scale: 30,
    maxPixels: 1e9, tileScale: 8
  }).get(s);
  print('  ' + s + ' : mean =', m, '| count =', c);
});