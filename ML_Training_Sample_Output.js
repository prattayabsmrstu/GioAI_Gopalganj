// ==================== ১. AOI ====================
var gopalganj = ee.FeatureCollection("FAO/GAUL/2025/level2")
    .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
    .filter(ee.Filter.eq('GAUL1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('GAUL2_NAME', 'Gopalganj'));

var aoi = gopalganj.geometry();
Map.centerObject(gopalganj, 10);
Map.setOptions('SATELLITE');

// ==================== ২. ফাংশনস ====================
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(1).and(scl.neq(3)).and(scl.neq(8))
      .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return image.select(['B2', 'B3', 'B4', 'B8', 'B11'])
      .divide(10000).updateMask(mask);
}

// ==================== ৩. ড্রাইভার ডাটা ====================
var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD");
var modisLST = ee.ImageCollection("MODIS/061/MOD11A2");
var srtm = ee.Image("USGS/SRTMGL1_003").select('elevation').rename('Elevation');

var distToWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
    .select('occurrence').gte(75).selfMask()
    .clip(aoi.buffer(50000))
    .fastDistanceTransform(256).sqrt().multiply(30).rename('DistToWater');

// ==================== ৪. লুপ ও এক্সপোর্ট পাইপলাইন ====================
var seasonKeys = ['DryWinter', 'PreMonsoon', 'Monsoon', 'PostMonsoon'];
var sampleList = [];

for (var yr = 2019; yr <= 2025; yr++) {
  seasonKeys.forEach(function(s) {
    var start = (s === 'DryWinter') ? (yr-1) + '-12-01' : (s === 'PreMonsoon' ? yr + '-03-01' : (s === 'Monsoon' ? yr + '-06-01' : yr + '-10-01'));
    var end = (s === 'DryWinter') ? yr + '-02-28' : (s === 'PreMonsoon' ? yr + '-05-31' : (s === 'Monsoon' ? yr + '-09-30' : yr + '-11-30'));
    var d0 = ee.Date(start);

    var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
        .map(maskS2clouds);

    if (col.size().getInfo() > 0) {
      var comp = col.median().clip(aoi);
      var ndti = comp.normalizedDifference(['B4', 'B3']).rename('NDTI');
      var mndwi = comp.normalizedDifference(['B3', 'B11']).rename('MNDWI');
      var ndwi = comp.normalizedDifference(['B3', 'B8']).rename('NDWI');
      
      var rainSeason = chirps.filterDate(start, end).sum().rename('Rainfall');
      var rain30 = chirps.filterDate(d0.advance(-30, 'day'), d0).sum().rename('Rain_30d');
      var rain7  = chirps.filterDate(d0.advance(-7, 'day'), d0).sum().rename('Rain_7d');
      var lst = modisLST.filterDate(start, end).select('LST_Day_1km').mean()
          .multiply(0.02).subtract(273.15).rename('LST');

      var stack = comp.addBands([ndti, mndwi, ndwi, rainSeason, rain30, rain7, lst, srtm, distToWater])
                      .addBands(ee.Image.constant(yr).rename('Year'));
      
      var water = mndwi.gt(0).or(ndwi.gt(0));
      
      var samp = stack.updateMask(water).sample({
        region: aoi, scale: 60, numPixels: 5000, geometries: true, tileScale: 16, seed: 42
      }).map(function(f) { return f.set('Season', s); });
      
      sampleList.push(samp);
    }
  });
}

var mlSamples = ee.FeatureCollection(sampleList).flatten();
print('Total samples:', mlSamples.size());

Export.table.toDrive({
  collection: mlSamples,
  description: 'Gopalganj_ML_Samples_FINAL_CORRECTED_V2',
  folder: 'Gopalganj_Turbidity_Data',
  fileFormat: 'CSV'
});
print('Task created. Open the Tasks tab and press Run.');