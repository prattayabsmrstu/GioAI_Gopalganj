// ==================== ১. AOI ====================
var gopalganj = ee.FeatureCollection("FAO/GAUL/2025/level2")
    .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
    .filter(ee.Filter.eq('GAUL1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('GAUL2_NAME', 'Gopalganj'));

var aoi = gopalganj.geometry();
Map.centerObject(gopalganj, 10);
Map.setOptions('SATELLITE');

// ==================== ২. ক্লাউড মাস্কিং ====================
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(1).and(scl.neq(3)).and(scl.neq(8))
      .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return image.select(['B2', 'B3', 'B4', 'B8', 'B11'])
      .divide(10000).updateMask(mask);
}

// ==================== ৩. ২৮টি কম্পোজিট (২০১৯-২০২৫) ====================
var seasonKeys = ['DryWinter', 'PreMonsoon', 'Monsoon', 'PostMonsoon'];
var kmBands = ['B2', 'B3', 'B4', 'B8', 'B11', 'NDTI'];
var owtImages = [];

for (var yr = 2019; yr <= 2025; yr++) {
  for (var i = 0; i < seasonKeys.length; i++) {
    var s = seasonKeys[i];
    var y = String(yr);
    var p = String(yr - 1);
    var start, end;

    if (s === 'DryWinter')       { start = p + '-12-01'; end = y + '-02-28'; }
    else if (s === 'PreMonsoon') { start = y + '-03-01'; end = y + '-05-31'; }
    else if (s === 'Monsoon')    { start = y + '-06-01'; end = y + '-09-30'; }
    else                         { start = y + '-10-01'; end = y + '-11-30'; }

    var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
        .map(maskS2clouds);

    var comp = col.median().clip(gopalganj);

    var ndti   = comp.normalizedDifference(['B4', 'B3']).rename('NDTI');
    var mndwi  = comp.normalizedDifference(['B3', 'B11']).rename('MNDWI');
    var ndwi   = comp.normalizedDifference(['B3', 'B8']).rename('NDWI');
    var gr     = comp.select('B3').divide(comp.select('B4')).rename('GR_Ratio');
    var redRef = comp.select('B4').rename('RedRefl');

    owtImages.push(
      comp.addBands([ndti, mndwi, ndwi, gr, redRef])
          .toFloat()
          .set('year', yr)
          .set('season', s)
    );
  }
}

var owtCollection = ee.ImageCollection(owtImages);
print('Total composites:', owtCollection.size());

// ==================== ৪. স্থায়ী জলভাগে গড় NDTI ====================
var persistent = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
    .select('occurrence').gte(75);

var ndtiWater = owtCollection.map(function (img) {
  var water = img.select('MNDWI').gt(0).or(img.select('NDWI').gt(0));
  return img.select('NDTI').updateMask(persistent).updateMask(water);
});

var meanNdtiImg = ndtiWater.mean().clip(gopalganj);

print('Mean NDTI on persistent water, 2019-2025:',
  meanNdtiImg.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi, scale: 60, maxPixels: 1e9, tileScale: 8
  }));

// ==================== ৫. K-Means ট্রেনিং ====================
var kmTraining = owtCollection.map(function (img) {
  var water = img.select('MNDWI').gt(0)
      .or(img.select('NDWI').gt(0)).rename('water');

  return img.select(kmBands).addBands(water).toFloat().sample({
    region: aoi, scale: 60, numPixels: 2000,
    geometries: false, tileScale: 8, seed: 42
  }).filter(ee.Filter.eq('water', 1)).select(kmBands);
}).flatten();

print('K-Means training samples:', kmTraining.size());

var kmeans = ee.Clusterer.wekaKMeans(4).train(kmTraining, kmBands);

// ==================== ৬. ক্লাস্টার ব্যাখ্যা ====================
var meanSpectra = owtCollection.select(kmBands).mean().clip(gopalganj);

var meanClusters = meanSpectra.cluster(kmeans).rename('cluster');
var statImage = meanNdtiImg.addBands(meanClusters);

var clusterStats = statImage.reduceRegion({
  reducer: ee.Reducer.mean().group(1, 'cluster'),
  geometry: aoi, scale: 60, maxPixels: 1e9, tileScale: 8
});

print('Cluster-wise mean NDTI:', clusterStats.get('groups'));

Map.addLayer(meanClusters, {
  min: 0, max: 3,
  palette: ['1f78b4', '33a02c', 'ff7f00', 'e31a1c']
}, 'K-Means Water Types 2019-2025');

Map.addLayer(meanNdtiImg, {
  min: -0.2, max: 0.3,
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000']
}, 'Mean NDTI 2019-2025', false);

// ==================== ৭. ML স্যাম্পল এক্সপোর্ট ====================
var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD");
var modisLST = ee.ImageCollection("MODIS/061/MOD11A2");
var srtm = ee.Image("USGS/SRTMGL1_003").select('elevation').rename('Elevation');

var mlSamples = owtCollection.map(function (img) {
  var yr = ee.Number(img.get('year'));
  var sn = img.get('season');
  var start = ee.Date.fromYMD(yr, 1, 1);

  var rainSeason = chirps.filterDate(start, start.advance(1, 'year'))
      .sum().rename('Rain_Season');
  var lst = modisLST.filterDate(start, start.advance(1, 'year'))
      .select('LST_Day_1km').mean()
      .multiply(0.02).subtract(273.15).rename('LST');

  var water = img.select('MNDWI').gt(0)
      .or(img.select('NDWI').gt(0)).rename('water');

  return img.addBands([rainSeason, lst, srtm]).addBands(water)
      .toFloat()
      .sample({
        region: aoi, scale: 60, numPixels: 1000,
        geometries: false, tileScale: 8, seed: 42
      })
      .filter(ee.Filter.eq('water', 1))
      .map(function (f) { return f.set('Year', yr, 'Season', sn); });
}).flatten();

print('Total ML samples:', mlSamples.size());

Export.table.toDrive({
  collection: mlSamples,
  description: 'Gopalganj_ML_Samples_2019_2025_V1',
  folder: 'Gopalganj_Turbidity_Data',
  fileFormat: 'CSV',
  selectors: ['NDTI', 'MNDWI', 'NDWI', 'GR_Ratio', 'RedRefl',
              'B2', 'B3', 'B4', 'B8', 'B11',
              'Rain_Season', 'LST', 'Elevation', 'Year', 'Season']
});

print('Export task created. Open the Tasks tab and press Run.');

// ১. AOI
var gopalganj = ee.FeatureCollection("FAO/GAUL/2025/level2")
    .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
    .filter(ee.Filter.eq('GAUL1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('GAUL2_NAME', 'Gopalganj'));
var aoi = gopalganj.geometry();

// ২. ডাটা সোর্স
var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD");
var modisLST = ee.ImageCollection("MODIS/061/MOD11A2");
var srtm = ee.Image("USGS/SRTMGL1_003").select('elevation').rename('Elevation');

// ৩. এক্সপোর্ট স্যাম্পল লিস্ট
var sampleList = [];

// ৪. লুপ শুরু (এখানেই সব ডাটা একসাথে প্রসেস হবে)
for (var yr = 2019; yr <= 2025; yr++) {
  var seasonKeys = ['DryWinter', 'PreMonsoon', 'Monsoon', 'PostMonsoon'];
  
  seasonKeys.forEach(function(s) {
    var start = (s === 'DryWinter') ? (yr-1) + '-12-01' : (s === 'PreMonsoon' ? yr + '-03-01' : (s === 'Monsoon' ? yr + '-06-01' : yr + '-10-01'));
    var end = (s === 'DryWinter') ? yr + '-02-28' : (s === 'PreMonsoon' ? yr + '-05-31' : (s === 'Monsoon' ? yr + '-09-30' : yr + '-11-30'));

    var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterBounds(aoi)
        .filterDate(start, end)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
        .map(function(img) { return img.divide(10000); });
        
    if (col.size().getInfo() > 0) {
      var comp = col.median().clip(gopalganj);
      var ndti = comp.normalizedDifference(['B4', 'B3']).rename('NDTI');
      var mndwi = comp.normalizedDifference(['B3', 'B11']).rename('MNDWI');
      var ndwi = comp.normalizedDifference(['B3', 'B8']).rename('NDWI');
      
      var rain = chirps.filterDate(start, end).sum().rename('Rainfall');
      var lst = modisLST.filterDate(start, end).select('LST_Day_1km').mean().rename('LST');
      
      var stack = comp.addBands([ndti, mndwi, ndwi, rain, lst, srtm])
                      .addBands(ee.Image.constant(yr).rename('Year'));
      
      var water = mndwi.gt(0).or(ndwi.gt(0));
      
      var samp = stack.updateMask(water).sample({
        region: aoi, scale: 60, numPixels: 500, geometries: true
      }).map(function(f) {
        return f.set('Season', s);
      });
      
      sampleList.push(samp);
    }
  });
}

var mlSamples = ee.FeatureCollection(sampleList).flatten();

// ৫. এক্সপোর্ট
Export.table.toDrive({
  collection: mlSamples,
  description: 'Gopalganj_ML_Samples_FINAL',
  folder: 'Gopalganj_Turbidity_Data',
  fileFormat: 'CSV'
});
print('Export ready. Check Tasks tab.');