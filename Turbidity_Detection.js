// ==================== ১. কনফিগারেশন ====================
// প্রতিবার রান করার সময় শুধু এই বছরটি পরিবর্তন করবেন
var year = 2019; 

var gopalganj = ee.FeatureCollection("FAO/GAUL/2025/level2")
    .filter(ee.Filter.eq('GAUL0_NAME', 'Bangladesh'))
    .filter(ee.Filter.eq('GAUL1_NAME', 'Dhaka'))
    .filter(ee.Filter.eq('GAUL2_NAME', 'Gopalganj'));

Map.centerObject(gopalganj, 10);
Map.setOptions('SATELLITE');

// ==================== ২. ফাংশন ====================
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(1).and(scl.neq(3)).and(scl.neq(8))
      .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return image.select(['B2', 'B3', 'B4', 'B8', 'B11'])
      .divide(10000).updateMask(mask);
}

// ==================== ৩. প্রসেসিং ও এক্সপোর্ট ====================
var seasons = [
  ['DryWinter', (year-1) + '-12-01', year + '-02-28'],
  ['PreMonsoon', year + '-03-01', year + '-05-31'],
  ['Monsoon', year + '-06-01', year + '-09-30'],
  ['PostMonsoon', year + '-10-01', year + '-11-30']
];

seasons.forEach(function(s) {
  var label = s[0] + '_' + year;
  var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(gopalganj)
      .filterDate(s[1], s[2])
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
      .map(maskS2clouds);

  // স্ট্যাটাস চেক
  col.size().evaluate(function(n) {
    if (n === 0) { print(label + ': No Data'); }
    else {
      print(label + ': ' + n + ' images');
      var comp = col.median().clip(gopalganj);
      var ndti = comp.normalizedDifference(['B4', 'B3']).rename('NDTI');
      var mndwi = comp.normalizedDifference(['B3', 'B11']).rename('MNDWI');
      var ndwi = comp.normalizedDifference(['B3', 'B8']).rename('NDWI');
      var water = mndwi.gt(0).or(ndwi.gt(0)).rename('water');
      
      var stack = comp.addBands([ndti, mndwi, ndwi, water]).toFloat();
      
      // স্ট্যাটাস প্রিন্ট
      var stats = ndti.updateMask(water).reduceRegion({
        reducer: ee.Reducer.mean(), geometry: gopalganj.geometry(), scale: 30, maxPixels: 1e9
      });
      stats.get('NDTI').evaluate(function(val) { print(label + ' Mean NDTI:', val); });
      
      // এক্সপোর্ট টাস্ক
      Export.image.toDrive({
        image: stack,
        description: 'Turbidity_' + label,
        folder: 'Gopalganj_Turbidity_Data',
        region: gopalganj.geometry(),
        scale: 30,
        maxPixels: 1e13
      });
    }
  });
});

// ভিজ্যুয়ালাইজেশন: একবারে একটি ঋতু
var vizSeason = 'Monsoon'; // DryWinter, PreMonsoon, Monsoon, PostMonsoon

seasons.forEach(function(s) {
  if (s[0] !== vizSeason) { return; }

  var col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
      .filterBounds(gopalganj)
      .filterDate(s[1], s[2])
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60))
      .map(maskS2clouds);

  var comp = col.median().clip(gopalganj);
  var ndti = comp.normalizedDifference(['B4', 'B3']).rename('NDTI');
  var mndwi = comp.normalizedDifference(['B3', 'B11']);
  var water = mndwi.gt(0);

  Map.addLayer(comp, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3},
      vizSeason + '_' + year + ' True Colour', false);

  Map.addLayer(ndti.updateMask(water),
      {min: -0.2, max: 0.3, palette: ['0000ff', '00ffff', 'ffff00', 'ff0000']},
      vizSeason + '_' + year + ' NDTI');
});