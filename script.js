Map.centerObject(roi, 10);
Map.setOptions('satellite');

Map.addLayer(roi.style({
  color: 'black', fillColor: 'ffffff00'
}), {}, 'ROI');

var startDate = '2023-01-01';
var endDate = '2023-12-31';

                     //*****Optical Data****//

var preImages = s2.filterDate(startDate, endDate)
              .filterBounds(roi)
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 0.1))
              .select(['B.*']);
  
var scaledImages = preImages.map(function(image) {
   var scaled = image.multiply(0.0001).copyProperties(image, image.propertyNames());
   
   return scaled});
   
var opticalIndices = scaledImages.map(function(scaled){
   var ndvi = scaled.normalizedDifference(['B8', 'B4']).rename('ndvi');
   var ndri = scaled.normalizedDifference(['B6', 'B5']).rename('ndri');
   var nsmi = scaled.normalizedDifference(['B11', 'B8']).rename('nsmi');
   var cire = (scaled.select(['B4']).divide(scaled.select(['B11']))).subtract(-1).rename('cire');
   var srre = scaled.select(['B7']).divide(scaled.select(['B5'])).rename('srre');
   var mtci = (scaled.select(['B6']).subtract(scaled.select(['B5']))).divide((scaled.select(['B5']).subtract(scaled.select(['B4'])))).rename('mtci');
   var pvi = scaled.normalizedDifference(['B8', 'B5']).rename('pvi');
   
  return scaled.addBands([ndvi, ndri, nsmi, cire, srre, mtci, pvi]);
}).median().clip(roi);
   
print(opticalIndices, 'Optical Indices');

Map.addLayer(scaledImages.median().clip(roi), {
  min: 0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2']
}, 'S2 Composite', false);

Map.addLayer(opticalIndices, {
  min: 0,
  max: 1
}, 'Optical Indices', false);

//LandCover Dataset from Dynamic World

var lulc = dynamicWorld.filterBounds(roi)
        .filterDate(startDate, endDate)
        .select(['label'])
        .median()
        .clip(roi);

//Satellite Embeddings from Alpha Earth

var embed = se.filterBounds(roi)
          .filterDate(startDate, endDate)
          .median()
          .clip(roi);

     // ***** Radar Data ***//

//Sentinel-1 ** C- Band **
var s1_sar = s1.filterBounds(roi)
        .filterDate(startDate, endDate)
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'));

var s1_indices = s1_sar.map(function(image) {
  
  var sigma = ee.Image.constant(10).pow(image.divide(10));
  
  var rvi = sigma.expression(
    '(8 * vh)/(vh + vv)', 
     { 'vh': sigma.select('VH'),
       'vv': sigma.select('VV')
     }).rename('rvi');
     
  var ndi = sigma.expression(
    'log10 (10 * vv * vh)',
    {'vh': sigma.select('VH'),
    'vv': sigma.select('VV')
    }).rename('ndi');
    
  var ri = sigma.select(['VV']).divide(sigma.select(['VH'])).rename('ri');
  
  var combined = ee.Image.cat([rvi, ndi, ri]);
  
  var smoothed = combined.focalMedian(30, 'square', 'meters');
  
  return smoothed.copyProperties(image, ['system:time_start']);
     
}).mean().clip(roi);

Map.addLayer(s1_sar.mean().clip(roi), {}, 'S1 dB', false);
Map.addLayer(s1_indices, {}, 'S1 RVI', false);

print(s1_indices, 'S1 Indices');

//ALOS PALSAR ** L- Band

var alos_sar = alos.filterBounds(roi)
          .filterDate(startDate, endDate)
          .median()
          .clip(roi);

var alos_ri = alos_sar.select(['HH']).add(alos_sar.select(['HV'])).rename('alos_ri');

var alos_rvi = alos_sar.expression(
  '(4 * hv) / (hh + hv)',
  {'hh': alos_sar.select(['HH']),
   'hv': alos_sar.select(['HV'])
}).rename('alos_rvi');

var alos_ndi = alos_sar.normalizedDifference(['HH', 'HV']).rename('alos_ndi');

var alos_indices = ee.Image.cat([alos_ri, alos_ndi, alos_rvi]);

Map.addLayer(alos_sar, {}, 'ALOS PALSAR', false);
Map.addLayer(alos_indices, {}, 'ALOS Indices', false);

print(alos_indices, 'ALOS Indices');

///*** Terrain Data ***///

var elevation = srtm.clip(roi).rename('elevation');
var slope = ee.Terrain.slope(elevation).rename('slope');
var aspect = ee.Terrain.aspect(elevation).rename('aspect');

//Merged Features

var featureImage = opticalIndices
      .addBands(s1_indices)
      .addBands(alos_indices)
      .addBands([elevation, slope, aspect])
      .addBands(embed)
      .addBands(lulc)
      .reproject('EPSG:4326', null, 10);
      
print(featureImage, 'Feature Images');

// GEDI L2A - Canopy Height Foot Prints

var gedi_rh = gedi.filterBounds(roi)
              .filterDate(startDate, endDate)
              .map(function(image) {
               return image.updateMask(image.select(['quality_flag']).eq(1))
                .updateMask(image.select(['degrade_flag']).eq(0));
              })
              .select(['rh98'])
              .mean()
              .clip(roi);

Map.addLayer(gedi_rh, {}, 'GEDI RH 98');

var samplePoints = gedi_rh.addBands(featureImage).sample({
  region: roi, 
  scale: 100, 
  projection: 'EPSG:4326', 
  numPixels: 200000, 
  tileScale: 12, 
  geometries: true
});

Export.table.toDrive({
  collection: samplePoints, 
  description: 'Sample Points', 
  fileFormat: 'SHP'
});

var split = samplePoints.randomColumn('random', 500);

var trainingData = split.filter(ee.Filter.lte('random', 0.7));
var validationData = split.filter(ee.Filter.gt('random', 0.7));

var trainingData = trainingData.map(function(feature) {
  
  var properties = feature.toDictionary();
  properties = properties.remove(['random']);
  return ee.Feature(feature.geometry(), properties);
});

print(trainingData.size(), 'Training Dataset');

var validationData = validationData.map(function(feature) {
  
  var properties = feature.toDictionary();
  properties = properties.remove(['random']);
  
  return ee.Feature(feature.geometry(), properties);
});

print(validationData.size(), 'Validation Dataset');

var predictionBands = featureImage.bandNames();
print(predictionBands, 'Prediction Bands');

var classifier = ee.Classifier.smileRandomForest(500).train({
  features: trainingData, 
  classProperty: 'rh98', 
  inputProperties: predictionBands
}).setOutputMode('REGRESSION');

var canopy_height = featureImage.classify(classifier).clip(roi).rename('canopy_height');

Map.addLayer(canopy_height, {}, 'Canopy Height', false);

Export.image.toDrive({
  image: canopy_height, 
  description: 'Canopy Height',
  region: roi, 
  scale: 10, 
  crs: 'EPSG:32643',
  maxPixels: 1e13
});

//******Model Evaluation*******//

var validation = canopy_height.sampleRegions({
  collection: validationData, 
  properties: ['rh98'], 
  scale: 100, 
  tileScale: 16, 
  geometries: true
}).filter(ee.Filter.notNull(['rh98']));

//Correlation value calculation

var actualVsPredicted = validation.reduceColumns({
  reducer: ee.Reducer.pearsonsCorrelation(), 
  selectors: ['rh98', 'canopy_height']
});

var correlation = ee.Number(actualVsPredicted.get('correlation'));
print(correlation, 'Correlation');
print(correlation.pow(2), 'R Square Value')

//calculate rmse

var rmse = validation.map(function(f) {
  var actual = ee.Number(f.get('rh98'));
  var predicted =ee.Number(f.get('canopy_height'));
  
  return ee.Feature(null, {
    'diff_squared': actual.subtract(predicted).pow(2)
  });
}).aggregate_mean('diff_squared');

var rmse = ee.Number(rmse).sqrt();
print(rmse, 'RMSE Value');

//Calculation of MAE

var mae = validation.map(function(f) {
   var actual = ee.Number(f.get('rh98'));
  var predicted = ee.Number(f.get('canopy_height'));
  return ee.Feature(null, {'abs_diff': actual.subtract(predicted).abs()});
}).aggregate_mean('abs_diff');

print(mae, 'MAE value')

var modelMetrics = ee.Feature(null, 
    {'RMSE': rmse,
    'MAE': mae,
    'Correlation': correlation,
    'R Square': correlation.pow(2)
});

print(modelMetrics);

Export.table.toDrive({
  collection: modelMetrics, 
  description: 'MODEL METRICS',
  fileFormat: 'CSV'
});









