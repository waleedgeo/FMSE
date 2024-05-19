
// User defined variables
var pre_startDate = '1990-01-01'
var pre_endDate = '1991-05-01'
var post_startDate = '2021-01-01'
var post_endDate = '2022-05-01'
// filter months data 
var startMonth = 1;
var endMonth = 4;
var cloud_lte = 30
//var aoi = ee.FeatureCollection("projects/sylhet-flood/assets/sylhet_div");
//var aoi1 = ee.FeatureCollection("projects/autowet/assets/aoi1_spain");

var aoi = geometry
// Fixed variables
Map.centerObject(aoi, 12)
// **** Map Style  *********
// Styling map

var mapStyle = [
    {
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#f5f5f5"
        }
      ]
    },
    {
      "elementType": "labels",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "elementType": "labels.icon",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#616161"
        }
      ]
    },
    {
      "elementType": "labels.text.stroke",
      "stylers": [
        {
          "color": "#f5f5f5"
        }
      ]
    },
    {
      "featureType": "administrative",
      "elementType": "geometry",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "administrative.land_parcel",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "administrative.land_parcel",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#bdbdbd"
        }
      ]
    },
    {
      "featureType": "administrative.neighborhood",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "poi",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "poi",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#eeeeee"
        }
      ]
    },
    {
      "featureType": "poi",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#757575"
        }
      ]
    },
    {
      "featureType": "poi.park",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#e5e5e5"
        }
      ]
    },
    {
      "featureType": "poi.park",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#9e9e9e"
        }
      ]
    },
    {
      "featureType": "road",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "road",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#ffffff"
        }
      ]
    },
    {
      "featureType": "road",
      "elementType": "labels.icon",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "road.arterial",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#757575"
        }
      ]
    },
    {
      "featureType": "road.highway",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#dadada"
        }
      ]
    },
    {
      "featureType": "road.highway",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#616161"
        }
      ]
    },
    {
      "featureType": "road.local",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#9e9e9e"
        }
      ]
    },
    {
      "featureType": "transit",
      "stylers": [
        {
          "visibility": "off"
        }
      ]
    },
    {
      "featureType": "transit.line",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#e5e5e5"
        }
      ]
    },
    {
      "featureType": "transit.station",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#eeeeee"
        }
      ]
    },
    {
      "featureType": "Wetlands",
      "elementType": "geometry",
      "stylers": [
        {
          "color": "#c9c9c9"
        }
      ]
    },
    {
      "featureType": "Wetlands",
      "elementType": "labels.text.fill",
      "stylers": [
        {
          "color": "#9e9e9e"
        }
      ]
    }
  ]
  
  Map.setOptions('mapStyle', {mapStyle: mapStyle});
//------------- Functions ---------------------------------------------
// Applies scaling factors for Landsat
function scalefactor (img) {
  var bands = img.bandNames()
  var opticalBands = img.select(bands).multiply(0.0000275).add(-0.2);
  return img.addBands(opticalBands, null, true)
}
//removing clouds
function mask_cloud(img) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = img.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return img.updateMask(mask);
}

// Function to rename Landsat bands to simple b,g,r, nir, swir names
var l8rename = function(img){
  var img = img.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'QA_PIXEL'])
  return img.rename(['b', 'g', 'r', 'nir', 'swir', 'QA_PIXEL'])
} 

var l75rename = function(img){
  var img = img.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'QA_PIXEL'])
  return img.rename(['b', 'g', 'r', 'nir', 'swir', 'QA_PIXEL'])
} 


// Spectral indices function
var si = function(img){
  var mndwi = img.normalizedDifference(['g', 'swir']).rename('mndwi');
  var ndbi = img.normalizedDifference(['swir', 'nir']).rename('ndbi')
  var ndvi = img.normalizedDifference(['nir', 'r']).rename('ndvi')
  var lswi = img.normalizedDifference(['nir', 'swir']).rename('lswi');
  var evi =  img.expression('2.5 * (nir - r) / (nir + 6.0 * r - 7.5 * b + 1.0)', {
                            nir: img.select('nir'), 
                            r: img.select('r'), 
                            b: img.select('b')
  
}).rename('evi');
  return img.addBands(mndwi).addBands(ndbi).addBands(ndvi).addBands(lswi).addBands(evi);
}

// Importing Landsat 8 for images between 2013-current, L7 for 2012-13,and Landsat 55 for 1984-2012
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
            .filterBounds(aoi)
            .filter(ee.Filter.calendarRange(startMonth, endMonth, 'month'))
            .map(l8rename);
var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
            .filterBounds(aoi)
            .filter(ee.Filter.calendarRange(startMonth, endMonth, 'month'))
            .map(l75rename);
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterDate('2012-05-05', '2013-03-18')
            .filterBounds(aoi)
            .filter(ee.Filter.calendarRange(startMonth, endMonth, 'month'))
            .map(l75rename);

// Merging all image collection into oone
var dataset = l5.merge(l7).merge(l8);
// Filtering collection (year/month/cloud_cover)
var pre_dataset = dataset.filterDate(pre_startDate, pre_endDate)
                    .filter(ee.Filter.lte('CLOUD_COVER', cloud_lte))

var pre_dataset = pre_dataset.map(mask_cloud)
var pre_dataset = pre_dataset.map(scalefactor)
var pre_dataset = pre_dataset.map(si)


var post_dataset = dataset.filterDate(post_startDate, post_endDate)
                    .filter(ee.Filter.lte('CLOUD_COVER', cloud_lte))
var post_dataset = post_dataset.map(mask_cloud)
var post_dataset = post_dataset.map(scalefactor)
var post_dataset = post_dataset.map(si)

var pre_dataset_size = pre_dataset.size()
print("Total images in Pre-Dataset: ", pre_dataset_size)
var post_dataset_size = post_dataset.size()
print("Total images in Post-Dataset: ", post_dataset_size)
var pre_image = pre_dataset.mean().clip(aoi)
var post_image = post_dataset.mean().clip(aoi)
// ****************** EXTRA  ************************************************************
// --------- Defining Visualization parameters  --------------

var palettes = require('users/gena/packages:palettes');
var vis_rgb = {
  bands: ['r', 'g', 'b'],
  min: 0.0,
  max: 0.3,

}
var vis_fcc = {
  bands: ['nir', 'r', 'g'],
  min: 0.0,
  max: 0.3,
}
var vis_mndwi = {
  min: -0.4,
  max: 0.6,
  palette:['#f1eef6','#d0d1e6','#a6bddb','#74a9cf','#2b8cbe','#045a8d']
  
}
var vis_ndbi = {
  min: -4,
  max: 0.14,
  palette: palettes.colorbrewer.OrRd[7]
}
var vis_ndvi = {
  min: -0.2,
  max: 0.7,
  palette: palettes.colorbrewer.RdYlGn[7]
}
// **************************************************
// defining SI individually (PRE)
var pre_mndwi = pre_image.select('mndwi')
var pre_ndbi = pre_image.select('ndbi')
var pre_ndvi = pre_image.select('ndvi')
var pre_lswi = pre_image.select('lswi')
var pre_evi = pre_image.select('evi')

// defining SI individually (Post)
var post_mndwi = post_image.select('mndwi')
var post_ndbi = post_image.select('ndbi')
var post_ndvi = post_image.select('ndvi')
var post_lswi = post_image.select('lswi')
var post_evi = post_image.select('evi')

//-----------------------------------------------------------------------
// ********** Creating Threshold Mask  *********************************
// ******* For Pre **************
// Wetlands Mask
var pre_pureWetlands = pre_mndwi.gt(pre_evi.lt(0.1)).or(pre_mndwi.gt(pre_ndvi));
var pre_Wetlands = pre_pureWetlands.updateMask(pre_pureWetlands);
Map.addLayer(pre_Wetlands, vis_mndwi, 'Pre-Wetlands', false)
// -------- veg mask
var pre_pureVeg = (pre_evi.gte(0.1), pre_ndvi.gte(0.3).and(pre_lswi.gt(0)));
var pre_veg = pre_pureVeg.updateMask(pre_pureVeg)
Map.addLayer(pre_veg, vis_ndvi, 'Pre-Veg', false)
// impervious mask
var pre_pureImp = pre_ndbi.gte(0.1)
var pre_impervious = pre_pureImp.updateMask(pre_pureImp)
Map.addLayer(pre_impervious, vis_ndbi, 'Pre-Impervious Surfaces', false)

// ******* For Post **************
// Wetlands Mask
var post_pureWetlands = post_mndwi.gt(post_evi.lt(0.1)).or(post_mndwi.gt(post_ndvi));
var post_Wetlands = post_pureWetlands.updateMask(post_pureWetlands);
Map.addLayer(post_Wetlands, vis_mndwi, 'post-Wetlands', false)
// -------- veg mask
var post_pureVeg = (post_evi.gte(0.1), post_ndvi.gte(0.3).and(post_lswi.gt(0)));
var post_veg = post_pureVeg.updateMask(post_pureVeg)
Map.addLayer(post_veg, vis_ndvi, 'post-Veg', false)
// impervious mask
var post_pureImp = post_ndbi.gte(0.1)
var post_impervious = post_pureImp.updateMask(post_pureImp)
Map.addLayer(post_impervious, vis_ndbi, 'post-Impervious Surfaces', false)

// ---------- Classification  ------------------


//**************** Pre  ***********************
// vectors
var pre_imp_c = pre_impervious
var pre_veg_c = pre_veg
var pre_wet_c = pre_Wetlands
var pre_imp_vector = pre_imp_c.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'pre_impervious',
  bestEffort: true,
  maxPixels: 1e12
  
});
var pre_veg_vector = pre_veg_c.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'pre_veg',
  bestEffort: true,
  maxPixels: 1e12
  
});
var pre_wet_vector = pre_wet_c.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'pre_wetlands',
  bestEffort: true,
  maxPixels: 1e12
  
});

//**************** post  ***********************
// vectors
var post_imp_c = post_impervious
var post_veg_c = post_veg
var post_wet_c = post_Wetlands
var post_imp_vector = post_imp_c.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'post_impervious',
  bestEffort: true,
  maxPixels: 1e12
  
});
var post_veg_vector = post_veg_c.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'post_veg',
  bestEffort: true,
  maxPixels: 1e12
  
});
var post_wet_vector = post_wet_c.reduceToVectors({
  geometry: aoi,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'post_wetlands',
  bestEffort: true,
  maxPixels: 1e12
  
});
var bands = ['b', 'g', 'r', 'nir', 'swir', 'mndwi', 'ndbi', 'ndvi', 'lswi', 'evi']

var pre_IMPERVIOUS = pre_image.select('b').gt(0).clip(pre_imp_vector).multiply(1).rename('LC').int();
var pre_VEG = pre_image.select('b').gt(0).clip(pre_veg_vector).multiply(2).rename('LC').int();
var pre_WET = pre_image.select('b').gt(0).clip(pre_wet_vector).multiply(3).rename('LC').int()


var pre_field = ee.ImageCollection([pre_IMPERVIOUS, pre_VEG, pre_WET]).mosaic();

/////////////  TAKE A STRATIFIED SAMPLE    ////////////

var pre_sample = pre_field.addBands(pre_image).stratifiedSample({
                      
                      numPoints: 5000, 
                      classBand: 'LC', 
                      region: aoi, 
                      scale: 100,
                      seed: 5, 
                      tileScale: 1.5,
                      geometries: true
                    
                  });

var pre_histogram = pre_sample.reduceColumns({reducer: ee.Reducer.frequencyHistogram(), selectors: ['LC']});
var pre_samplesize = pre_histogram.get('histogram').getInfo()
print("Pre-Sampled points per class: ", pre_samplesize)

/////////////  RANDOMIZE TRAINING AND VALIDATION DATA    ////////////
var pre_sample = pre_sample.randomColumn({seed: 1});
var pre_training = pre_sample.filter(ee.Filter.lt('random', 0.7));
var pre_validation = pre_sample.filter(ee.Filter.gte('random', 0.7));

var pre_training = pre_image.select(bands).sampleRegions({
  
                collection: pre_training, 
                properties: ['LC'], 
                scale: 30
  
});

var pre_classifier = ee.Classifier.smileRandomForest(115).train({
                  features: pre_training, 
                  classProperty: 'LC', 
                  inputProperties: bands
  
});

var pre_LULC = pre_image.select(bands).classify(pre_classifier);

Map.addLayer(pre_LULC, {min: 1, max: 3, palette: ['C33232', '23F03F', '233CF0']}, 'Pre-LULC Classified', false);

// ********* Classification Legend  *************

// *************************** Legend  *******************
var lulc_names = ['Impervious', 'Vegetation', 'Wetlands']
var lulc_colors = ['#C33232', '#23F03F', '#233CF0']
// Setting Legend
// Create the panel for the legend items.
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Create and add the legend title.
var legendTitle = ui.Label({
  value: 'Land-use (Wetlands)',// enter your title here
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
  }
});
legend.add(legendTitle);

var loading = ui.Label('Loading legend...', {margin: '2px 0 4px 0'});
legend.add(loading);

// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
  // Create the label that is actually the colored box.
  var colorBox = ui.Label({
    style: {
      backgroundColor: color,
      // Use padding to give the box height and width.
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });

  // Create the label filled with the description text.
  var description = ui.Label({
    value: name,
    style: {margin: '0 0 4px 6px'}
  });

  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};

// Get the list of palette colors and class names from the image.
pre_LULC.toDictionary().select(['classification' + ".*"]).evaluate(function(result) {
  var palette = lulc_colors;// classes colour
  var names = lulc_names;// classes name
  loading.style().set('shown', false);

  for (var i = 0; i < names.length; i++) {
    legend.add(makeRow(palette[i], names[i]));
  }
  });

Map.add(legend);
// ************** post classification  *********************
var post_IMPERVIOUS = post_image.select('b').gt(0).clip(post_imp_vector).multiply(1).rename('LC').int();
var post_VEG = post_image.select('b').gt(0).clip(post_veg_vector).multiply(2).rename('LC').int();
var post_WET = post_image.select('b').gt(0).clip(post_wet_vector).multiply(3).rename('LC').int()


var post_field = ee.ImageCollection([post_IMPERVIOUS, post_VEG, post_WET]).mosaic();

/////////////  TAKE A STRATIFIED SAMPLE    ////////////

var post_sample = post_field.addBands(post_image).stratifiedSample({
                      
                      numPoints: 5000, 
                      classBand: 'LC', 
                      region: aoi, 
                      scale: 100,
                      seed: 5, 
                      tileScale: 1.5,
                      geometries: true
                    
                  });

var post_histogram = post_sample.reduceColumns({reducer: ee.Reducer.frequencyHistogram(), selectors: ['LC']});

var post_samplesize = post_histogram.get('histogram').getInfo()
print("Post-Sampled points per class: ", post_samplesize)

/////////////  RANDOMIZE TRAINING AND VALIDATION DATA    ////////////
var post_sample = post_sample.randomColumn({seed: 1});
var post_training = post_sample.filter(ee.Filter.lt('random', 0.7));
var post_validation = post_sample.filter(ee.Filter.gte('random', 0.7));

var post_training = post_image.select(bands).sampleRegions({
  
                collection: post_training, 
                properties: ['LC'], 
                scale: 30
  
});

var post_classifier = ee.Classifier.smileRandomForest(115).train({
                  features: post_training, 
                  classProperty: 'LC', 
                  inputProperties: bands
  
});

var post_LULC = post_image.select(bands).classify(post_classifier);

Map.addLayer(post_LULC, {min: 1, max: 3, palette: ['C33232', '23F03F', '233CF0']}, 'post-LULC Classified', false);

// Explain classifier
//print("Pre Classifier Explain: ", pre_classifier.explain())
//print("Post Classifier Explain: ", post_classifier.explain())
//  ************************** Accuracy************

// validation classifier using testing samples - Pre
var pre_valc = pre_validation.classify(pre_classifier);
// accuracy matrix
var pre_em = pre_valc.errorMatrix('LC','classification');// Error Matrix
var pre_oa = pre_em.accuracy(); // Overall accuracy
var pre_ks = pre_em.kappa(); // Kappa statistic
var pre_ua = pre_em.consumersAccuracy().project([1]); // Consumer's accuracy
var pre_pa = pre_em.producersAccuracy().project([0]); // Producer's accuracy
var pre_f1 = (pre_ua.multiply(pre_pa).multiply(2.0)).divide(pre_ua.add(pre_pa)); // F1-statistic

// Create a Feature with null geometry and the value we want to export.
// Use .array() to convert Confusion Matrix to an Array so it can be
// exported in a CSV file
var pre_emfc = ee.FeatureCollection([
  ee.Feature(null, {
    'accuracy': pre_oa,
    'Kappa': pre_ks,
    'F1s': pre_f1,
    'UA': pre_ua,
    'PA': pre_pa,
    'ErrorMatrix': pre_em.array()
  })
  ])
print('pre_export accuracy', pre_emfc)

// validation classifier using testing samples - post
var post_valc = post_validation.classify(post_classifier);
// accuracy matrix
var post_em = post_valc.errorMatrix('LC','classification');// Error Matrix
var post_oa = post_em.accuracy(); // Overall accuracy
var post_ks = post_em.kappa(); // Kappa statistic
var post_ua = post_em.consumersAccuracy().project([1]); // Consumer's accuracy
var post_pa = post_em.producersAccuracy().project([0]); // Producer's accuracy
var post_f1 = (post_ua.multiply(post_pa).multiply(2.0)).divide(post_ua.add(post_pa)); // F1-statistic

// preparing error matrix to export
var post_emfc = ee.FeatureCollection([
  ee.Feature(null, {
    'accuracy': post_oa,
    'Kappa': post_ks,
    'F1s': post_f1,
    'UA': post_ua,
    'PA': post_pa,
    'ErrorMatrix': post_em.array()
  })
  ])
print('post_export accuracy', post_emfc)  
//************ exporting LULC images - temp ***
/*
Export.image.toAsset({
  
      image: pre_LULC,
      description: "aoi1_Pre_LULC", 
      region: aoi,
      scale: 30,
      maxPixels: 1e12
  
})

// Post
Export.image.toAsset({
  
      image: post_LULC,
      description: "aoi1_post_LULC", 
      region: aoi,
      scale: 30,
      maxPixels: 1e12
})


Export.table.toDrive({
  
        collection: pre_emfc, 
        description: 'Pre_EM_AOI1', 
        folder:'autowet',
        fileNamePrefix: 'pre_accuracy_aoi1',
        fileFormat: 'CSV', 
  
});

Export.table.toDrive({
  
        collection: post_emfc, 
        description: 'Post_EM_AOI1', 
        folder:'autowet',
        fileNamePrefix: 'post_accuracy_aoi1',
        fileFormat: 'CSV', 
  
});

*/
// ***********************  Area Calculation Part *************************

//********* Pre ********* *
var pre_names = ['pre_impervious', 'pre_veg','pre_wetlands'];
var pre_renamed = pre_LULC.eq([1, 2, 3]).rename(pre_names);

var pre_area = pre_renamed.multiply(ee.Image.pixelArea().divide(10000));
var pre_classarea = pre_area.reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale:30, maxPixels: 1e13});
var pre_areatot = ee.Number(pre_classarea);
var a = ee.Array(pre_classarea.get('pre_impervious'));
var b = ee.Array(pre_classarea.get('pre_veg'));
var c = ee.Array(pre_classarea.get('pre_wetlands'));
var pre_array = ee.List([a, b, c]);

var pre_chart = ui.Chart.array.values(pre_array, 0, pre_names).setChartType('PieChart').setOptions({is3D: true, title:'Pre Area (km2)'}); 
print(pre_chart);

//********* post *********
var post_names = ['post_impervious', 'post_veg','post_wetlands'];
var post_renamed = post_LULC.eq([1, 2, 3]).rename(post_names);

var post_area = post_renamed.multiply(ee.Image.pixelArea().divide(10000));
var post_classarea = post_area.reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale:30, maxPixels: 1E13});
  
var post_areatot = ee.Number(post_classarea);
var w = ee.Array(post_classarea.get('post_impervious'));
var x = ee.Array(post_classarea.get('post_veg'));
var y = ee.Array(post_classarea.get('post_wetlands'));
var post_array = ee.List([w, x, y]);

var post_chart = ui.Chart.array.values(post_array, 0, post_names).setChartType('PieChart').setOptions({is3D: true, title:'Post Area (km2)'}); 
print(post_chart);


// ***************  Change Detection  *******************************

var wetlandsAreaC = ee.Number(post_classarea.get('post_wetlands')).subtract(ee.Number(pre_classarea.get('pre_wetlands')))
var wetlandsAreaC = wetlandsAreaC.getInfo().toFixed(2)
print('Wetlands Area Change Km2: ', wetlandsAreaC)

var merged = pre_LULC.multiply(100).add(post_LULC).rename('transitions')

/*
101 = impervious to impervious
102 = impervious to vegetation
103 = impervious to Wetlands
201 = Vegetation to Impervious
202 = Vegetation to Vegetation
203 = Vegetation to Wetlands
301 = Wetlands to impervious
302 = Wetlands to vegetation
303 = Wetlands to Wetlands

*/
// define masks for LULC CD with focus on Wetlands

var imp_Wetlands = merged.eq(103)
var imp_Wetlands = ee.Image(imp_Wetlands.updateMask(imp_Wetlands).multiply(1)).rename('CD').toInt8()
var veg_Wetlands = merged.eq(203)
var veg_Wetlands = ee.Image(veg_Wetlands.updateMask(veg_Wetlands).multiply(2)).rename('CD').toInt8()
var Wetlands_imp = merged.eq(301)
var Wetlands_imp = ee.Image(Wetlands_imp.updateMask(Wetlands_imp).multiply(3)).rename('CD').toInt8()
var Wetlands_veg = merged.eq(302)
var Wetlands_veg = ee.Image(Wetlands_veg.updateMask(Wetlands_veg).multiply(4)).rename('CD').toInt8()
var Wetlands_unch = merged.eq(303)
var Wetlands_unch = ee.Image(Wetlands_unch.updateMask(Wetlands_unch).multiply(5)).rename('CD').toInt8()

var imageC_cd = ee.ImageCollection.fromImages([imp_Wetlands, veg_Wetlands, Wetlands_imp, Wetlands_veg, Wetlands_unch])
var final_merge = imageC_cd.mosaic()
var cd_colors = ['#4daf4a','#ff7f00','#e41a1c','#984ea3','#377eb8']
var cd_names = ['Imp to Wetlands', 'Veg to Wetlands', 'Wetlands to Imp', 'Wetlands to Veg', 'Wetlands Unchanged']
Map.addLayer(final_merge, {min:1, max:5, palette:cd_colors}, 'Change Transition')


// *************************** legend2  *******************

// Setting legend2
// Create the panel for the legend2 items.
var legend2 = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Create and add the legend2 title.
var legend2Title = ui.Label({
  value: 'Change Transition',// enter your title here
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
  }
});
legend2.add(legend2Title);

var loading = ui.Label('Loading legend2...', {margin: '2px 0 4px 0'});
legend2.add(loading);

// Creates and styles 1 row of the legend2.
var makeRow = function(color, name) {
  // Create the label that is actually the colored box.
  var colorBox = ui.Label({
    style: {
      backgroundColor: color,
      // Use padding to give the box height and width.
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });

  // Create the label filled with the description text.
  var description = ui.Label({
    value: name,
    style: {margin: '0 0 4px 6px'}
  });

  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};

var lulc_name = cd_names;

var lulc_color = cd_colors
// Get the list of palette colors and class names from the image.
final_merge.toDictionary().select(['CD' + ".*"]).evaluate(function(result) {
  var palette = lulc_color;// classes colour
  var names = lulc_name;// classes name
  loading.style().set('shown', false);

  for (var i = 0; i < names.length; i++) {
    legend2.add(makeRow(palette[i], names[i]));
  }
  });

Map.add(legend2);


// ********** Evaluate area of Change Transitions  **********

var cd_names = cd_names
var cd_renamed = final_merge.eq([1, 2, 3, 4, 5]).rename(cd_names);

var cd_area = cd_renamed.multiply(ee.Image.pixelArea().divide(10000));
var cd_classarea = cd_area.reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi, scale:30, maxPixels: 1e13});
  
var cd_areatot = ee.Number(cd_classarea);
//['Imp to Wetlands', 'Veg to Wetlands', 'Wetlands to Imp', 'Wetlands to Veg', 'Wetlands Unchanged']
var a = ee.Array(cd_classarea.get('Imp to Wetlands'));
var b = ee.Array(cd_classarea.get('Veg to Wetlands'));
var c = ee.Array(cd_classarea.get('Wetlands to Imp'));
var d = ee.Array(cd_classarea.get('Wetlands to Veg'));
var e = ee.Array(cd_classarea.get('Wetlands Unchanged'));
var cd_array = ee.List([a, b, c, d, e]);

var cd_chart = ui.Chart.array.values(cd_array, 0, cd_names).setChartType('ColumnChart').setOptions({is3D: true, title:'Change Transitions (km2)'}); 
print(cd_chart);
