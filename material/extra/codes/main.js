
var aoi = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[67.98392410136336, 26.049909335428502],
          [67.98392410136336, 25.42892423506662],
          [68.59778518534773, 25.42892423506662],
          [68.59778518534773, 26.049909335428502]]], null, false);


var startDate_pre = '2022-03-01'
var endDate_pre = '2022-05-01'


var startDate_post = '2022-08-20'
var endDate_post = '2022-09-10'

//Function to convert from dB
function toNatural(img) {
return ee.Image(10.0).pow(img.select(0).divide(10.0));
}

//Function to convert to dB
function toDB(img) {
return ee.Image(img).log10().multiply(10.0);
}

//Apllying a Refined Lee Speckle filter as coded in the SNAP 3.0 S1TBX:
//https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/RefinedLee.java
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels
   
  // convert to natural.. do not apply function on dB!
  var myimg = toNatural(img);
   
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);
   
  var mean3 = myimg.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = myimg.reduceNeighborhood(ee.Reducer.variance(), kernel3);
   
  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);
   
  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);
   
  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel);
  var sample_var = variance3.neighborhoodToBands(sample_kernel);
   
  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
   
  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());
   
  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);
   
  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);
   
  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));
   
  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);
   
  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());
   
  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
   
  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
   
  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
   
  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0],
  [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);
   
  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);
   
  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = myimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = myimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
   
  dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
   
  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
  dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
  dir_var = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
  dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  dir_var = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }
   
  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());
   
  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
   
  var b = varX.divide(dir_var);
   
  var result = dir_mean.add(b.multiply(myimg.subtract(dir_mean)));
  //return(result);
  return(img.select([]).addBands(ee.Image(toDB(result.arrayGet(0))).rename("VH")));
}


var calc_zscore = function(s1_collection_t1, s1_collection_t2, direction) {
  
  var base_mean = s1_collection_t1
    .filter(ee.Filter.equals('orbitProperties_pass', direction))
    .mean()
  
  var anom = s1_collection_t2
    .filter(ee.Filter.equals('orbitProperties_pass', direction))
    .mean()
    .subtract(base_mean)
    .set({'system:time_start': s1_collection_t2.get('system:time_start')});
  
  var base_sd = s1_collection_t1
    .filter(ee.Filter.equals('orbitProperties_pass', direction))
    .reduce(ee.Reducer.stdDev())
    .rename(['VV', 'VH', 'angle']);
      
  return anom.divide(base_sd)
    .set({'system:time_start': anom.get('system:time_start')});
}


//----------------------
var s1_t1 = ee.ImageCollection('COPERNICUS/S1_GRD')
            //.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
            //.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
            //.filter(ee.Filter.eq('instrumentMode', 'IW'))
            ////.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
            .filterBounds(aoi)
            //.select(['VH'])
            //.map(RefinedLee)
            .filterDate(startDate_pre, endDate_pre)

var s1_t2 = ee.ImageCollection('COPERNICUS/S1_GRD')
            //.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
            //.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
            //.filter(ee.Filter.eq('instrumentMode', 'IW'))
            //.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
            //.filterBounds(aoi)
            //.select(['VH'])
            .map(RefinedLee)
            .filterDate(startDate_post, endDate_post)

var zscore_des = calc_zscore(s1_t1, s1_t2, 'DESCENDING')
var zscore_asc = calc_zscore(s1_t1, s1_t2, 'ASCENDING')
var zscore = ee.ImageCollection.fromImages([zscore_des, zscore_asc]).mean()


Map.centerObject(aoi)
Map.addLayer(zscore.select('VV'), {min:-7, max:7, palette:['red', 'white', 'blue']}, 'ZScore')


// flood mask

var floods = mapFloods.mapFloods(z, parseInt(init_zvv_thd), parseInt(init_zvh_thd), 
    parseInt(init_pow_thd), parseInt(init_elev_thd), parseInt(init_slp_thd));


// this flood threshold mapping algorithm is remodified from Global Flood Mapper Model

var mapFloods = function(
  z, // ee.Image of z-score with bands "VV" and "VH"
  zvv_thd, // VV Z-score threshold
  zvh_thd, // VH Z-score threshold
  pow_thd,
  elev_thd,
  slp_thd) // Open water threshold (%)

  {
  // defaults
  if(!zvv_thd) {
    zvv_thd = -3;
  }
  if(!zvh_thd) {
    zvh_thd = -3;
  }
  if(!pow_thd) {
    pow_thd = 75;
  }
  if(!elev_thd) {
    elev_thd = 800;
  }
  if(!slp_thd) {
    slp_thd = 15;
  }
  
  // JRC water mask
  var jrc = ee.ImageCollection("JRC/GSW1_1/MonthlyHistory")
    .filterDate('2016-01-01', '2022-01-01')
  var jrcvalid = jrc.map(function(x) {return x.gt(0)}).sum();
  var jrcwat = jrc.map(function(x) {return x.eq(2)}).sum().divide(jrcvalid).multiply(100);
  var jrcmask = jrcvalid.gt(0);
  var ow = jrcwat.gte(ee.Image(pow_thd));
  
  // add elevation and slope masking
  // using fabdem
  var elevation = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")
                  .mosaic()
                  .setDefaultProjection('EPSG:3857',null,30)
  var slope = ee.Terrain.slope(elevation);

  // Classify floods
  var vvflag = z.select('VV').lte(ee.Image(zvv_thd));
  var vhflag = z.select('VH').lte(ee.Image(zvh_thd));

  var flood_class = ee.Image(0)
    .add(vvflag) 
    .add(vhflag.multiply(2))
    .where(ow.eq(1), 4)
    .rename('flood_class')
    //.updateMask(jrcmask)
    .where(elevation.gt(elev_thd).multiply(ow.neq(1)), 0)
    .where(slope.gt(slp_thd).multiply(ow.neq(1)), 0);

  return flood_class;
};


// Mask z score threshold
var fmask = mapFloods(zscore)

Map.addLayer(fmask, {}, 'flood mask')


var img = s1_t2.mean()

var field = fmask


var sample = field.addBands(image).stratifiedSample({
                      
    numPoints: 5000, 
    classBand: 'LC', 
    region: aoi, 
    scale: 10,
    seed: 5, 
    tileScale: 1.5,
    geometries: true
  
});

var histogram = sample.reduceColumns({reducer: ee.Reducer.frequencyHistogram(), selectors: ['LC']});
var samplesize = histogram.get('histogram').getInfo()
print("Sampled points per class: ", samplesize)



var sample = sample.randomColumn({seed: 1});

var training = sample.filter(ee.Filter.lt('random', 0.7));
var validation = sample.filter(ee.Filter.gte('random', 0.7));

var training = image.select(bands).sampleRegions({
  
                collection: training, 
                properties: ['LC'], 
                scale: 10
  
});

var classifier = ee.Classifier.smileRandomForest(115).train({
                  features: training, 
                  classProperty: 'LC', 
                  inputProperties: bands
  
});

var results = img.select(bands).classify(classifier);

Map.addLayer(results, {min: 1, max: 3, palette: ['C33232', '23F03F', '233CF0']}, 'Flood Classified', false);
