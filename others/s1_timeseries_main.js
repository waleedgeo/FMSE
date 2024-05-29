// inputs
var aoi = /* color: #98ff00 */ee.Geometry.Point([68.2634891104009, 27.51666880244842]);

var startDate = '2022-01-01'
var endDate = '2023-01-01'

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
        //.filter(ee.Filter.eq('instrumentMode', 'IW'))
        //.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .select(['VV', 'VH'])
        .map(function(image) {
          var edge = image.lt(-30.0);
          var maskedImage = image.mask().and(edge.not());
          return image.updateMask(maskedImage);
        });

var filtered = s1
  .filter(ee.Filter.bounds(aoi))
  .filter(ee.Filter.date(startDate, endDate))


// Display a time-series chart
var chart = ui.Chart.image.series({
  imageCollection: filtered.select('VV'),
  region: aoi,
  reducer: ee.Reducer.mean(),
  scale: 10
}).setOptions({
      lineWidth: 1,
      title: 'S1-VV Time Series',
      interpolateNulls: true,
      vAxis: {title: 'Backscatter (VV)'},
      hAxis: {title: '', format: 'YYYY-MMM'}
    })
print(chart);

// Display a time-series chart for VH
var chart2 = ui.Chart.image.series({
    imageCollection: filtered.select('VH'),
    region: aoi,
    reducer: ee.Reducer.mean(),
    scale: 10
  }).setOptions({
        lineWidth: 1,
        title: 'S1-VH Time Series',
        interpolateNulls: true,
        vAxis: {title: 'Backscatter (VH)'},
        hAxis: {title: '', format: 'YYYY-MMM'}
        
      })
print(chart2);


// ----------


// Calculate the mean backscatter image
var meanImage = filtered.mean();

// Calculate the standard deviation of the time series
var stdImage = filtered.reduce(ee.Reducer.stdDev());

// Calculate the backscatter threshold as a fraction of the standard deviation
var backscatterThreshold = stdImage.multiply(-2); // Adjust the multiplier as needed

// Find dates with backscatter lower than the threshold
var potentialFloodDates = filtered
  .map(function(image) {
    var diffVV = image.select('VV').subtract(meanImage.select('VV'));
    var diffVH = image.select('VH').subtract(meanImage.select('VH'));
    var floodedVV = diffVV.lt(backscatterThreshold.select('VV_stdDev'));
    var floodedVH = diffVH.lt(backscatterThreshold.select('VH_stdDev'));
    var flooded = floodedVV.and(floodedVH);
    return image.set('flooded', flooded);
  })
  .filter(ee.Filter.eq('flooded', true))
  .aggregate_array('system:time_start');

// Print the dates with potential flooding
print('Dates with potential flooding:', potentialFloodDates);