# Part-1: Flood Mapping

import ee
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
import numpy as np
import matplotlib.pyplot as plt



def boundary(feature):
    '''
    This function takes a feature and returns the boundary of the feature
    
    Args:
    feature: ee.Feature
    returns:
    boundary: ee.Feature
    
    '''
    
    bbox = feature.geometry().bounds()
    boundary = ee.Geometry.Polygon(bbox.coordinates().get(0))
    
    return boundary
    


def get_s1_col(date, days, aoi):
    """
    Fetch Sentinel-1 Image Collection based on the given date and filters.
    
    Parameters:
    date (ee.Date): The starting date for filtering the images.
    days (int): Number of days for filtering the images.
    aoi (ee.Geometry): Area of Interest for filtering the images.
    
    Returns:
    ee.ImageCollection: Filtered Sentinel-1 Image Collection.
    """
    filters = [
        ee.Filter.listContains("transmitterReceiverPolarisation", "VV"),
        ee.Filter.listContains("transmitterReceiverPolarisation", "VH"),
        ee.Filter.Or(ee.Filter.equals("instrumentMode", "IW"), ee.Filter.equals("instrumentMode", "SM")),
        ee.Filter.bounds(aoi),
        ee.Filter.eq('resolution_meters', 10),
        ee.Filter.date(date, date.advance(days + 1, 'day'))
    ]
    
    return ee.ImageCollection('COPERNICUS/S1_GRD').filter(filters)

def calc_zscore(s1_pre, s1_post, direction):
    """
    Calculate Z-score for the given direction (ascending/descending).

    Parameters:
    s1_pre (ee.ImageCollection): Pre-flood image collection.
    s1_post (ee.ImageCollection): Post-flood image collection.
    direction (str): Orbit direction (ASCENDING or DESCENDING).

    Returns:
    ee.Image: Z-score image.
    """
    base_mean = s1_pre.filter(ee.Filter.equals('orbitProperties_pass', direction)).mean()
    anom = s1_post.filter(ee.Filter.equals('orbitProperties_pass', direction)).mean().subtract(base_mean)
    base_sd = s1_pre.filter(ee.Filter.equals('orbitProperties_pass', direction)).reduce(ee.Reducer.stdDev()).rename(['VV', 'VH'])
    return anom.divide(base_sd).set({'system:time_start': s1_post.get('system:time_start')})

def calculate_zscore(s1_pre, s1_post, aoi):
    """
    Calculate combined Z-score for both ascending and descending orbits.

    Parameters:
    s1_pre (ee.ImageCollection): Pre-flood image collection.
    s1_post (ee.ImageCollection): Post-flood image collection.
    aoi (ee.Geometry): Area of Interest.

    Returns:
    ee.Image: Combined Z-score image.
    """
    asc = ee.Filter.eq("orbitProperties_pass", "ASCENDING")
    des = ee.Filter.eq("orbitProperties_pass", "DESCENDING")
    
    cond_asc = s1_pre.filter(asc).size().gt(0).And(s1_post.filter(asc).size().gt(0))
    cond_des = s1_pre.filter(des).size().gt(0).And(s1_post.filter(des).size().gt(0))

    if cond_asc.getInfo():
        return calc_zscore(s1_pre, s1_post, 'ASCENDING')
    elif cond_des.getInfo():
        return calc_zscore(s1_pre, s1_post, 'DESCENDING')
    else:
        zscore_des = calc_zscore(s1_pre, s1_post, 'DESCENDING')
        zscore_asc = calc_zscore(s1_pre, s1_post, 'ASCENDING')
        return ee.ImageCollection.fromImages([zscore_des, zscore_asc]).mean().clip(aoi)

def map_floods(z, aoi, zvv_thd, zvh_thd, pow_thd, elev_thd, slp_thd, under_estimate):
    
    """
    Generate flood mask based on Z-score and various thresholds.

    Parameters:
    z (ee.Image): Z-score image.
    aoi (ee.Geometry): Area of Interest.
    zvv_thd (float): Threshold for VV band Z-score.
    zvh_thd (float): Threshold for VH band Z-score.
    pow_thd (float): Threshold for open water percentage.
    elev_thd (float): Elevation threshold.
    slp_thd (float): Slope threshold.

    Returns:
    tuple: Flood class and flood layer images.
    """
    
    #  default values
    if zvv_thd is None:
        zvv_thd = -3
    if zvh_thd is None:
        zvh_thd = -3
    if pow_thd is None:
        pow_thd = 75
    if elev_thd is None:
        elev_thd = 800
    if slp_thd is None:
        slp_thd = 10
    
    if under_estimate is None:
        under_estimate = False
    
    # JRC water mask
    jrc = ee.ImageCollection("JRC/GSW1_4/MonthlyHistory").filterDate('2016-01-01', '2022-01-01')
    jrcvalid = jrc.map(lambda x: x.gt(0)).sum()
    jrcwat = jrc.map(lambda x: x.eq(2)).sum().divide(jrcvalid).multiply(100)
    jrcmask = jrcvalid.gt(0)
    ow = jrcwat.gte(ee.Image(pow_thd))

    # Elevation and slope masking using FABDEM
    elevation = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM").mosaic().setDefaultProjection('EPSG:3857', None, 30).clip(aoi)
    slope = ee.Terrain.slope(elevation).clip(aoi)

    # Classify floods
    vvflag = z.select('VV').lte(ee.Image(zvv_thd))
    vhflag = z.select('VH').lte(ee.Image(zvh_thd))
    flood_class = ee.Image(0).add(vvflag).add(vhflag.multiply(2)).where(ow.eq(1), 4).rename('flood_class')
    flood_class = flood_class.where(elevation.gt(elev_thd).multiply(ow.neq(1)), 0).where(slope.gt(slp_thd).multiply(ow.neq(1)), 0)

    # Combine flood classes into a single layer
    

    if under_estimate==True:
        # lowest probability vv+vh 
        flood_layer = flood_class.where(flood_class.eq(3), 1).where(flood_class.neq(3), 2)
        flood_layer = flood_layer.selfMask().rename('label')
    else:
        # highest probability combining all vv vh vv+vh flooded classes
        flood_layer = flood_class.where(flood_class.eq(1), 1).where(flood_class.eq(2), 1).where(flood_class.eq(3), 1).where(flood_class.eq(4), 2)
        flood_layer = flood_layer.where(jrcmask.eq(0), 2).where(flood_class.eq(0), 2).selfMask().rename('label')
    

    
    return flood_class.clip(aoi), flood_layer.clip(aoi)


# masking flood done

# Generate distance rasters
def distance_to_feature(feature_collection, crs, scale, aoi):
    
    '''	
    Generate distance rasters for the given feature collection.
    feature_collection (ee.FeatureCollection): The feature collection for which distance rasters are to be generated.
    crs (str): The CRS to reproject the distance raster.
    scale (int): The scale for reprojection.
    aoi (ee.Geometry): Area of Interest for clipping the distance raster.
    Returns:
    ee.Image: Distance raster for the feature collection.
    
    '''
    # Convert the FeatureCollection to a FeatureCollection
    feature_collection = ee.FeatureCollection(feature_collection).filterBounds(aoi)

    # Use the distance function to generate a distance raster
    distance_raster = feature_collection.distance()

    distance_raster = distance_raster.reproject(crs=crs, scale=scale)

    # Clip the raster to the AOI
    distance_raster = distance_raster.clip(aoi)

    return distance_raster


def prepare_datasets(aoi, projection='EPSG:4326', scale=30):
    """
    Prepare DEM, slope, and aspect datasets.

    Parameters:
    aoi (ee.Geometry): Area of Interest.
    projection (str): Projection to reproject the images. Default is 'EPSG:4326'.
    scale (int): Scale for reprojection. Default is 30.

    Returns:
    tuple: DEM, slope, aspect, and dtriver images reprojected to the specified projection.
    """
    dem_proj = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")\
        .filterBounds(aoi)\
        .mosaic()\
        .clip(aoi)\
        .setDefaultProjection('EPSG:3857', None, 30)

    slope_proj = ee.Terrain.slope(dem_proj)
    aspect_proj = ee.Terrain.aspect(dem_proj)

    dem = dem_proj.reproject(crs=projection, scale=scale)
    slope = slope_proj.reproject(crs=projection, scale=scale)
    aspect = aspect_proj.reproject(crs=projection, scale=scale)
    
    shoreline = ee.FeatureCollection('projects/sat-io/open-datasets/shoreline/mainlands')\
        .merge(ee.FeatureCollection('projects/sat-io/open-datasets/shoreline/big_islands'))\
        .filterBounds(aoi)
    
    rivers = ee.FeatureCollection("projects/sat-io/open-datasets/HydroAtlas/RiverAtlas_v10")\
        .filterBounds(aoi)

    # Combine rivers and shoreline into a single FeatureCollection
    rivers_and_shoreline = rivers.merge(shoreline)

    # Generate distance rasters for roads, and rivers+shoreline
    dtriver = distance_to_feature(rivers_and_shoreline, projection, scale, aoi)
    #rivers_and_shoreline_distance = distance_to_feature(rivers_and_shoreline, 30)

    return dem, slope, aspect, dtriver




def label_non_flooded(flood_binary_layer, aoi):
    """
    Label non-flooded pixels as 2 while keeping flooded pixels as 1.

    Parameters:
    flood_binary_layer (ee.Image): Binary flood layer with flooded pixels (1) and masked non-flooded pixels.
    aoi (ee.Geometry): Area of Interest.

    Returns:
    ee.Image: Layer with flooded pixels as 1 and non-flooded pixels as 2.
    """
    # Create a constant image with value 2 for the whole AOI
    non_flooded_layer = ee.Image.constant(2).clip(aoi)
    
    # Combine the flood binary layer with the non-flooded layer
    combined_layer = non_flooded_layer.where(flood_binary_layer.unmask().neq(0), flood_binary_layer)
    
    return combined_layer.rename('label')

def create_sample_feature_collection(flood_layer, flood_unmask, aoi, num_samples, class_band='label', scale=30):
    """
    Create a stratified sample feature collection from the flood layer.

    Parameters:
    flood_layer (ee.Image): The flood layer image.
    flood_unmask (bool): Flag to unmask non-flooded pixels. Default is False.
    aoi (ee.Geometry): Area of Interest.
    num_samples (int): Number of sample points.
    class_band (str): The band name containing the class labels. Default is 'label'.
    scale (int): Scale for sampling. Default is 20.

    Returns:
    ee.FeatureCollection: Sampled feature collection with updated labels.
    """
    
    if flood_unmask==True:
        flood_layer = label_non_flooded(flood_layer, aoi)
    
    sample = flood_layer.stratifiedSample(
        numPoints=num_samples,
        classBand=class_band,
        region=aoi,
        scale=scale,
        seed=5,
        tileScale=1.5,
        geometries=True
    )

    def update_feature(feature):
        value = feature.get(class_band)
        updated_value = ee.Algorithms.If(ee.Algorithms.IsEqual(value, ee.Number(2)), ee.Number(0), value)
        return feature.set(class_band, updated_value)
    label = sample.map(update_feature)
    return label

def prepare_s1_image(s1_post, additional_bands, aoi):
    """
    Prepare the image for classification by adding additional bands.

    Parameters:
    image (ee.Image): Base image to add additional bands.
    additional_bands (list): List of additional bands to add to the image.
    aoi (ee.Geometry): Area of Interest.

    Returns:
    ee.Image: Prepared image with added bands.
    """
    image = s1_post.mean().clip(aoi).toFloat()
    for band in additional_bands:
        image = image.addBands(band)
    return image

def create_training_and_validation_samples(image, label_fc, split, scale=30):
    """
    Create training and validation samples from the prepared image and label.

    Parameters:
    image (ee.Image): Prepared image with added bands.
    label (ee.FeatureCollection): Label feature collection.
    split (float): Ratio to split the samples into training and validation sets.
    scale (int): Scale for sampling. Default is 20.

    Returns:
    tuple: Training and validation feature collections.
    """
    sample_all = image.sampleRegions(
        collection=label_fc,
        properties=['label'],
        scale=scale
    ).randomColumn()

    training = sample_all.filter(ee.Filter.lt('random', split))
    validation = sample_all.filter(ee.Filter.gte('random', split))

    return training, validation

def train_classifier(training, bandNames):
    """
    Train a RandomForest classifier on the training samples.

    Parameters:
    training (ee.FeatureCollection): Training feature collection.
    band_names (list): List of band names used for classification.

    Returns:
    ee.Classifier: Trained RandomForest classifier.
    """
    return ee.Classifier.smileRandomForest(100).train(
        features=training,
        classProperty='label',
        inputProperties=bandNames
    )

def classify_image(image, classifier, bandNames):
    """
    Classify the image using the trained classifier.

    Parameters:
    image (ee.Image): Prepared image with added bands.
    classifier (ee.Classifier): Trained classifier.
    probability (bool): Flag to return probability values. Default is False.

    Returns:
    ee.Image: Classified image with binary flood mapping.
    OR
    ee.Image: Classified image with probability values.
    """
    classified_image = image.select(bandNames).classify(classifier)
    classified_image = classified_image.gt(0).selfMask().rename('flooded')
    
    return classified_image



def calculate_accuracy_metrics(training, validation, classifier):
    """
    Calculate accuracy metrics and plot ROC curve for the flood mapping.
    
    Parameters:
    training (ee.FeatureCollection): Training feature collection.
    validation (ee.FeatureCollection): Validation feature collection.
    classifier (ee.Classifier): Trained classifier.
    
    Returns:
    dict: A dictionary containing accuracy metrics and ROC-AUC score.
    """
    # Classify the training and validation data
    training_classified = training.classify(classifier)
    validation_classified = validation.classify(classifier)
    
    # Compute the error matrix for training and validation
    train_accuracy = training_classified.errorMatrix('label', 'classification')
    validation_accuracy = validation_classified.errorMatrix('label', 'classification')
    
    # Extract true and predicted labels from validation data
    accuracy = validation_accuracy.accuracy().getInfo()

    # Calculate overall precision and recall
    precision_list = validation_accuracy.consumersAccuracy().getInfo()
    recall_list = validation_accuracy.producersAccuracy().getInfo()
    precision = precision_list[0]
    recall = recall_list[0] # 0 for flooded class
    
    
    overall_precision = np.mean(precision_list)
    overall_recall = np.mean(recall_list)
 
    f1_overall = 2 * (overall_precision * overall_recall) / (overall_precision + overall_recall)
    f1 = validation_accuracy.fscore().getInfo()[0]
    # Calculate additional metrics using sklearn


    # create pd.DataFrame for results
    
    accuracy_dict = {
        'accuracy': accuracy,
        'precision': precision[0],
        'recall': recall[0],
        'f1_score': f1,
        'precision_mean': overall_precision,
        'recall_mean': overall_recall,
        'f1_mean': f1_overall
    }
    

    #covert accuracy_dict to pd.DataFrame
    results = pd.DataFrame.from_dict(accuracy_dict, orient='index')

    return results


def flood_mapping(aoi, s1_post, flood_layer, num_samples, split, city):
    # Prepare datasets
    dem, slope, aspect, dtriver = prepare_datasets(aoi)
    print('Done with preparing datasets...')
    # Create sample feature collection
    label_fc = create_sample_feature_collection(flood_layer,False, aoi, num_samples)
    print('Done with creating sample feature collection...')
    # Prepare image for classification
    additional_bands = [dem]#, slope, aspect, dtriver]
    image = prepare_s1_image(s1_post, additional_bands, aoi)

    # Create training and validation samples
    training, validation = create_training_and_validation_samples(image, label_fc, split)
    print('Done with creating training and validation samples...')
    # Train the classifier
    bandNames = image.bandNames().getInfo()
    classifier = train_classifier(training, bandNames)

    # Classify the image
    model_output = classify_image(image, classifier, bandNames)
    print('Done with classification...')
    # Calculate and print accuracy metrics
    results = calculate_accuracy_metrics(training, validation, classifier)
    results = results.round(2)
    results.to_csv(f'{city}_flood_mapping_accuracy.csv', index=True)
    print('Flood Mapping Accuracy Results: ',results)
    print('Done ...')
    return model_output





# Susceptibility

def apply_scale_factors(image):
    """
    Apply scaling factors to Landsat images.

    Parameters:
    image (ee.Image): Landsat image.

    Returns:
    ee.Image: Scaled Landsat image.
    """
    optical_bands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermal_bands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    return image.addBands(optical_bands, None, True).addBands(thermal_bands, None, True)

def prepare_landsat_images(aoi, endDate):
    
    """
    Prepare Landsat 8 and 9 images, combine them, and apply scaling factors.

    Parameters:
    aoi (ee.Geometry): Area of Interest.
    endDate (str): End date for filtering the images.

    Returns:
    ee.Image: Combined and processed Landsat image for the specified year.
    """
    year = ee.Date(endDate).get('year')

    l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')\
            .filterBounds(aoi)\
            .filterDate('2013-01-01', '2021-12-31')

    l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')\
            .filterBounds(aoi)


    landsat_combined = l8.merge(l9)

    landsat_filtered = landsat_combined.filter(ee.Filter.calendarRange(year, year, 'year'))\
                                       .filter(ee.Filter.lt('CLOUD_COVER', 20))\
                                       .map(apply_scale_factors)\
                                       .median()\
                                       .clip(aoi)
    
    return landsat_filtered


def prepare_datasets_for_susceptibility(aoi, landsat_filtered):
    
    """
    Prepare the necessary datasets for susceptibility analysis.

    Parameters:
    aoi (ee.Geometry): Area of Interest.
    landsat_filtered (ee.Image): Processed Landsat image.

    Returns:
    ee.Image: Combined image with all the necessary bands for susceptibility analysis.
    """
    
    dem_proj = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM")\
            .filterBounds(aoi)\
            .mosaic()\
            .clip(aoi).setDefaultProjection('EPSG:3857', None, 30).rename('elevation')

    slope_proj = ee.Terrain.slope(dem_proj)
    aspect_proj = ee.Terrain.aspect(dem_proj)
    
    dem = dem_proj.addBands(slope_proj).addBands(aspect_proj).reproject(crs='EPSG:4326', scale=30)

    ndvi = landsat_filtered.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI')
    ndwi = landsat_filtered.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI')
    ndbi = landsat_filtered.normalizedDifference(['SR_B6', 'SR_B5']).rename('NDBI')
    #ndsi = landsat_filtered.normalizedDifference(['SR_B2', 'SR_B5']).rename('NDSI')

    #rainfall = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')\
    #            .filterDate('2018-01-01', '2023-01-01')\
    #            .filterBounds(aoi)\
    #            .map(lambda image: image.gt(10).selfMask())\
    #            .map(lambda image: image.clip(aoi))\
    #            .sum().rename('rainfall')\
    #            .reproject(crs='EPSG:4326', scale=30)
    print('Done with preparing datasets for susceptibility analysis...')
    
    image_sus = landsat_filtered.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])\
                                .addBands([ndvi, ndwi, ndbi, dem])\
                                .clip(aoi)\
                                .setDefaultProjection('EPSG:4326')
    
    return image_sus


def train_susceptibility_model(image_sus, label, split, city):
    """
    Train a susceptibility model and calculate accuracy metrics.

    Parameters:
    image_sus (ee.Image): Combined image with all the necessary bands.
    label (ee.FeatureCollection): Feature collection with labels for training.
    split (float): Ratio to split the samples into training and validation sets.

    Returns:
    ee.Image: Image classified by the trained susceptibility model.
    """
    bands_sus = image_sus.bandNames().getInfo()
    
    print('Bands for susceptibility analysis:', bands_sus)
    
    sample_all_sus = image_sus.select(bands_sus).sampleRegions(
        collection=label,
        properties=['label'],
        scale=30,
        tileScale=1.5
    ).randomColumn()

    training_sus = sample_all_sus.filter(ee.Filter.lt('random', split))
    validation_sus = sample_all_sus.filter(ee.Filter.gte('random', split))

    classifier_sus = train_classifier(training_sus, bands_sus)
    
    classifier_prob = classifier_sus.setOutputMode('PROBABILITY')

    flood_prob = image_sus.classify(classifier_prob)

    results = calculate_accuracy_metrics(training_sus, validation_sus, classifier_sus)
    results = results.round(2)
    results.to_csv(f'{city}_flood_susceptibility_accuracy.csv', index=True)
    print('Flood Susceptibility Mapping Accuracy Results: ', results)
    
    
    return flood_prob

# Example usage for susceptibility analysis
def susceptibility_analysis(aoi, endDate, flood_binary, num_samples, split, city):
    # Prepare Landsat images
    landsat_filtered = prepare_landsat_images(aoi, endDate)

    # Prepare datasets
    image_sus = prepare_datasets_for_susceptibility(aoi, landsat_filtered)

    # Create sample feature collection
    label_new = create_sample_feature_collection(flood_binary.rename('label'), True, aoi, num_samples)
    # Train susceptibility model
    flood_prob = train_susceptibility_model(image_sus, label_new, split, city)

    return flood_prob


def quantile_based_categorization(susceptibility_layer, aoi):
    """
    Convert a continuous flood susceptibility layer into five categorical classes based on quantiles.

    Parameters:
    susceptibility_layer (ee.Image): Continuous flood susceptibility layer with values between 0 and 1.
    aoi (ee.Geometry): Area of Interest for calculating quantiles.

    Returns:
    ee.Image: Categorical flood susceptibility layer with values from 1 to 5.
    """
    # Calculate quantiles for the susceptibility layer
    quantiles = susceptibility_layer.reduceRegion(
        reducer=ee.Reducer.percentile([20, 40, 60, 80], ['p20', 'p40', 'p60', 'p80']),
        geometry=aoi,
        scale=30,
        maxPixels=1e9
    )
    #print('Quantiles:', quantiles.getInfo())
    # Extract quantile values
    q20 = quantiles.get('classification_p20')
    q40 = quantiles.get('classification_p40')
    q60 = quantiles.get('classification_p60')
    q80 = quantiles.get('classification_p80')

    # Apply quantile thresholds to create categorical classes
    very_low = susceptibility_layer.lte(ee.Number(q20)).multiply(1)
    low = susceptibility_layer.gt(ee.Number(q20)).And(susceptibility_layer.lte(ee.Number(q40))).multiply(2)
    moderate = susceptibility_layer.gt(ee.Number(q40)).And(susceptibility_layer.lte(ee.Number(q60))).multiply(3)
    high = susceptibility_layer.gt(ee.Number(q60)).And(susceptibility_layer.lte(ee.Number(q80))).multiply(4)
    very_high = susceptibility_layer.gt(ee.Number(q80)).multiply(5)

    # Combine all classes into a single layer
    categorical_layer = very_low.add(low).add(moderate).add(high).add(very_high).rename('susceptibility')

    return categorical_layer



# Exposure Analysis

def calculate_flood_exposure(flood_binary, aoi, export=False, city=None):
    """
    Calculate the population exposed to the flood based on the flood layer.

    Parameters:
    flood_layer (ee.Image): Flood layer with values 1 for flooded and 2 for non-flooded.
    population_dataset (str): Path to the population dataset.
    aoi (ee.Geometry): Area of Interest.

    Returns:
    ee.Number: Total population exposed to the flood.
    ee.Number: Total population in the study area.
    """
    
    if city is None:
        city = 'city'
    # Load population dataset
    flood_binary = flood_binary.reproject(crs='EPSG:3395', scale=10)
    
    population = ee.ImageCollection('WorldPop/GP/100m/pop')\
                    .filter(ee.Filter.eq('year', 2020))\
                    .mosaic()\
                    .clip(aoi).reproject(crs='EPSG:3395', scale=92.77)\
                    .rename('population')
                    
    # Mask non-flooded areas
    flood_exposure = population.updateMask(flood_binary.eq(1))
    
    # Calculate total exposed population
    total_exposed_population = flood_exposure.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=aoi,
        scale=92.7,
        maxPixels=1e13
    ).get('population')
    
     # Calculate total population in the study area
    total_population = population.reduceRegion(
        reducer=ee.Reducer.sum(),
        geometry=aoi,
        scale=92.7,
        maxPixels=1e13
    ).get('population')
    
    
    exposed_population = total_exposed_population.getInfo()
    total_pop = total_population.getInfo()
    non_exposed_population = total_pop - exposed_population
    
    labels = ['Exposed Population', 'Non-exposed Population']
    sizes = [exposed_population, non_exposed_population]
    colors = ['#fc4d4c', '#d8d8d8']  # Tomato color for exposed population, light green for non-exposed population
    
    plt.figure(figsize=(6, 6))
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
    plt.axis('equal')
    plt.title('Population Exposed to Flood Hazard')
    if export==True:
        plt.savefig(f'{city}_exposed_population_flood_hazard.png', dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    plt.show()
    
    
    return total_exposed_population, total_population





def calculate_exposure_df(susceptibility_layer, aoi, crs='EPSG:3395', flood_map=False, export=False, city=None):
    """
    Calculate exposure for population, nighttime light, and land cover for each susceptibility level or flood map.

    Parameters:
    susceptibility_layer (ee.Image): Flood susceptibility layer or flood map.
    population_dataset (str): Path to the population dataset.
    nightlight_dataset (str): Path to the nighttime light dataset.
    landcover_dataset (str): Path to the land cover dataset.
    aoi (ee.Geometry): Area of Interest.
    flood_map (bool): Flag to indicate if susceptibility categories or flood map should be used.

    Returns:
    pd.DataFrame: Dataframe with exposure information.
    """
    if city is None:
        city = 'city'
    
    # if aoi is ee.Feature, convert it to ee.Geometry
    if isinstance(aoi, ee.Feature):
        aoi = aoi.geometry()
    
    aoi = aoi.simplify(maxError=100)
    
        # Define remap function for landcover
    def remapper(image):
        return image.remap([1, 2, 4, 5, 7, 8, 9, 10, 11], [1, 2, 3, 4, 5, 6, 7, 8, 9])
    
    #susceptibility_layer = susceptibility_layer.reproject(crs=crs, scale=30)
    # Load datasets
    population = ee.ImageCollection('WorldPop/GP/100m/pop')\
                    .filter(ee.Filter.eq('year', 2020))\
                    .mosaic()\
                    .clip(aoi).rename('b1')\
                    #.reproject(crs=crs, scale=100)
                    
    nightlight = ee.Image('projects/sat-io/open-datasets/npp-viirs-ntl/LongNTL_2022').clip(aoi).rename('b1')#.reproject(crs=crs, scale=500)
    landcover = ee.ImageCollection('projects/sat-io/open-datasets/landcover/ESRI_Global-LULC_10m_TS')\
                    .filterDate('2022-01-01', '2022-12-31')\
                    .map(remapper)\
                    .mosaic().clip(aoi).rename('b1')\
                    #.reproject(crs=crs, scale=500)
    

    # Initialize results dictionary with descriptive land cover names
    results = {
        'Susceptibility Level': [],
        'Category': [],
        'Exposed Population': [],
        'Exposed Nighttime Light': [],
        'lulc_water': [],
        'lulc_trees': [],
        'lulc_flooded_vegetation': [],
        'lulc_crops': [],
        'lulc_built_area': [],
        'lulc_bare_ground': [],
        'lulc_snow_ice': [],
        'lulc_clouds': [],
        'lulc_rangeland': []
        
    }
    
    if flood_map:
        susceptibility_levels = [1]
        susceptibility_layer = susceptibility_layer.gt(0).selfMask()
        category_names = ['Flooded']
    else:
        # Define quantile-based categories
        susceptibility_levels = range(1, 6)
        category_names = ['Very Low', 'Low', 'Moderate', 'High', 'Very High']
    
    # Calculate exposure for each susceptibility level
    for level, category in zip(susceptibility_levels, category_names):
        # Mask areas that do not match the current susceptibility level
        level_mask = susceptibility_layer.eq(level)
        level_population = population.updateMask(level_mask)
        level_nightlight = nightlight.updateMask(level_mask)
        level_landcover = landcover.updateMask(level_mask)

        # Calculate total exposed population
        exposed_population = level_population.reduceRegion(
            reducer=ee.Reducer.sum(),
            geometry=aoi,
            scale=100,
            maxPixels=1e13
        ).get('b1').getInfo()

        # Calculate total exposed nighttime light
        exposed_nightlight = level_nightlight.reduceRegion(
            reducer=ee.Reducer.sum(),
            geometry=aoi,
            scale=500,
            maxPixels=1e13
        ).get('b1').getInfo()

        # Calculate land cover area for each class
        landcover_areas = level_landcover.reduceRegion(
            reducer=ee.Reducer.frequencyHistogram(),
            geometry=aoi,
            scale=100,
            maxPixels=1e13
        ).get('b1').getInfo()
        
        # Calculate area for each land cover class in square meters
        pixel_area = 100 * 100  # 100m resolution, so each pixel is 10000 square meters
        landcover_areas_m2 = {str(k): v * pixel_area for k, v in landcover_areas.items()}

        # Append results to the dictionary
        results['Susceptibility Level'].append(level)
        results['Category'].append(category)
        results['Exposed Population'].append(exposed_population)
        results['Exposed Nighttime Light'].append(exposed_nightlight)
        landcover_classes = {
            'lulc_water': 1,
            'lulc_trees': 2,
            'lulc_flooded_vegetation': 3,
            'lulc_crops': 4,
            'lulc_built_area': 5,
            'lulc_bare_ground': 6,
            'lulc_snow_ice': 7,
            'lulc_clouds': 8,
            'lulc_rangeland': 9
        }
        
        for lulc_name, lulc_code in landcover_classes.items():
            results[lulc_name].append(landcover_areas_m2.get(str(lulc_code), 0))
        
    # Convert results to dataframe
    exposure_df = pd.DataFrame(results)
    
    
    # Convert population values to integers
    exposure_df['Exposed Population'] = exposure_df['Exposed Population'].astype(int)

    # Convert LULC areas from m2 to km2
    lulc_columns = [
        'lulc_water', 'lulc_trees', 'lulc_flooded_vegetation', 'lulc_crops',
        'lulc_built_area', 'lulc_bare_ground', 'lulc_snow_ice', 'lulc_clouds', 'lulc_rangeland'
    ]

    for column in lulc_columns:
        exposure_df[column] = exposure_df[column] / 10000 # Convert m2 to hectare

    # Calculate the total values for population, nighttime light, and LULC areas
    total_population = exposure_df['Exposed Population'].sum()
    total_ntl = exposure_df['Exposed Nighttime Light'].sum()
    total_lulc = exposure_df[lulc_columns].sum()

    # Add percentage columns
    exposure_df['population_p'] = (exposure_df['Exposed Population'] / total_population) * 100
    exposure_df['ntl_p'] = (exposure_df['Exposed Nighttime Light'] / total_ntl) * 100

    for column in lulc_columns:
        percentage_column = column + '_p'
        exposure_df[percentage_column] = (exposure_df[column] / total_lulc[column]) * 100

    # Round the percentage columns to 2 decimal places
    percentage_columns = ['population_p', 'ntl_p'] + [col + '_p' for col in lulc_columns]
    exposure_df[percentage_columns] = exposure_df[percentage_columns].round(2)

    # round lulc columns to 2 decimal places
    exposure_df[lulc_columns] = exposure_df[lulc_columns].round(2)
    
    # round nighttime light column to 2 decimal places
    exposure_df['Exposed Nighttime Light'] = exposure_df['Exposed Nighttime Light'].round(2)

    # Fill NaN values with 0
    exposure_df.fillna(0, inplace=True)
    
    if export==True:
        exposure_df.to_csv(f'{city}_exposure_df.csv', index=False)
    
    return exposure_df




def visualize_exposure(df, export=False):
    """
    Generate visualizations for exposure analysis:
    1. A donut chart for population exposure to susceptibility per class.
    2. A donut chart for economic activity (nighttime light) exposure per susceptibility class.
    3. Five combined subplot charts for land cover area exposure for each flood susceptibility class.

    Parameters:
    df (pd.DataFrame): Dataframe containing exposure information.
    """
    # Colors for each susceptibility class
    colors = ['#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61', '#d7191c']
    labels = df['Category'].tolist()
    
    # Apply the ggplot theme
    sns.set_theme(style="whitegrid")
    sns.set_context("talk")

    
    # Function to create donut charts
    def create_donut_chart(data, title, labels, name, export):
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(aspect="equal"))
        wedges, texts, autotexts = ax.pie(data, labels=None, autopct='', startangle=140,
                                          colors=colors, pctdistance=0.85,
                                          wedgeprops={'edgecolor': 'black'})
        
        for i, a in enumerate(autotexts):
            percentage = f"{data[i]/sum(data)*100:.1f}%"
            label = labels[i]
            autotexts[i].set_text(f"{percentage}\n{label}")
            autotexts[i].set_fontsize(12)

        centre_circle = plt.Circle((0,0),0.70,fc='white')
        fig.gca().add_artist(centre_circle)
        plt.title(title, fontsize=16, fontweight='bold')
        plt.axis('equal')
        plt.tight_layout()
            
        if export==True:
            plt.savefig(f'{city}_{name}_exposure_sus.png', dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    
        plt.show()

    # Donut chart for population exposure
    create_donut_chart(df['Exposed Population'], 'Population Exposure to Flood Susceptibility', labels, 'pop', export)

    # Donut chart for economic activity (nighttime light) exposure
    create_donut_chart(df['Exposed Nighttime Light'], 'Economic Activity Exposure to Flood Susceptibility', labels, 'ntl', export)

    # Land cover columns
    landcover_columns = [
        'lulc_water', 'lulc_trees', 'lulc_flooded_vegetation', 'lulc_crops',
        'lulc_built_area', 'lulc_bare_ground', 'lulc_snow_ice', 'lulc_clouds', 'lulc_rangeland'
    ]

    # Create subplots for land cover area exposure
    fig, axes = plt.subplots(3, 2, figsize=(15, 20))
    axes = axes.flatten()

    for i, (category, row) in enumerate(df.iterrows()):
        landcover_data = row[landcover_columns] / 1e6  # Convert from m^2 to km^2
        landcover_data_sorted = landcover_data.sort_values(ascending=False)
        ax = axes[i]
        sns.barplot(x=landcover_data_sorted.index, y=landcover_data_sorted.values, palette="YlOrRd_r", ax=ax)

        # Adjust the positioning of the text labels on the bars
        for bar in ax.patches:
            yval = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, yval + 0.02 * landcover_data_sorted.max(), f'{yval:.2f}', ha='center', va='bottom', fontsize=10)

        ax.set_title(f'Land Cover Area Exposed to {row["Category"]} Susceptibility', fontsize=16, fontweight='bold')
        ax.set_ylabel('Area (sq km)', fontsize=14)
        ax.set_xlabel('Land Cover Type', fontsize=14)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12)
        ax.set_ylim(0, landcover_data_sorted.max() * 1.1)
        ax.set_facecolor('white')

    # Remove the unused subplot
    fig.delaxes(axes[-1])

    plt.tight_layout()
    
    if export==True:
        plt.savefig(f'{city}_lulc_exposure_sus.png', dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    
    plt.show()




# export to asset

def export_image_to_asset(image, description, asset_id, region, scale):
    """
    Export an ee.Image to an Earth Engine asset and monitor the upload progress.

    Parameters:
    image (ee.Image): The image to export.
    description (str): A description for the export task.
    asset_id (str): The destination asset ID.
    region (ee.Geometry): The region to export.
    scale (int): The scale (in meters) to export the image.

    Returns:
    None
    """
    # Create the export task
    export_task = ee.batch.Export.image.toAsset(
        image=image,
        description=description,
        assetId=asset_id,
        region=region,
        scale=scale,
        maxPixels=1e13
    )

    # Start the export task
    export_task.start()

    # Function to monitor the export task
    def monitor_task(task):
        while task.active():
            print('Task status:', task.status())
            time.sleep(10)
        print('Task finished:', task.status())

    # Monitor the export task
    monitor_task(export_task)


# export to gdrive
def export_layers(aoi, city, flood_layer, flood_class, flood_mapped, susceptibility_layer, susceptibility_category_layer,
                  export_flood_layer=True, export_flood_class=True, export_flood_mapped=True,
                  export_susceptibility_layer=True, export_susceptibility_category_layer=True):
    """
    Export specified layers to Google Drive.

    Parameters:
    aoi (ee.Geometry): Area of Interest.
    flood_binary (ee.Image): Flood binary layer.
    flood_class (ee.Image): Flood class layer.
    flood_mapped (ee.Image): Flood mapped layer.
    susceptibility_layer (ee.Image): Susceptibility layer.
    susceptibility_category_layer (ee.Image): Susceptibility category layer.
    export_flood_binary (bool): Flag to export flood binary layer.
    export_flood_class (bool): Flag to export flood class layer.
    export_flood_mapped (bool): Flag to export flood mapped layer.
    export_susceptibility_layer (bool): Flag to export susceptibility layer.
    export_susceptibility_category_layer (bool): Flag to export susceptibility category layer.
    """

    tasks = []
    aoi = aoi.geometry()
    if export_flood_layer:
        flood_binary_task = ee.batch.Export.image.toDrive(
            image=flood_layer,
            description=f'{city}_flood_mask_layer',
            folder='FMSE',
            scale=10,
            region=aoi,
            maxPixels=1e13
        )
        flood_binary_task.start()
        tasks.append(flood_binary_task)

    if export_flood_class:
        flood_class_task = ee.batch.Export.image.toDrive(
            image=flood_class,
            description=f'{city}_flood_class_layer',
            folder='FMSE',
            scale=10,
            region=aoi,
            maxPixels=1e13
        )
        flood_class_task.start()
        tasks.append(flood_class_task)

    if export_flood_mapped:
        flood_mapped_task = ee.batch.Export.image.toDrive(
            image=flood_mapped,
            description=f'{city}_flood_mapped_layer',
            folder='FMSE',
            scale=10,
            region=aoi,
            maxPixels=1e13
        )
        flood_mapped_task.start()
        tasks.append(flood_mapped_task)

    if export_susceptibility_layer:
        susceptibility_task = ee.batch.Export.image.toDrive(
            image=susceptibility_layer,
            description=f'{city}_flood_susceptibility_layer',
            folder='FMSE',
            scale=30,
            region=aoi,
            maxPixels=1e13
        )
        susceptibility_task.start()
        tasks.append(susceptibility_task)

    if export_susceptibility_category_layer:
        susceptibility_category_task = ee.batch.Export.image.toDrive(
            image=susceptibility_category_layer,
            description=f'{city}_flood_susceptibility_category_layer',
            folder='FMSE',
            scale=30,
            region=aoi,
            maxPixels=1e13
        )
        susceptibility_category_task.start()
        tasks.append(susceptibility_category_task)

    def monitor_tasks(tasks):
        while any([task.status()['state'] in ['READY', 'RUNNING'] for task in tasks]):
            for task in tasks:
                status = task.status()
                description = status['description']
                state = status['state']
                print(f'Task {description} is {state}')
            time.sleep(30)  # Check every 30 seconds

    # Monitor the export tasks
    monitor_tasks(tasks)
