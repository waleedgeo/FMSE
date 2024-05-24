# Part-1: Flood Mapping

import ee

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

def create_sample_feature_collection(flood_layer, flood_unmask, aoi, num_samples, class_band='label', scale=20):
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

def create_training_and_validation_samples(image, label_fc, split, scale=20):
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
    return ee.Classifier.smileRandomForest(115).train(
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
    Calculate and print accuracy metrics for training and validation samples.

    Parameters:
    training (ee.FeatureCollection): Training feature collection.
    validation (ee.FeatureCollection): Validation feature collection.
    classifier (ee.Classifier): Trained classifier.
    """
    def print_am(error_matrix, dataset_name):
        overall_accuracy = error_matrix.accuracy().getInfo()
        kappa = error_matrix.kappa().getInfo()
        producers_accuracy = error_matrix.producersAccuracy().getInfo()
        consumers_accuracy = error_matrix.consumersAccuracy().getInfo()
        f1_score = error_matrix.fscore().getInfo()[1]  # F1 score for the flood class (1)

        print(f"{dataset_name} Accuracy Metrics:")
        print(f"Overall Accuracy: {overall_accuracy}")
        print(f"Kappa: {kappa}")
        print(f"Producer's Accuracy: {producers_accuracy}")
        print(f"User's Accuracy: {consumers_accuracy}")
        print(f"F1 Score: {f1_score}")

    train_accuracy = training.classify(classifier).errorMatrix('label', 'classification')
    validation_accuracy = validation.classify(classifier).errorMatrix('label', 'classification')

    print_am(train_accuracy, "Training")
    print_am(validation_accuracy, "Validation")


def flood_mapping(aoi, s1_post, flood_layer, num_samples, split):
    # Prepare datasets
    dem, slope, aspect, dtriver = prepare_datasets(aoi)
    print('Done with preparing datasets...')
    # Create sample feature collection
    label_fc = create_sample_feature_collection(flood_layer,False, aoi, num_samples)
    print('Done with creating sample feature collection...')
    # Prepare image for classification
    additional_bands = [dem, slope, aspect, dtriver]
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
    calculate_accuracy_metrics(training, validation, classifier)
    print('Done with accuracy metrics...')
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


def train_susceptibility_model(image_sus, label, split):
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
    print('Training sus first: ',training_sus.first().getInfo())
    classifier_sus = ee.Classifier.smileRandomForest(115).train(
        features=training_sus,
        classProperty='label',
        inputProperties=bands_sus
    )

    classifier_prob = classifier_sus.setOutputMode('PROBABILITY')

    flood_prob = image_sus.classify(classifier_prob)

    def calculate_metrics(validation, classifier):
        validation_classified = validation.classify(classifier)
        validation_accuracy = validation_classified.errorMatrix('label', 'classification')

        f1_score = validation_accuracy.fscore().getInfo()[1] if len(validation_accuracy.fscore().getInfo()) > 1 else None
        producer_accuracy = validation_accuracy.producersAccuracy().getInfo()[1] if len(validation_accuracy.producersAccuracy().getInfo()) > 1 else None
        consumer_accuracy = validation_accuracy.consumersAccuracy().getInfo()[1] if len(validation_accuracy.consumersAccuracy().getInfo()) > 1 else None

        return f1_score, producer_accuracy, consumer_accuracy

    f1_score, producer_accuracy, consumer_accuracy = calculate_metrics(validation_sus, classifier_sus)

    print("Validation F1 Score:", f1_score)
    print("Validation Producer Accuracy:", producer_accuracy)
    print("Validation Consumer Accuracy:", consumer_accuracy)

    return flood_prob

# Example usage for susceptibility analysis
def susceptibility_analysis(aoi, endDate, flood_binary, num_samples, split):
    # Prepare Landsat images
    landsat_filtered = prepare_landsat_images(aoi, endDate)

    # Prepare datasets
    image_sus = prepare_datasets_for_susceptibility(aoi, landsat_filtered)

    # Create sample feature collection
    label_new = create_sample_feature_collection(flood_binary.rename('label'), True, aoi, num_samples)
    print('Size of label collection:', label_new.size().getInfo())

    # Train susceptibility model
    flood_prob = train_susceptibility_model(image_sus, label_new, split)

    return flood_prob




