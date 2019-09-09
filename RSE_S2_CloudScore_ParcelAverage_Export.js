var nuts = ee.FeatureCollection("users/matthieutaymans/NutsRapeseed"),
    parcel_MVBRG = ee.FeatureCollection("users/matthieutaymans/MVBRG/nuts_classif/MVBRG_classif_allB_simpVectors"),
    parcel_B = ee.FeatureCollection("users/matthieutaymans/Bavaria/nuts_classif/B_classif_allB_simpVectors_22-03_2_co"),
    parcel_B_manual = ee.FeatureCollection("users/matthieutaymans/Bavaria/fields_manual/B_fields_final"),
    parcel_MV_manual = ee.FeatureCollection("users/matthieutaymans/MVBRG/fields_manual/MV_fields_final");
    
///////////////////////////////////////////
//   INPUTS
///////////////////////////////////////////

// 1. Date

// var start_date='2018-01-01';   // to be defined 

var start_date='2018-01-01';
var end_date='2018-08-01';   // to be defined

//  2. Cloud score

var max_cloud_percent= 80;    // remove cloudy images but i would suggest using any acquisition
var sort_ascending = false; 

//Parameters for sentinel cloud score
var cloudThresh = 20;//Ranges from 1-100.Lower value will mask more pixels out. Generally 10-30 works well with 20 being used most commonly 
var cloudHeights = ee.List.sequence(200,10000,500);//Height of clouds to use to project cloud shadows
var irSumThresh = 0.35;//Sum of IR bands to include as shadows within TDOM and the shadow shift method (lower number masks out less)
var dilatePixels = 5; //Pixels to dilate around clouds
var contractPixels = 1;//Pixels to reduce cloud mask and dark shadows by to reduce inclusion of single-pixel comission errors

//buffer around the clouds
var cloud_buffer = 50;

// 2. Parcels

// BAVARIA
// var parcel = parcel_B 
var parcel = parcel_B_manual


// MV
// var parcel = parcel_MVBRG
// var parcel = parcel_MV_manual


 var parcel = parcel.map(function(f) { return f.buffer(-10)})
// var parcel = parcel.map(function(f) { return f.set({'bufferedarea': f.area()}) })
// print(parcel.limit(10))

// Map.addLayer(parcel,  {color: '0000FF'},'parcel');

// var parcel =parcel.limit(10) // For dvpt, limit the number of parcel

// 3. Filter bound on NUTS

// BAVARIA
var extent=['DE2']; 
 // MV
// var extent=['DE4','DE8']; 

var aoi = nuts.filter(ee.Filter.inList('NUTS_ID',extent));


// 4. Output names, dir and tasks


// Choose locality
var loc='B';//or MV B
var delineationMethod='manual';//or manual'

// Output file name
var outputInit = 'StudyArea_20190401'
var outputtask = outputInit.concat('_',delineationMethod,'Fields_',loc);
var outputfolder = '/EarthEngine/';
var outputfilename = outputtask


///////////////////////////////////////////
// A / FUNCTIONS
///////////////////////////////////////////

// 1 . Cloud score

var rescale = function(img, exp, thresholds) {
    return img.expression(exp, {img: img})
        .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
  };

////////////////////////////////////////
// Cloud masking algorithm for Sentinel2
//Built on ideas from Landsat cloudScore algorithm
//Currently in beta and may need tweaking for individual study areas
function sentinelCloudScore(img) {
  

  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1);
  
  // Clouds are reasonably bright in the blue and cirrus bands.
  score = score.min(rescale(img, 'img.blue', [0.1, 0.5]));
  score = score.min(rescale(img, 'img.cb', [0.1, 0.3]));
  score = score.min(rescale(img, 'img.cb + img.cirrus', [0.15, 0.2]));
  
  // Clouds are reasonably bright in all visible bands.
  score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]));

  
  //Clouds are moist
  var ndmi = img.normalizedDifference(['nir','swir1']);
  score=score.min(rescale(ndmi, 'img', [-0.1, 0.1]));
  
  // However, clouds are not snow.
  var ndsi = img.normalizedDifference(['green', 'swir1']);
  score=score.min(rescale(ndsi, 'img', [0.8, 0.6]));
  
  score = score.multiply(100).byte();
 
  return img.addBands(score.rename('cloudScore'));
}


function projectShadows(cloudMask,image,cloudHeights){
  var meanAzimuth = image.get('MEAN_SOLAR_AZIMUTH_ANGLE');
  var meanZenith = image.get('MEAN_SOLAR_ZENITH_ANGLE');

  
  //Find dark pixels
  var darkPixels = image.select(['nir','swir1','swir2']).reduce(ee.Reducer.sum()).lt(irSumThresh)
    .focal_min(contractPixels).focal_max(dilatePixels)
  ;//.gte(1);
    
  
  //Get scale of image
  var nominalScale = cloudMask.projection().nominalScale();

  //Find where cloud shadows should be based on solar geometry
  //Convert to radians
  var azR =ee.Number(meanAzimuth).multiply(Math.PI).divide(180.0).add(ee.Number(0.5).multiply(Math.PI ));
  var zenR  =ee.Number(0.5).multiply(Math.PI ).subtract(ee.Number(meanZenith).multiply(Math.PI).divide(180.0));
  
  //Find the shadows
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    
    var shadowCastedDistance = zenR.tan().multiply(cloudHeight);//Distance shadow is cast
    var x = azR.cos().multiply(shadowCastedDistance).divide(nominalScale).round();//X distance of shadow
    var y = azR.sin().multiply(shadowCastedDistance).divide(nominalScale).round();//Y distance of shadow
    return cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y));
  });

  var shadowMask = ee.ImageCollection.fromImages(shadows).max();
 
  //Create shadow mask
  shadowMask = shadowMask.and(cloudMask.not());
  shadowMask = shadowMask.and(darkPixels);
  
  var cloudShadowMask = shadowMask.or(cloudMask);
  
  image = image.updateMask(cloudShadowMask.not()).addBands(shadowMask.rename(['cloudShadowMask']));
  return image;
}

//Function for wrapping the entire process to be applied across collection
function wrapIt(img){
  img = sentinelCloudScore(img);
  var cloudMask = img.select(['cloudScore']).gt(cloudThresh)
    .focal_min(contractPixels).focal_max(dilatePixels)

  img = projectShadows(cloudMask,img,cloudHeights);
  
  return img;
}


//  2 . Add indices

// ADD A NDVI BAND
function addNdyi(img) {
  var nd = img.normalizedDifference(['green', 'blue']);
  return img.addBands(nd.float().rename('NDYI'));
}

// ADD A NDYI BAND
function addNdvi(img) {
  var nd = img.normalizedDifference(['nir', 'red']);
  return img.addBands(nd.float().rename('NDVI'));
}

// ADD A BSI BAND
function addBSI(img) {
  var bsi = img.expression('((swir1+red)-(nir+blue))/(((swir1+red)+(nir+blue))*1.0)', {
      'swir1': img.select('swir1'),
      'red': img.select('red'),
      'nir': img.select('nir'),
      'blue': img.select('blue')
}
);
    return img.addBands(bsi.float().rename('BSI'));
}


//Function to add date
var addDate = function(image){
   var date = image.date().format('YYYYMMdd');
   return image.set('date', ee.Number.parse(date));
};



/* Sentinel-2 produces multiple images, resulting sometimes 4x more images than the actual size. 
  This is bad for any statistical analysis.  This function mosaics images by time. */
var mosaicByTime = function(images) {
  var TIME_FIELD = 'date'; //system:index

  var distinct = images.distinct([TIME_FIELD]);

  var filter = ee.Filter.equals({ leftField: TIME_FIELD, rightField: TIME_FIELD });
  var join = ee.Join.saveAll('matches');
  var results = join.apply(distinct, images, filter);

  // mosaic
  results = results.map(function(i) {
    var mosaic = ee.ImageCollection.fromImages(i.get('matches')).sort('date').mosaic(); //system:index
    
    return mosaic.copyProperties(i).set(TIME_FIELD, i.get(TIME_FIELD));
  });
  
  return ee.ImageCollection(results);
};

//create mosaics

///////////////////////////////////////////
// B  S2 LOAD AND CLOUD SCREENING
///////////////////////////////////////////
 var s2 = ee.ImageCollection('COPERNICUS/S2').filterDate(start_date,end_date)  //
                 .filterMetadata('CLOUDY_PIXEL_PERCENTAGE',"less_than",max_cloud_percent)
                 .filterBounds(aoi)
                 .sort('system:time_start');

var s2 = s2.map(addDate);
// print(s2)

var collectionS2 = mosaicByTime(s2)
// print('Mosaic by Date',collectionS2)


var collectionS2_noCloud = collectionS2.map(function(img){
    var t = img.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000);//Rescale to 0-1
    t = t.addBands(img.select(['QA60']));
    var out = t.copyProperties(img).copyProperties(img,['system:time_start']);
       return out;
          })
    .select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']
    ,['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2'])
     .map(wrapIt);

// print(collectionS2_noCloud)



///////////////////////////////////////////
// C ADD INDICES
///////////////////////////////////////////

// ADD NDVI BANDS
var collectionS2_noCloud = collectionS2_noCloud.map(addNdvi);

// ADD BSI BAND
var collectionS2_noCloud = collectionS2_noCloud.map(addNdyi);



///////////////////////////////////////////
// F / Export the parcel means  to a table
///////////////////////////////////////////

// FOR CLASSIFIED FIELDS
// var fc = collectionS2_noCloud.map(function(img) {
//   return img.reduceRegions({collection: parcel, reducer:ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true).combine(ee.Reducer.count(), '', true), scale: 10});
// })

// Export.table.toDrive({collection:ee.FeatureCollection(fc.flatten()).filter(ee.Filter.notNull(['NDVI_mean'])),
// description:outputtask,folder:outputfolder,selectors:['system:index','NDVI_mean','NDVI_stdDev','NDVI_count', 'NDYI_mean','NDYI_stdDev']}) //

// FOR VISUAL DELINEATED FIELDS
var fc = collectionS2_noCloud.map(function(img) {
  return img.reduceRegions({collection: parcel, reducer:ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true).combine(ee.Reducer.count(), '', true), scale: 10});
})
Export.table.toDrive({collection:ee.FeatureCollection(fc.flatten()).filter(ee.Filter.notNull(['NDVI_mean'])),
description:outputtask,folder:outputfolder,selectors:['system:index','field_id','NDVI_mean','NDVI_stdDev','NDVI_count', 'NDYI_mean','NDYI_stdDev']}) //

