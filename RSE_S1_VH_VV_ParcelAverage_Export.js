var parcel_B = ee.FeatureCollection("users/matthieutaymans/Bavaria/nuts_classif/B_classif_allB_simpVectors_22-03_2_co"),
    parcel_MVBRG = ee.FeatureCollection("users/matthieutaymans/MVBRG/nuts_classif/MVBRG_classif_allB_simpVectors");
    
    
// var parcels = ee.FeatureCollection("users/matthieutaymans/Bavaria/fields_manual/B_fields_final");
// var parcels = ee.FeatureCollection("users/matthieutaymans/MVBRG/fields_manual/MV_fields_final");


// var parcels=parcel_MVBRG
var parcels=parcel_B

var s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT').filterDate('2018-01-01', '2018-08-01').
  filterMetadata('instrumentMode', 'equals', 'IW').filterBounds(parcels)
  
// var parcels = parcels.map(function(f) { return f.set({'id': f.id(), 'area': f.area(), 'perimeter': f.perimeter()}) })
// var parcels = parcels.filterMetadata('area', 'greater_than', 10000.0)
// print(parcels)

var fset = s1.map(function(img) {
  var parcelsInside = parcels.filterBounds(img.geometry())
  var imgExtract = parcelsInside.map(function(f) {
    return f.set(img.reduceRegion(ee.Reducer.mean(), f.buffer(-10).geometry(), 10).set('tstamp', img.date()))
  })
  return imgExtract
})

print(fset.flatten().limit(10))


print(fset.flatten().select(ee.List(['VH','VV']), null, false).limit(10))
// var fset=fset.limit(1)
// Export.table.toDrive(fset.flatten(),
//   "S1_export_B_20190410_full",
//   "/EarthEngine/",
//   'ExportS1_Bavaria_Classif_20190405_full')
  
Export.table.toDrive(fset.flatten().select(ee.List(['VH','VV']), null, false),
  "S1_export_B_20190419",
  "/EarthEngine/",
  'ExportS1_Bavaria_Classif_20190419_full')

// Export.table.toDrive(fset.flatten().select(ee.List(['VH','VV']), null, false),
//   "S1_export_MVBRG_20190405",
//   "/EarthEngine/",
//   'ExportS1_MVBRG_Classif_20190405_full')
