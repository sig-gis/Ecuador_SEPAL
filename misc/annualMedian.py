import ee
from landsat import *
ee.Initialize()

func = functions()
p = [25,75]

geom =  ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete").geometry().bounds().getInfo()
collection = ee.ImageCollection("projects/Sacha/PreprocessedData/L8_Biweekly_V5_Mosaiced_2")
# remove bands of biweekly percentiles
newNames = ee.List(['p25_green','p25_nir','p25_blue','p25_red','p25_swir1','p25_swir2','p75_green','p75_nir','p75_blue','p75_red','p75_swir1','p75_swir2'])
badNames = ee.List(['green_p25','nir_p25','blue_p25','red_p25','swir1_p25','swir2_p25','green_p75','nir_p75','blue_p75','red_p75','swir1_p75','swir2_p75'])
badStdDev = ee.List(['count','blue_stdDev','red_stdDev','green_stdDev','nir_stdDev','swir1_stdDev','swir2_stdDev','ND_green_swir1','ND_nir_red_stdDev','ND_nir_swir2_stdDev','svvi'])

bn = collection.first().bandNames()
sn = bn.removeAll(newNames.cat(badStdDev))
collection = collection.select(sn)

exportLive = False
start = 0
sy = 2009
for i in range(2,3,1):

	year = ee.Date("2009-01-01")
	startDay = 0
	endDay = 364

	startDate = year.advance(startDay,'day').advance(i,'year') 
	endDate = year.advance(endDay,'day').advance(i,'year')

	tcollection = collection.filterDate(startDate,endDate)

	medianMosaic = func.medianMosaic(tcollection)
	medianPercentiles = func.medianPercentiles(tcollection,p)
	medianPercentiles = medianPercentiles.select(badNames,newNames)
	stdDevBands = func.addSTDdev(tcollection)

	median = medianMosaic.addBands(medianPercentiles).addBands(stdDevBands)
	median = median.set('system:time_start',startDate)
	nameYr = sy+i
	
	name = 'LS_AN_' + str(nameYr) +'000365'

	if exportLive == True:
		task_ordered= ee.batch.Export.image.toAsset(image=median,
									  description = name,
									  assetId= "projects/Sacha/PreprocessedData/L8_Annual_V5/" + name,
									  region=geom['coordinates'],
									  maxPixels=1e13,
									  crs="epsg:32717",
									  scale=30)
		print(name)	
		task_ordered.start()
	else:
		print(str(median.getInfo()['bands']).split("u'id"),name)