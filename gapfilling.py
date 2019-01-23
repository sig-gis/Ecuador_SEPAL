
import ee
ee.Initialize()

def calculateGaps(collection):
	
	# make a gap 1 and data 0
	def unmaskNoData(img):  
		im = img.select("red").gt(9999).rename(["gaps"])
		gaps = ee.Image(im.unmask(1).set('system:time_start', img.get('system:time_start'))).float()
		return gaps
  
	# iterate over the collection and count the number of consequetive gaps
	def gapsCumulative(img,mylist):
		previous = ee.Image(ee.List(mylist).get(-1));
		gaps = img.add(previous)
		gaps = gaps.where(gaps.eq(previous),0).set('system:time_start', img.get('system:time_start'))
		return ee.List(mylist).add(gaps)
	
	gapCollection = ee.ImageCollection(collection.map(unmaskNoData))

	# count the number of gaps
	first = ee.List([ee.Image(gapCollection.first())]);
	previous = ee.Image(ee.List(first).get(-1));
	gaps = ee.ImageCollection(ee.List(gapCollection.iterate(gapsCumulative, first)));

	return ee.ImageCollection(gaps);
  

def fillCollectionYear(gaps,inCollection,year):
	
	stDay = year*365+8
	endDay = year*365-8
	# appy gap filling using the year before
	def fillGaps(img):
		date = ee.Date(img.get("system:time_start"))
		startDate = date.advance(-8,"day")
		endDate = date.advance(8,"day")
		thisyear = ee.Image(inCollection.filterDate(startDate,endDate).first())  
		startDate = date.advance(-stDay,"day")
		endDate = date.advance(-endDay,"day")
		previousYear = ee.ImageCollection(ic.filterDate(startDate,endDate))
		
		def returnMask(img,previousYear):
			mask = img.gt(2)
			previousYear =previousYear.updateMask(mask)
			img = thisyear.unmask(previousYear)
			return img
		
		img = ee.Algorithms.If(previousYear.size().gt(0),returnMask(img,ee.Image(previousYear.first())),thisyear)
		return img
	collection = gaps.map(fillGaps)
  
	return ee.ImageCollection(collection)



ic = ee.ImageCollection("projects/Sacha/PreprocessedData/L8_Biweekly_V4_Mosaiced")

myGaps = calculateGaps(ic)
myCollection = ee.ImageCollection(fillCollectionYear(myGaps,ic,1))

myGaps = calculateGaps(myCollection)
myCollection = fillCollectionYear(myGaps,myCollection,2)

myGaps = calculateGaps(myCollection)
myCollection = fillCollectionYear(myGaps,myCollection,3)

myList = myCollection.toList(500)

im = ee.Image(myList.get(150))

studyArea =  ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
studyArea = studyArea.geometry().bounds()

geom = studyArea.getInfo()

task_ordered= ee.batch.Export.image.toAsset(image=im, 
							  description = "image", 
							  assetId= "projects/Sacha/PreprocessedData/test_gapfilled/test1",
							  region=geom['coordinates'], 
							  maxPixels=1e13,
							  crs="epsg:4326",
							  scale=30)
	
task_ordered.start()




