
import ee
ee.Initialize()

filez = ["LS_BW_200901402","LS_BW_200904205"] #,"LS_BW_200928029","LS_BW_200929430","LS_BW_200930832","LS_BW_200932233","LS_BW_200933634","LS_BW_200935036","LS_BW_200936401","LS_BW_201001302","LS_BW_201002704","LS_BW_201004105","LS_BW_201005506","LS_BW_201006908","LS_BW_201008309","LS_BW_201009711","LS_BW_201011112","LS_BW_201012513","LS_BW_201013915","LS_BW_201015316","LS_BW_201016718","LS_BW_201018119","LS_BW_201019520","LS_BW_201020922","LS_BW_201022323","LS_BW_201023725","LS_BW_201025126","LS_BW_201026527","LS_BW_201027929","LS_BW_201029330","LS_BW_201030732","LS_BW_201032133","LS_BW_201033534","LS_BW_201034936","LS_BW_201036301","LS_BW_201101202","LS_BW_201102603","LS_BW_201104005","LS_BW_201105406","LS_BW_201106808","LS_BW_201108209","LS_BW_201109610","LS_BW_201111012","LS_BW_201112413","LS_BW_201113815","LS_BW_201115216","LS_BW_201116617","LS_BW_201118019","LS_BW_201119420","LS_BW_201120822","LS_BW_201122223","LS_BW_201123624","LS_BW_201125026","LS_BW_201126427","LS_BW_201127829","LS_BW_201129230","LS_BW_201130631","LS_BW_201132033","LS_BW_201133434","LS_BW_201134836","LS_BW_201136201","LS_BW_201201102","LS_BW_201202503","LS_BW_201203905","LS_BW_201205306","LS_BW_201206708","LS_BW_201208109","LS_BW_201209510","LS_BW_201210912","LS_BW_201212313","LS_BW_201213715","LS_BW_201215116","LS_BW_201216517","LS_BW_201217919","LS_BW_201219320","LS_BW_201220722","LS_BW_201222123","LS_BW_201223524","LS_BW_201224926","LS_BW_201226327","LS_BW_201227729","LS_BW_201229130","LS_BW_201230531","LS_BW_201231933","LS_BW_201233334","LS_BW_201234736","LS_BW_201236100","LS_BW_201300902","LS_BW_201302303","LS_BW_201303705","LS_BW_201305106","LS_BW_201306507","LS_BW_201307909","LS_BW_201309310","LS_BW_201310712","LS_BW_201312113","LS_BW_201313514","LS_BW_201314916","LS_BW_201316317","LS_BW_201317719","LS_BW_201319120","LS_BW_201320521","LS_BW_201321923","LS_BW_201323324","LS_BW_201324726","LS_BW_201326127","LS_BW_201327528","LS_BW_201328930","LS_BW_201330331","LS_BW_201331733","LS_BW_201333134","LS_BW_201334535","LS_BW_201335900","LS_BW_201400802","LS_BW_201402203","LS_BW_201403604","LS_BW_201405006","LS_BW_201406407","LS_BW_201407809","LS_BW_201409210","LS_BW_201410611","LS_BW_201412013","LS_BW_201413414","LS_BW_201414816","LS_BW_201416217","LS_BW_201417618","LS_BW_201419020","LS_BW_201420421","LS_BW_201421823","LS_BW_201423224","LS_BW_201424625","LS_BW_201426027","LS_BW_201427428","LS_BW_201428830","LS_BW_201430231","LS_BW_201431632","LS_BW_201433034","LS_BW_201434435","LS_BW_201435800","LS_BW_201500702","LS_BW_201502103","LS_BW_201503504","LS_BW_201504906","LS_BW_201506307","LS_BW_201507709","LS_BW_201509110","LS_BW_201510511","LS_BW_201511913","LS_BW_201513314","LS_BW_201514716","LS_BW_201516117","LS_BW_201517518","LS_BW_201518920","LS_BW_201520321","LS_BW_201521723","LS_BW_201523124","LS_BW_201524525","LS_BW_201525927","LS_BW_201527328","LS_BW_201528730","LS_BW_201530131","LS_BW_201531532","LS_BW_201532934","LS_BW_201534335","LS_BW_201535700","LS_BW_201600601","LS_BW_201602003","LS_BW_201603404","LS_BW_201604806","LS_BW_201606207","LS_BW_201607608","LS_BW_201609010","LS_BW_201610411","LS_BW_201611813","LS_BW_201613214","LS_BW_201614615","LS_BW_201616017","LS_BW_201617418","LS_BW_201618820","LS_BW_201620221","LS_BW_201621622","LS_BW_201623024","LS_BW_201624425","LS_BW_201625827","LS_BW_201627228","LS_BW_201628629","LS_BW_201630031","LS_BW_201631432","LS_BW_201632834","LS_BW_201634235","LS_BW_201635600","LS_BW_201700401","LS_BW_201701803","LS_BW_201703204","LS_BW_201704605","LS_BW_201706007","LS_BW_201707408","LS_BW_201708810","LS_BW_201710211","LS_BW_201711612","LS_BW_201713014","LS_BW_201714415","LS_BW_201715817","LS_BW_201717218","LS_BW_201718619","LS_BW_201720021","LS_BW_201721422","LS_BW_201722824","LS_BW_201724225","LS_BW_201725626","LS_BW_201727028","LS_BW_201728429","LS_BW_201729831","LS_BW_201731232","LS_BW_201732633","LS_BW_201734035"]

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



ic = ee.ImageCollection("projects/Sacha/PreprocessedData/L8_Biweekly_V5_Mosaiced_2")

namez = []


myGaps = calculateGaps(ic)
myCollection = ee.ImageCollection(fillCollectionYear(myGaps,ic,1))

myGaps = calculateGaps(myCollection)
myCollection = fillCollectionYear(myGaps,myCollection,2)

myGaps = calculateGaps(myCollection)
myCollection = fillCollectionYear(myGaps,myCollection,3)

myList = myCollection.toList(500)


bandNames = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'swir2','blue_stdDev', 'red_stdDev', 'green_stdDev', 'nir_stdDev', 'swir1_stdDev', 'swir2_stdDev', 'ND_green_swir1_stdDev', 'ND_nir_red_stdDev', 'ND_nir_swir2_stdDev', 'svvi', 'thermal', 'date', 'year', 'cloudMask', 'TDOMMask', 'pixel_qa'])


studyArea =  ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
studyArea = studyArea.geometry().bounds()

geom = studyArea.getInfo()


for i in range(0,100,1):
	im = ee.Image(myList.get(i))
	name = filez[i]
	print name
	im = im.select(bandNames)
	task_ordered= ee.batch.Export.image.toAsset(image=im, 
								  description = name, 
								  assetId= "projects/Sacha/PreprocessedData/L8_Biweekly_V5_Mosaiced_gapfilled/"+ name,
								  region=geom['coordinates'], 
								  maxPixels=1e13,
								  crs="epsg:4326",
								  scale=30)
		
	task_ordered.start()





