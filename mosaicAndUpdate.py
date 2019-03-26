from rfClassification import *
import ee


def mosaicAndUpdate(collection,year,exportPath,reduceClass):
	def clipRegion(img):
	  aoi = img.get('regionName');
	  studyRegion = ecoregions.filterMetadata('PROVINCIA','equals',aoi);
	  return img.clip(studyRegion);
	
	rf = randomForest()
	
	ecoregions = ee.FeatureCollection('projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Buffered')

	lcToUpdate = collection.filter(ee.Filter.calendarRange(year,year,'year'))

	mmu = lcToUpdate.filterMetadata('MMU','equals','1 ha')
	print(mmu.size().getInfo())
	if mmu.size().getInfo() < 6:
		no_mmu = lcToUpdate.filterMetadata('MMU','equals','None')
		def sieveImg(img):
			img = rf.sieve(img,11)
			return img.rename('Mode')
		no_mmu = no_mmu.map(sieveImg)


		lcToUpdate = ee.ImageCollection(no_mmu).merge(ee.ImageCollection(mmu)).map(clipRegion).mosaic()
	else:
		lcToUpdate = mmu.map(clipRegion).mosaic()

	print(lcToUpdate.bandNames().getInfo())
	# Hansen should be updated to relevent year
	hansen = ee.Image("UMD/hansen/global_forest_change_2018_v1_6");
	# // MAE layers for updating LC
	rice = ee.FeatureCollection("projects/Sacha/RICE2013")
	shrimp = ee.FeatureCollection("projects/Sacha/shrimpFarm2013")

	# // MAE layers for updating LC
	aw1 = ee.FeatureCollection("projects/Sacha/oriente_a_inudaciont")
	aw2 = ee.FeatureCollection("projects/Sacha/sierra_a_inudacion")
	aw3 = ee.FeatureCollection("projects/Sacha/costa_a_inudacion")
	inf = ee.FeatureCollection("projects/Sacha/Infrastructure")

	# v3 land cover of sierra -used for paramo updating
	v3Sierra = ee.Image('projects/Sacha/Production/LandCover/LC_Nivel2_2016v3/SIERRA_LC_mmu1ha2016')


	# // paint the FCs
	emptImg = ee.Image().byte();
	artificalW = aw1.merge(aw2).merge(aw3);
	riceAg = artificalW.filterBounds(rice)
	shrimpAw = inf.filterBounds(shrimp)

	riceImg = emptImg.paint(riceAg);
	shrimpImg = emptImg.paint(shrimpAw)
	artificalWImg = emptImg.paint(artificalW);
	infra = emptImg.paint(inf);

	forestloss = hansen.select('loss');
	# update using artifical water, forestloss, and infrastructure 
	lcToUpdate = lcToUpdate.where(lcToUpdate.eq(44),5).where(forestloss.eq(1).And(lcToUpdate.eq(40)),43).where(artificalWImg.eq(0),30).where(infra.eq(0),10);
	# //update paramo from v3, shrimp, and rice
	lcToUpdate = lcToUpdate.where(v3Sierra.eq(60),60).where(riceImg.eq(0),5).where(shrimpImg.eq(0),30)
	lcToUpdate = lcToUpdate.set('system:time_start',ee.Date(year))
	print(lcToUpdate.get('system:time_start').getInfo())
	task_ordered= ee.batch.Export.image.toAsset(image=lcToUpdate,
    	description = 'LandCover_' + str(year),
    	assetId = exportPath + 'LandCover_' + str(year),
    	region = ecoregions.geometry().bounds().getInfo()['coordinates'], 
    	maxPixels = 1e13,
    	crs = "EPSG:32717",
    	scale = 30)

	task_ordered.start()

	if reduceClass:
		def reduceClasses(img):
  			img = img.where(img.eq(10).Or(img.eq(1)), 1).where(img.eq(20).Or(img.eq(21)), 2)\
  			.where(img.eq(30).Or(img.eq(3)), 3).where(img.eq(40).Or(img.eq(41)).Or(img.eq(42))\
  			.Or(img.eq(43)), 4).where(img.eq(50).Or(img.eq(51)).Or(img.eq(52)).Or(img.eq(53)), 5)\
  			.where(img.eq(60).Or(img.eq(61)).Or(img.eq(62)), 6);
  			
  			mask = img.gt(6);
			
			return img.mask(mask.Not());
		
		lcToUpdate = reduceClasses(lcToUpdate)
		
		task_ordered= ee.batch.Export.image.toAsset(image=lcToUpdate,
	    	description = 'LandCover_ReduceClasses_' + str(year),
	    	assetId = exportPath + 'LandCover_ReduceClasses_' + str(year),
	    	region = ecoregions.geometry().bounds().getInfo()['coordinates'], 
	    	maxPixels = 1e13,
	    	crs = "EPSG:32717",
	    	scale = 30)

		# task_ordered.start()





if __name__ == '__main__':
	ee.Initialize()
	year = 2018
	collection = ee.ImageCollection('projects/Sacha/Production/LandCover/LC_Nivel2_V8')
	exportPath = 'users/TEST/'
	reduceClass = True

	mosaicAndUpdate(collection, year, exportPath, reduceClass)