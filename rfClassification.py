import ee
from addCovariates import *

class randomForest():

	def __init__(self):
		ee.Initialize()

		self.exportPath = 'users/TEST/'
		self.epsg = "EPSG:32717"
		# set number of trees for the random forest classifier
		self.trees = 40
		
		# set the version number of the product
		self.versionNumber = 9;

		# specify the training data set
		self.TrainingID = 'CompiledData02102019'

		# // specify the Landsat and Sentinel 1 data set
		self.LandsatID = 'TOA_composites_V2'
		self.DualS1_ID = 'S1Dual_Annual_V2'
		#  specify the sample dataset
		self.trainingData = ee.FeatureCollection("projects/Sacha/AncillaryData/RefData/CompiledData02102019")
		# // specify the Sentinel 1 and Landsat S2 data set
		self.Dual_Annual = ee.ImageCollection("projects/Sacha/PreprocessedData/S1Dual_Annual_V2")
		self.S2LSTOA = ee.ImageCollection("projects/Sacha/PreprocessedData/TOA_composites_V2")
		# Specify the regions feature collection
		self.ecoregions = ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
		
		# todo need to convert covariates into py
		self.covariates = addCovariates()
		self.bandNames = ee.List(["red", "nir", "swir1", "swir2", 

							"ND_nir_red", "ND_nir_swir1", "ND_nir_swir2", 
							"ND_red_swir1", "ND_red_swir2", 
							"ND_swir1_swir2",    

							"p25_nir", "p25_red", "p25_swir1", "p25_swir2", 
							"p25_ND_nir_red", "p25_ND_nir_swir1", "p25_ND_nir_swir2",  
							"p25_ND_red_swir1", "p25_ND_red_swir2", 
							"p25_ND_swir1_swir2", 

							"p75_nir", "p75_red", "p75_swir1", "p75_swir2",
							"p75_ND_nir_red", "p75_ND_nir_swir1", "p75_ND_nir_swir2", 
							"p75_ND_red_swir1", "p75_ND_red_swir2", 
							"p75_ND_swir1_swir2", 

							"VH", "VH_amplitude_1", "VH_phase_1", "VV", "VV_amplitude_1","VV_phase_1","ratio",

							"aspect", "eastness", "elevation", "northness", "slope", 

							"change_abs", "change_norm", "max_extent","occurrence", "transition", "seasonality",
							'svvi','blue_stdDev','red_stdDev','green_stdDev','nir_stdDev','swir1_stdDev','swir2_stdDev',
							'ND_green_swir1_stdDev','ND_nir_red_stdDev','ND_nir_swir2_stdDev','clusters']);


	def main(self,year,ProcessingRegion,**kwargs):

		studyArea = self.ecoregions.filter(ee.Filter.eq("PROVINCIA", ProcessingRegion)).geometry().buffer(1000);
		sampleData = self.trainingData.filterBounds(studyArea)

		imageLandsatS2 = image if kwargs.get('image') else self.S2LSTOA.filterDate(str(year)).first() 
		trianingImg  = self.S2LSTOA.filterDate('2016').first()

		analysisYr = self.addBandsToComposite(imageLandsatS2)
		trainingYr = self.addBandsToComposite(trianingImg)


		data = trainingYr.sampleRegions(sampleData,["niv2_int"],20);
		classifier = ee.Classifier.randomForest(self.trees,0).train(ee.FeatureCollection(data),"niv2_int",self.bandNames);
		classification = analysisYr.classify(classifier,'Mode');

		# todo think if we want to include this if its not automated
		classification = classification.set({
		    'system:time_start': ee.Date.fromYMD(year,6,1).millis(),
		    'version': self.versionNumber,
		    'TrainingData': self.TrainingID,
		    'LandsatID': self.LandsatID,
		    'S1_Dual_ID': self.DualS1_ID,
		    'regionName': ProcessingRegion,
		});

		classificationMMU = self.sieve(classification,11)
		
		self.mmuEx =  kwargs.get('mmu')
		self.exportMap(classification,studyArea,"") if kwargs.get('exportImg') and self.mmuEx == False or self.mmuEx == 3 else None
		self.exportMap(classificationMMU,studyArea,"_MMU") if kwargs.get('exportImg') and self.mmuEx == True or self.mmuEx == 3 else None

		print(classification.bandNames().getInfo())
		return classification, classificationMMU

	def exportMap(self,img,studyArea, suffix):
	    img = img
	    ed = str(year) + suffix
	    regionName = ProcessingRegion.replace(" ",'_')

	    task_ordered= ee.batch.Export.image.toAsset(image=img.clip(studyArea), 
	                              description = regionName + '_LandCover_' + ed, 
	                              assetId = self.exportPath + regionName + '_LandCover_' + ed,
	                              region = studyArea.bounds().getInfo()['coordinates'], 
	                              maxPixels = 1e13,
	                              crs = self.epsg,
	                              scale = 30)

	    task_ordered.start()
	    print('Export Started: ',self.exportPath + regionName + '_LandCover_' + ed)

	# in functions
	def addBandsToComposite(self,composite):
		year = ee.Date(composite.get('system:time_start')).get('year')
		
		stdBands = composite.select(['svvi','blue_stdDev','red_stdDev','green_stdDev',
	                                 'nir_stdDev','swir1_stdDev','swir2_stdDev',
	                                 'ND_green_swir1_stdDev','ND_nir_red_stdDev','ND_nir_swir2_stdDev'])
		img25 = composite.select(['p25_blue','p25_green','p25_red','p25_nir','p25_swir1','p25_swir2'])
		img75 = composite.select(['p75_blue','p75_green','p75_red','p75_nir','p75_swir1','p75_swir2'])

		composite = composite.select(['blue','green','red','nir','swir1','swir2'])
		# adding abnds to the composite

		addJRCAndTopo = self.covariates.addJRCAndTopo(composite)
		# gets nd for selected bands
		img = addJRCAndTopo.addBands(self.covariates.addCovariates("",composite));
		
		# // gets nd of 25 and 75percentile bands
		img = img.addBands(self.covariates.addCovariates("p75_",img75)).addBands(self.covariates.addCovariates("p25_",img25));

		bands = img.bandNames().removeAll(['blue_1','green_1','red_1','nir_1','swir1_1','swir2_1']);

		img = img.select(bands).addBands(stdBands)

		Seg = ee.Algorithms.Image.Segmentation.SNIC(image = img.select(['ND_nir_red']),
													size = 20,
													compactness = 0,
													connectivity = 8);

		img = img.addBands(Seg.select(['clusters']))
		# the dual pol SAR data is added here
		DualSar = ee.Image(self.Dual_Annual.filter(ee.Filter.calendarRange(year,year,'year')).first()).divide(10000);
		img = img.addBands(DualSar)

		return img

	def mosaicUniqueDates(self,collection):
		allDates = ee.List(collection.aggregate_array('system:time_start')).distinct()
		def getUDate(i):
			uD = collection.filterDate(i).mosaic()
			return uD
		out = allDates.map(getUDate)

		return ee.ImageCollection(out) 

	def sieve(self,image,mmu):
		"""function to filter isolated pixels be designated MMU.
			MMU value is the pixel count. """
		connected = image.connectedPixelCount(mmu+20);
		elim = connected.gt(mmu);
		mode = image.focal_mode(mmu/2,'circle');
		mode = mode.mask(image.mask());
		filled = image.where(elim.Not(),mode);
		return filled.rename('Mode_MMU')

if __name__ == "__main__":
	
	# // EcoRegions include: AMAZONIA NOROCCIDENTAL, ANDES DEL NORTE, GALAPAGOS, SIERRA,
	# // CHOCO, PACIFICO ECUATORIAL
	ProcessingRegion = 'AMAZONIA NOROCCIDENTAL'
	year = 2018;

	''' main() takes 2 paramters -year,ProcessingRegion- and returns a tulep of classified images (0: no mmu, 1: mmu)
	optionally you can choose to export by setting exportImg to true and selecting which product you want with 
	mmu. By default main will use the Sentinel2/Landsat Collection, but you can pass in another compostie by setting 
	image = ee.Image('path/to/asset')
	mmu can be 0: exports classification without mmu
				 1: exports classification with mmu
				 3: exports both
	 example: using default image collection and not exporting
	 randomForest().main(2016,'AMAZONIA NOROCCIDENTAL')
	
	randomForest().main(2016,'AMAZONIA NOROCCIDENTAL',exportImg=True,mmu=1)
	exporting mmu
	
	providing your own image and exporting
	myImg = ee.Image('projects/Sacha/.../LSS2_ECUADOR_ANNUAL_MEDIAN_CLDMX80_2016_000365_ST')
	randomForest().main(2016,'AMAZONIA NOROCCIDENTAL',exportImg=False,mmu=True, image=myImg) '''
	randomForest().main(year,ProcessingRegion,exportImg=True,mmu=3)