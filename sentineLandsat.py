# Sentinel-2 package
from paramsTemplate import *
import ee
from Py6S import *
import math
import datetime
import os, sys
from utils import *
import sun_angles
import view_angles
import time

class env(object):

	def __init__(self):
		"""Initialize the environment."""

		# Initialize the Earth Engine object, using the authentication credentials.
		ee.Initialize()
		
		self.dem = ee.Image("JAXA/ALOS/AW3D30_V1_1").select(["AVE"])
		self.epsg = "EPSG:4326"
				
		##########################################
		# variable for the landsat data request #
		##########################################
		self.metadataCloudCoverMax = 30;

		##########################################
		# Export variables		  		         #
		##########################################		

		#self.assetId ="projects/Sacha/PreprocessedData/L8_Biweekly_V4/"
		self.assetId ="projects/Sacha/PreprocessedData/TOA_composites"
		self.name = "LSS2_ECUADOR_ANNUAL_MEDIAN_2018_000365" 
		self.exportScale = 20		
	
		##########################################
		# variable for the shadowMask  algorithm #
		##########################################

		# zScoreThresh: Threshold for cloud shadow masking- lower number masks out 
		# less. Between -0.8 and -1.2 generally works well
		self.zScoreThresh = -0.9

		# shadowSumThresh: Sum of IR bands to include as shadows within TDOM and the 
		# shadow shift method (lower number masks out less)
		self.shadowSumThresh = 0.30;
		
		# contractPixels: The radius of the number of pixels to contract (negative buffer) clouds and cloud shadows by. Intended to eliminate smaller cloud 
		#    patches that are likely errors (1.5 results in a -1 pixel buffer)(0.5 results in a -0 pixel buffer)
		# (1.5 or 2.5 generally is sufficient)
		self.contractPixels = 1.5; 
		
		# dilatePixels: The radius of the number of pixels to dilate (buffer) clouds 
		# and cloud shadows by. Intended to include edges of clouds/cloud shadows 
		# that are often missed (1.5 results in a 1 pixel buffer)(0.5 results in a 0 pixel buffer)
		# (2.5 or 3.5 generally is sufficient)
		self.dilatePixels = 3.25;	
		
		
		##########################################
		# variable for cloudScore  algorithm     #
		##########################################	
		
		# 9. Cloud and cloud shadow masking parameters.
		# If cloudScoreTDOM is chosen
		# cloudScoreThresh: If using the cloudScoreTDOMShift method-Threshold for cloud 
		#    masking (lower number masks more clouds.  Between 10 and 30 generally works best)
		self.cloudScoreThresh = 30;
		
		# Percentile of cloud score to pull from time series to represent a minimum for 
		# the cloud score over time for a given pixel. Reduces commission errors over 
		# cool bright surfaces. Generally between 5 and 10 works well. 0 generally is a bit noisy	
		self.cloudScorePctl = 8 	
		
		##########################################
		# variable for terrain  algorithm        #
		##########################################		
		
		self.terrainScale = 600
		
		##########################################
		# variable band selection  		         #
		##########################################		
		
		self.divideBandsLandsat = ee.List(['blue','green','red','nir','swir1','swir2'])
		self.bandNamesLandsat = ee.List(['blue','green','red','nir','swir1','thermal','swir2','qa'])
		self.sensorBandDictLandsatSR = ee.Dictionary({'L8' : ee.List([1,2,3,4,5,9,6,10]),\
                                                      'L7' : ee.List([0,1,2,3,4,5,8,9])})
		

		self.s2BandsIn = ee.List(['QA60','B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12'])
		self.s2BandsOut = ee.List(['QA60','cb','blue','green','red','re1','re2','re3','nir','nir2','waterVapor','cirrus','swir1','swir2'])
		self.divideBands = ee.List(['blue','green','red','re1','re2','re3','nir','nir2','cb','cirrus','swir1','swir2','waterVapor'])
		self.medianIncludeBands = ee.List(['blue','green','red','re1','re2','re3','nir','nir2','cb','cirrus','swir1','swir2','waterVapor'])

		
		##########################################
		# enable / disable modules 		         #
		##########################################		  
		
		self.cloudMask = True
		self.QAcloudMask = True
		self.shadowMask = True
		self.brdfCorrect = False
		self.terrainCorrection = False



class Landsat():      
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 
	
	def main(self,studyArea,startDate,endDate,startDay,endDay,week,regionName):
		
		self.env.startDate = startDate
		self.env.endDate = endDate


		self.env.startDoy = startDay
		self.env.endDoy = endDay
		self.env.regionName = regionName
		
		self.studyArea = studyArea
	
		landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterDate(self.env.startDate,self.env.endDate).filterBounds(studyArea)
		landsat8 = landsat8.filterMetadata('CLOUD_COVER','less_than',self.env.metadataCloudCoverMax)
		landsat8 = landsat8.select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat)

		landsat7 =  ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterDate(self.env.startDate,self.env.endDate).filterBounds(studyArea)
		landsat7 = landsat7.filterMetadata('CLOUD_COVER','less_than',self.env.metadataCloudCoverMax)
		landsat7 = landsat7.select(self.env.sensorBandDictLandsatSR.get('L7'),self.env.bandNamesLandsat)
		
		landsat7 = landsat7.map(self.HarmonizeLandsat7)
		landsat8 = landsat8.map(self.HarmonizeLandsat8)
	
		#landsat = landsat7.merge(landsat8)
		landsat = landsat8
		  
		if landsat.size().getInfo() > 0:
			
			# mask clouds using the QA band
			#if self.env.maskSR == True:
				#print "removing clouds" 
				#landsat = landsat.map(self.CloudMaskSRL8)    
			
			# mask clouds using cloud mask function
			if self.env.shadowMask == True:
				#print "shadow masking"
				self.fullCollection = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterBounds(studyArea).filter(ee.Filter.lt("CLOUD_COVER",30))\
																				   .select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat) 
																				   
				l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA').filterBounds(studyArea).filter(ee.Filter.lt("CLOUD_COVER",30))\
																				   .select(self.env.sensorBandDictLandsatSR.get('L7'),self.env.bandNamesLandsat) 																		

				self.fullCollection = self.fullCollection.merge(l7)																			
				landsat = self.maskShadows(landsat)
			
			#landsat = landsat.map(self.scaleLandsat).map(self.addDateYear)
			landsat = landsat.map(self.addDateYear)

			# mask clouds using cloud mask function
			if self.env.cloudMask == True:
				#print "removing some more clouds"
				landsat = landsat.map(self.maskClouds)

			if self.env.brdfCorrect == True:
				landsat = landsat.map(self.brdf)

						
			if self.env.terrainCorrection == True:
				landsat = ee.ImageCollection(landsat.map(self.terrain))

		
			return landsat

	def addDateYear(self,img):
		#add a date and year band
		date = ee.Date(img.get("system:time_start"))
		
		day = date.getRelative('day','year').add(1);
		yr = date.get('year');
		mk = img.mask().reduce(ee.Reducer.min());
	  	
	  	img = img.addBands(ee.Image.constant(day).mask(mk).uint16().rename(['date']));
	  	img = img.addBands(ee.Image.constant(yr).mask(mk).uint16().rename(['year']));
		
		return img;
	
	def CloudMaskSRL8(self,img):
		"""apply cf-mask Landsat""" 
		QA = img.select("qa")
		
		cloud =  QA.bitwiseAnd(4).neq(0);
		
		return img.updateMask(cloud.Not()).copyProperties(img)		
         
	def scaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
		thermal = img.select(ee.List(['thermal'])).multiply(0.1)
		scaled = ee.Image(img).select(self.env.divideBandsLandsat).multiply(ee.Number(0.0001))
		
		return img.select(['TDOMMask']).addBands(scaled).addBands(thermal)
		
	def reScaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
		noScaleBands = ee.List(['date','year','TDOMMask','cloudMask','count'])
		noScale = ee.Image(img).select(noScaleBands)
		thermalBand = ee.List(['thermal'])
		thermal = ee.Image(img).select(thermalBand).multiply(10)
                
		otherBands = ee.Image(img).bandNames().removeAll(thermalBand).removeAll(noScaleBands)
		scaled = ee.Image(img).select(otherBands).divide(0.0001)
        
		image = ee.Image(scaled.addBands([thermal,noScale])).int16()
        
		return image.copyProperties(img)


	def maskClouds(self,img):
		"""
		Computes spectral indices of cloudyness and take the minimum of them.
		
		Each spectral index is fairly lenient because the group minimum 
		is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)
		originally written by Matt Hancher for Landsat imageryadapted to Sentinel by Chris Hewig and Ian Housman
		"""
		
		score = ee.Image(1.0);
		# Clouds are reasonably bright in the blue band.
		blue_rescale = img.select('blue').subtract(ee.Number(0.1)).divide(ee.Number(0.3).subtract(ee.Number(0.1)))
		score = score.min(blue_rescale);

		# Clouds are reasonably bright in all visible bands.
		visible = img.select('red').add(img.select('green')).add(img.select('blue'))
		visible_rescale = visible.subtract(ee.Number(0.2)).divide(ee.Number(0.8).subtract(ee.Number(0.2)))
		score = score.min(visible_rescale);

		# Clouds are reasonably bright in all infrared bands.
		infrared = img.select('nir').add(img.select('swir1')).add(img.select('swir2'))
		infrared_rescale = infrared.subtract(ee.Number(0.3)).divide(ee.Number(0.8).subtract(ee.Number(0.3)))
		score = score.min(infrared_rescale);

		# Clouds are reasonably cool in temperature.
		temp_rescale = img.select('thermal').subtract(ee.Number(300)).divide(ee.Number(290).subtract(ee.Number(300)))
		score = score.min(temp_rescale);

		# However, clouds are not snow.
		ndsi = img.normalizedDifference(['green', 'swir1']);
		ndsi_rescale = ndsi.subtract(ee.Number(0.8)).divide(ee.Number(0.6).subtract(ee.Number(0.8)))
		score =  score.min(ndsi_rescale).multiply(100).byte();
		mask = score.lt(self.env.cloudScoreThresh).rename(['cloudMask']);
		img = img.updateMask(mask).addBands([mask]);
        
		return img;
        
	def maskShadows(self,collection):

		def TDOM(image):
			zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
			irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
			TDOMMask = zScore.lt(self.env.zScoreThresh).reduce(ee.Reducer.sum()).eq(2)\
				.And(irSum.lt(self.env.shadowSumThresh)).Not()
			TDOMMask = TDOMMask.focal_min(self.env.contractPixels).focal_max(self.env.dilatePixels).rename(['TDOMMask'])
			
			image = image.addBands([TDOMMask])

			return image.updateMask(TDOMMask)
			
		shadowSumBands = ['nir','swir1']

		# Get some pixel-wise stats for the time series
		irStdDev = self.fullCollection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
		irMean = self.fullCollection.select(shadowSumBands).reduce(ee.Reducer.mean())

		# Mask out dark dark outliers
		collection_tdom = collection.map(TDOM)

		return collection_tdom

	def terrain(self,img):   
		
		
		degree2radian = 0.01745;
		otherBands = img.select(['thermal','date','year','TDOMMask','cloudMask'])

		def topoCorr_IC(img):
			
			dem = ee.Image("USGS/SRTMGL1_003")
			
			
			# Extract image metadata about solar position
			SZ_rad = ee.Image.constant(ee.Number(img.get('SOLAR_ZENITH_ANGLE'))).multiply(degree2radian).clip(img.geometry().buffer(10000)); 
			SA_rad = ee.Image.constant(ee.Number(img.get('SOLAR_AZIMUTH_ANGLE'))).multiply(degree2radian).clip(img.geometry().buffer(10000)); 
			
				
			# Creat terrain layers
			slp = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000));
			slp_rad = ee.Terrain.slope(dem).multiply(degree2radian).clip(img.geometry().buffer(10000));
			asp_rad = ee.Terrain.aspect(dem).multiply(degree2radian).clip(img.geometry().buffer(10000));
  
  
			
			# Calculate the Illumination Condition (IC)
			# slope part of the illumination condition
			cosZ = SZ_rad.cos();
			cosS = slp_rad.cos();
			slope_illumination = cosS.expression("cosZ * cosS", \
												{'cosZ': cosZ, 'cosS': cosS.select('slope')});
			
			
			# aspect part of the illumination condition
			sinZ = SZ_rad.sin(); 
			sinS = slp_rad.sin();
			cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
			aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", \
                                           {'sinZ': sinZ, \
                                            'sinS': sinS, \
                                            'cosAziDiff': cosAziDiff});
			
			# full illumination condition (IC)
			ic = slope_illumination.add(aspect_illumination);
			
			

			# Add IC to original image
			img_plus_ic = ee.Image(img.addBands(ic.rename(['IC'])).addBands(cosZ.rename(['cosZ'])).addBands(cosS.rename(['cosS'])).addBands(slp.rename(['slope'])));
			
			return ee.Image(img_plus_ic);
 
		def topoCorr_SCSc(img):
			img_plus_ic = img;
			mask1 = img_plus_ic.select('nir').gt(-0.1);
			mask2 = img_plus_ic.select('slope').gte(5) \
                            .And(img_plus_ic.select('IC').gte(0)) \
                            .And(img_plus_ic.select('nir').gt(-0.1));

			img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));

			bandList = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']; # Specify Bands to topographically correct


			def applyBands(image):
				blue = apply_SCSccorr('blue').select(['blue'])
				green = apply_SCSccorr('green').select(['green'])
				red = apply_SCSccorr('red').select(['red'])
				nir = apply_SCSccorr('nir').select(['nir'])
				swir1 = apply_SCSccorr('swir1').select(['swir1'])
				swir2 = apply_SCSccorr('swir2').select(['swir2'])
				return replace_bands(image, [blue, green, red, nir, swir1, swir2])

			def apply_SCSccorr(band):
				method = 'SCSc';
						
				out = ee.Image(1).addBands(img_plus_ic_mask2.select('IC', band)).reduceRegion(reducer= ee.Reducer.linearRegression(2,1), \
  																	   geometry= ee.Geometry(img.geometry().buffer(-5000)), \
																		scale= self.env.terrainScale, \
																		bestEffort =True,
																		maxPixels=1e10)
																		
				
				#out_a = ee.Number(out.get('scale'));
				#out_b = ee.Number(out.get('offset'));
				#out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')));
				
				
				fit = out.combine({"coefficients": ee.Array([[1],[1]])}, False);

				#Get the coefficients as a nested list, 
				#cast it to an array, and get just the selected column
				out_a = (ee.Array(fit.get('coefficients')).get([0,0]));
				out_b = (ee.Array(fit.get('coefficients')).get([1,0]));
				out_c = out_a.divide(out_b)

					
				# apply the SCSc correction
				SCSc_output = img_plus_ic_mask2.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
																'image': img_plus_ic_mask2.select([band]),
																'ic': img_plus_ic_mask2.select('IC'),
																'cosB': img_plus_ic_mask2.select('cosS'),
																'cosZ': img_plus_ic_mask2.select('cosZ'),
																'cvalue': out_c });
		  
				return ee.Image(SCSc_output);
																	  
			#img_SCSccorr = ee.Image([apply_SCSccorr(band) for band in bandList]).addBands(img_plus_ic.select('IC'));
			img_SCSccorr = applyBands(img).select(bandList).addBands(img_plus_ic.select('IC'))
		
			bandList_IC = ee.List([bandList, 'IC']).flatten();
			
			img_SCSccorr = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList);
  			
			return img_SCSccorr.unmask(img_plus_ic.select(bandList)) 
	
		
		
		img = topoCorr_IC(img)
		img = topoCorr_SCSc(img)
		
		return img.addBands(otherBands)

	def defringe(self,img):
		
		# threshold for defringing landsat5 and 7
		fringeCountThreshold = 279

		k = ee.Kernel.fixed(41, 41, 
                                [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]);


		m = ee.Image(img).mask().reduce(ee.Reducer.min())
		sum = m.reduceNeighborhood(ee.Reducer.sum(), k, 'kernel')
		mask = sum.gte(fringeCountThreshold)        
        
		return img.updateMask(mask)

 	def brdf(self,img):   
		
		import sun_angles
		import view_angles

	
		def _apply(image, kvol, kvol0):
			blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
			green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
			red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
			nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
			swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
			swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
			return replace_bands(image, [blue, green, red, nir, swir1, swir2])


		def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
			"""fiso + fvol * kvol + fgeo * kgeo"""
			iso = ee.Image(f_iso)
			geo = ee.Image(f_geo)
			vol = ee.Image(f_vol)
			pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
			pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
			cfac = pred0.divide(pred).rename(['cfac'])
			corr = image.select(band_name).multiply(cfac).rename([band_name])
			return corr


		def _kvol(sunAz, sunZen, viewAz, viewZen):
			"""Calculate kvol kernel.
			From Lucht et al. 2000
			Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""
			
			relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
			pa1 = viewZen.cos() \
				.multiply(sunZen.cos())
			pa2 = viewZen.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle1 = pa1.add(pa2)
			phase_angle = phase_angle1.acos()
			p1 = ee.Image(PI().divide(2)).subtract(phase_angle)
			p2 = p1.multiply(phase_angle1)
			p3 = p2.add(phase_angle.sin())
			p4 = sunZen.cos().add(viewZen.cos())
			p5 = ee.Image(PI().divide(4))

			kvol = p3.divide(p4).subtract(p5).rename(['kvol'])

			viewZen0 = ee.Image(0)
			pa10 = viewZen0.cos() \
				.multiply(sunZen.cos())
			pa20 = viewZen0.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle10 = pa10.add(pa20)
			phase_angle0 = phase_angle10.acos()
			p10 = ee.Image(PI().divide(2)).subtract(phase_angle0)
			p20 = p10.multiply(phase_angle10)
			p30 = p20.add(phase_angle0.sin())
			p40 = sunZen.cos().add(viewZen0.cos())
			p50 = ee.Image(PI().divide(4))

			kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])

			return (kvol, kvol0)
         
 		date = img.date()
		footprint = determine_footprint(img)
		(sunAz, sunZen) = sun_angles.create(date, footprint)
		(viewAz, viewZen) = view_angles.create(footprint)
		(kvol, kvol0) = _kvol(sunAz, sunZen, viewAz, viewZen)
		return _apply(img, kvol.multiply(PI()), kvol0.multiply(PI()))
	
	def medoidMosaic(self,collection):
		""" medoid composite with equal weight among indices """
		nImages = ee.ImageCollection(collection).select([0]).count().rename('count')
		bandNames = ee.Image(collection.first()).bandNames()
		otherBands = bandNames.removeAll(self.env.divideBands)

		others = collection.select(otherBands).reduce(ee.Reducer.mean()).rename(otherBands);
		
		collection = collection.select(self.env.divideBands)

		bandNumbers = ee.List.sequence(1,self.env.divideBands.length());

		median = ee.ImageCollection(collection).median()
        
		def subtractmedian(img):
			diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
			return diff.reduce('sum').addBands(img);
        
		medoid = collection.map(subtractmedian)
  
		medoid = ee.ImageCollection(medoid).reduce(ee.Reducer.min(self.env.divideBands.length().add(1))).select(bandNumbers,self.env.divideBands);
  
		return medoid.addBands(others).addBands(nImages);		

	def medianMosaic(self,collection):
		
		""" median composite """ 
		median = collection.select(medianIncludeBands).median();
		othersBands = bandNames.removeAll(medianIncludeBands);
		others = collection.select(otherBands).mean();
    
		return median.addBands(others)

	def setMetaData(self,img):

			img = ee.Image(img).set({'regionName': str(self.env.regionName),
									 'system:time_start':ee.Date(self.env.startDate).millis(),
									 'startDOY':str(self.env.startDoy),
									 'endDOY':str(self.env.endDoy),
									 'assetId':str(self.env.assetId),
									 'compositingMethod':'medoid',
									 'toaOrSR':'SR',
									 'epsg':str(self.env.epsg),
									 'exportScale':str(self.env.exportScale),
									 'shadowSumThresh':str(self.env.shadowSumThresh),
									 'maskSR':str(self.env.maskSR),
									 'cloudMask':str(self.env.cloudMask),
									 'hazeMask':str(self.env.hazeMask),
									 'shadowMask':str(self.env.shadowMask),
									 'brdfCorrect':str(self.env.brdfCorrect),
									 'terrainCorrection':str(self.env.terrainCorrection),
									 'contractPixels':str(self.env.contractPixels),
									 'dilatePixels':str(self.env.dilatePixels),
									 'zScoreThresh':str(self.env.zScoreThresh),
									 'metadataCloudCoverMax':str(self.env.metadataCloudCoverMax),
									 'cloudScorePctl':str(self.env.cloudScorePctl),
									 'hazeThresh':str(self.env.hazeThresh),
									 'terrainScale':str(self.env.terrainScale)})

			return img



	def HarmonizeLandsat7(self,img):
		
		t = img.get("system:time_start")
		bandNames = ee.List(['blue','green','red','nir','swir1','swir2']);
		slopes = [1.10601,0.99091,1.05681,1.0045,1.03611,1.04011]
		intercepts = [-0.0139,0.00411,-0.0024,-0.0076,0.00411,0.00861]	

		otherBands = ee.Image(img).bandNames().removeAll(bandNames)
		others = img.select(otherBands)

		img = ee.Image(img).select(bandNames).multiply(slopes).add(intercepts).float();
		
		return img.addBands(others).copyProperties(img).set("system:time_start",t)
	
	def HarmonizeLandsat8(self,img):
		
		t = img.get("system:time_start")
		bandNames = ee.List(['blue','green','red','nir','swir1','swir2']);
		slopes = [1.0946,1.0043,1.0524,0.8954,1.0049,1.0002];
		intercepts = [-0.0107,0.0026,-0.0015,0.0033,0.0065,0.0046];		
		
		otherBands = ee.Image(img).bandNames().removeAll(bandNames)
		others = img.select(otherBands)
		
		img = ee.Image(img).select(bandNames).multiply(slopes).add(intercepts).float();
				
		return img.addBands(others).copyProperties(img).set("system:time_start",t)




class sentinel2():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 
	
	
	def main(self,studyArea,startDate,endDate,startDay,endDay,week,regionName):
		
		self.env.regionName = regionName
		self.env.startDate = startDate
		self.env.endDate = endDate
		
		self.env.startDoy = startDay
		self.env.endDoy = endDay
		
		s2 = self.getSentinel2(startDate,endDate,studyArea);
		print(s2.size().getInfo())

		if s2.size().getInfo() > 0:		
			
			
			s2 = s2.map(self.scaleS2)
			
			# masking the shadows
			if self.env.shadowMask == True:
				s2 = self.maskShadows(s2,studyArea)
			
			self.collectionMeta = s2.getInfo()['features']
			
			s2 = s2.select(self.env.s2BandsIn,self.env.s2BandsOut).map(self.addDateYear)
			
			if self.env.QAcloudMask == True:
				s2 = s2.map(self.QAMaskCloud)
	
			if self.env.cloudMask == True:
				s2 = s2.map(self.sentinelCloudScore)
				s2 = self.cloudMasking(s2)

			if self.env.brdfCorrect == True:
				s2 = s2.map(self.brdf)
					
			if self.env.terrainCorrection == True:
				s2 = s2.map(self.getTopo)
				corrected = s2.filter(ee.Filter.gt("slope",20))
				notCorrected = s2.filter(ee.Filter.lt("slope",20))
				s2 = corrected.map(self.terrain).merge(notCorrected)			
					
			return s2


	def getSentinel2(self,start,end,studyArea):
	
		s2s = ee.ImageCollection('COPERNICUS/S2').filterDate(start,end) \
	                                             .filterBounds(studyArea) \
												 .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',self.env.metadataCloudCoverMax)) \
												 .filter(ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT',self.env.metadataCloudCoverMax))\
		
		return s2s
	
	def addDateYear(self,img):
		#add a date and year band
		date = ee.Date(img.get("system:time_start"))
		
		day = date.getRelative('day','year').add(1);
		yr = date.get('year');
		mk = img.mask().reduce(ee.Reducer.min());
	  	
		img = img.addBands(ee.Image.constant(day).mask(mk).uint16().rename('date'));
		img = img.addBands(ee.Image.constant(yr).mask(mk).uint16().rename('year'));
		
		return img;

	def maskShadows(self,collection,studyArea):

		def TDOM(image):
			zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
			irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
			TDOMMask = zScore.lt(self.env.zScoreThresh).reduce(ee.Reducer.sum()).eq(2)\
							 .And(irSum.lt(self.env.shadowSumThresh)).Not()
			TDOMMask = TDOMMask.focal_min(self.env.dilatePixels)
			return image.addBands(TDOMMask.rename(['TDOMMask']))

		def mask(image):
			outimg = image.updateMask(image.select(['TDOMMask']))
			return outimg

		shadowSumBands = ['B8','B11']

		allCollection = ee.ImageCollection('COPERNICUS/S2').filterBounds(studyArea).filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE",30))
	                                           
		# Get some pixel-wise stats for the time series
		irStdDev = allCollection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
		irMean = allCollection.select(shadowSumBands).reduce(ee.Reducer.mean())

		# Mask out dark dark outliers
		collection_tdom = collection.map(TDOM)

		return collection_tdom.map(mask)

	def getTopo(self,img):
		''' funtion to filter for areas with terrain and areas without'''
		dem = self.env.dem.unmask(0)
		geom = ee.Geometry(img.get('system:footprint')).bounds()
		slp_rad = ee.Terrain.slope(dem).clip(geom);
		
		slope = slp_rad.reduceRegion(reducer= ee.Reducer.percentile([80]), \
									 geometry= geom,\
									 scale= 100 ,\
									 maxPixels=10000000)
		return img.set('slope',slope.get('slope'))

	def scaleS2(self,img):
		
		divideBands = ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12']
		bandNames = img.bandNames()
		otherBands = bandNames.removeAll(divideBands)
		
		others = img.select(otherBands)
		out = img.select(divideBands).divide(10000)
		return out.addBands(others).copyProperties(img,['system:time_start','system:footprint','MEAN_SOLAR_ZENITH_ANGLE','MEAN_SOLAR_AZIMUTH_ANGLE']).set("centroid",img.geometry().centroid());

	def reScaleS2(self,img):
		
		bandNames = img.bandNames()
		otherBands = bandNames.removeAll(self.env.divideBands)
		others = img.select(otherBands)
				
		t = img.select(self.env.divideBands);
		t = t.multiply(10000)
						
		out = ee.Image(t.copyProperties(img).copyProperties(img,['system:time_start'])).addBands(others).int16()

		return out;

	def pixelArea(self,img):
		geom = ee.Geometry(img.get('system:footprint')).bounds()

		area = img.select(['red']).gt(0).reduceRegion(reducer= ee.Reducer.sum(),\
							     geometry= geom,\
							     scale= 100,\
							     maxPixels=10000000)
		
		return img.set("pixelArea",area.get("red"))


	# Function to mask clouds using the Sentinel-2 QA band.
	def QAMaskCloud(self,img):
		
		bandNames = img.bandNames()
		otherBands = bandNames.removeAll(self.env.divideBands)
		others = img.select(otherBands)
		
		
		qa = img.select('QA60').int16();		
		
		img = img.select(self.env.divideBands)
		
				
		# Bits 10 and 11 are clouds and cirrus, respectively.
		cloudBitMask = int(math.pow(2, 10));
		cirrusBitMask = int(math.pow(2, 11));
		
		# Both flags should be set to zero, indicating clear conditions.
		mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0));
		
		img = img.updateMask(mask).addBands(others)
		
		# Return the masked and scaled data.
		return img

	def sentinelCloudScore(self,img):
		"""
		Computes spectral indices of cloudyness and take the minimum of them.
    
		Each spectral index is fairly lenient because the group minimum 
		is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)
		
		originally written by Matt Hancher for Landsat imagery
		adapted to Sentinel by Chris Hewig and Ian Housman
		"""

        
		def rescale(img, thresholds):
			"""
			Linear stretch of image between two threshold values.
			"""
			return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])
        
		# cloud until proven otherwise
		score = ee.Image(1)
		blueCirrusScore = ee.Image(0)
		
		# clouds are reasonably bright
		blueCirrusScore = blueCirrusScore.max(rescale(img.select(['blue']), [0.1, 0.5]))
		blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cb']), [0.1, 0.5]))
		blueCirrusScore = blueCirrusScore.max(rescale(img.select(['cirrus']), [0.1, 0.3]))
		score = score.min(blueCirrusScore)
		
		score = score.min(rescale(img.select(['red']).add(img.select(['green'])).add(img.select('blue')), [0.2, 0.8]))
		score = score.min(rescale(img.select(['nir']).add(img.select(['swir1'])).add(img.select('swir2')), [0.3, 0.8]))

		# clouds are moist
		ndsi = img.normalizedDifference(['green','swir1'])
		score=score.min(rescale(ndsi, [0.8, 0.6]))
		score = score.multiply(100).byte();
		score = score.clamp(0,100);
  
		return img.addBands(score.rename(['cloudScore']))
				
	def cloudMasking(self,collection):

		def maskClouds(img):
	
			cloudMask = img.select(['cloudScore']).lt(self.env.cloudScoreThresh)\
												  .focal_min(self.env.dilatePixels) \
												  .focal_max(self.env.contractPixels) \
												  .rename(['cloudMask'])    
            
			bandNames = img.bandNames()
			otherBands = bandNames.removeAll(self.env.divideBands)
			others = img.select(otherBands)
			
			img = img.select(self.env.divideBands).updateMask(cloudMask)
            
			return img.addBands(cloudMask).addBands(others);
		
		
		# Find low cloud score pctl for each pixel to avoid comission errors
		minCloudScore = collection.select(['cloudScore']).reduce(ee.Reducer.percentile([self.env.cloudScorePctl]));
		
		collection = collection.map(maskClouds)
		
		return collection

	def brdf(self,img):   
		

	
		def _apply(image, kvol, kvol0):
			blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
			green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
			red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
			re1 = _correct_band(image, 're1', kvol, kvol0, f_iso=0.2085, f_geo=0.0256, f_vol=0.0845)
			re2 = _correct_band(image, 're2', kvol, kvol0, f_iso=0.2316, f_geo=0.0273, f_vol=0.1003)
			re3 = _correct_band(image, 're3', kvol, kvol0, f_iso=0.2599, f_geo=0.0294, f_vol=0.1197)
			nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
			re4 = _correct_band(image, 're4', kvol, kvol0, f_iso=0.2907, f_geo=0.0410, f_vol=0.1611)
			swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
			swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
			return replace_bands(image, [blue, green, red,re1,re2,re3, nir,re4, swir1, swir2])


		def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
			"""fiso + fvol * kvol + fgeo * kgeo"""
			iso = ee.Image(f_iso)
			geo = ee.Image(f_geo)
			vol = ee.Image(f_vol)
			pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
			pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
			cfac = pred0.divide(pred).rename(['cfac'])
			corr = image.select(band_name).multiply(cfac).rename([band_name])
			return corr


		def _kvol(sunAz, sunZen, viewAz, viewZen):
			"""Calculate kvol kernel.
			From Lucht et al. 2000
			Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""
			
			relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
			pa1 = viewZen.cos() \
				.multiply(sunZen.cos())
			pa2 = viewZen.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle1 = pa1.add(pa2)
			phase_angle = phase_angle1.acos()
			p1 = ee.Image(PI().divide(2)).subtract(phase_angle)
			p2 = p1.multiply(phase_angle1)
			p3 = p2.add(phase_angle.sin())
			p4 = sunZen.cos().add(viewZen.cos())
			p5 = ee.Image(PI().divide(4))

			kvol = p3.divide(p4).subtract(p5).rename(['kvol'])

			viewZen0 = ee.Image(0)
			pa10 = viewZen0.cos() \
				.multiply(sunZen.cos())
			pa20 = viewZen0.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle10 = pa10.add(pa20)
			phase_angle0 = phase_angle10.acos()
			p10 = ee.Image(PI().divide(2)).subtract(phase_angle0)
			p20 = p10.multiply(phase_angle10)
			p30 = p20.add(phase_angle0.sin())
			p40 = sunZen.cos().add(viewZen0.cos())
			p50 = ee.Image(PI().divide(4))

			kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])

			return (kvol, kvol0)
         
		date = img.date()
		footprint =  ee.List(img.geometry().bounds().bounds().coordinates().get(0));
		(sunAz, sunZen) = sun_angles.create(date, footprint)
		(viewAz, viewZen) = view_angles.create(footprint)
		(kvol, kvol0) = _kvol(sunAz, sunZen, viewAz, viewZen)
		
		bandNames = img.bandNames()
		otherBands = bandNames.removeAll(self.env.divideBands)
		others = img.select(otherBands)
						
		img = ee.Image(_apply(img, kvol.multiply(PI()), kvol0.multiply(PI())))
				
		return img

	def terrain(self,img):
		degree2radian = 0.01745;

		geom = ee.Geometry(img.get('system:footprint')).bounds().buffer(10000) 
		dem = self.env.dem
	
		bandNames = img.bandNames()
		otherBands = bandNames.removeAll(self.env.divideBands)
		others = img.select(otherBands)

		bandList = ['blue','green','red','re1','re2','re3','nir','re4','cb','cirrus','swir1','swir2','waterVapor']
		
		def topoCorr_IC(img):
		
			# Extract image metadata about solar position
			SZ_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_ZENITH_ANGLE'))).multiply(degree2radian).clip(geom); 
			SA_rad = ee.Image.constant(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))).multiply(degree2radian).clip(geom); 

			
			# Creat terrain layers
			slp = ee.Terrain.slope(dem).clip(geom);
			slp_rad = ee.Terrain.slope(dem).multiply(degree2radian).clip(geom);
			asp_rad = ee.Terrain.aspect(dem).multiply(degree2radian).clip(geom);
  			

			# Calculate the Illumination Condition (IC)
			# slope part of the illumination condition
			cosZ = SZ_rad.cos();
			cosS = slp_rad.cos();
			slope_illumination = cosS.expression("cosZ * cosS", \
												{'cosZ': cosZ, 'cosS': cosS.select('slope')});
	
			# aspect part of the illumination condition
			sinZ = SZ_rad.sin(); 
			sinS = slp_rad.sin();
			cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
			aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", \
											 {'sinZ': sinZ, \
                                              'sinS': sinS, \
                                              'cosAziDiff': cosAziDiff});
			
			# full illumination condition (IC)
			ic = slope_illumination.add(aspect_illumination);
			
			# Add IC to original image
			img_plus_ic = ee.Image(img.addBands(ic.rename(['IC'])).addBands(cosZ.rename(['cosZ'])).addBands(cosS.rename(['cosS'])).addBands(slp.rename(['slope'])));
		
			return ee.Image(img_plus_ic);
 			
			
		def topoCorr_SCSc(img):
			img_plus_ic = img;
			mask1 = img_plus_ic.select('nir').gt(-0.1);
			mask2 = img_plus_ic.select('slope').gte(5) \
							   .And(img_plus_ic.select('IC').gte(0)) \
							   .And(img_plus_ic.select('nir').gt(-0.1));

			img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));

			
		
			def apply_SCSccorr(band):
				method = 'SCSc';
						
				out = ee.Image(1).addBands(img_plus_ic_mask2.select('IC', band)).reduceRegion(reducer= ee.Reducer.linearRegression(2,1), \
  																	   geometry= ee.Geometry(img.geometry().buffer(-5000)), \
																		scale= 300, \
																		bestEffort =True,
																		maxPixels=1e10)
																					
				
				fit = out.combine({"coefficients": ee.Array([[1],[1]])}, False);

				#Get the coefficients as a nested list, 
				#cast it to an array, and get just the selected column
				out_a = (ee.Array(fit.get('coefficients')).get([0,0]));
				out_b = (ee.Array(fit.get('coefficients')).get([1,0]));
				out_c = out_a.divide(out_b)

					
				# apply the SCSc correction
				SCSc_output = img_plus_ic_mask2.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
																'image': img_plus_ic_mask2.select([band]),
																'ic': img_plus_ic_mask2.select('IC'),
																'cosB': img_plus_ic_mask2.select('cosS'),
																'cosZ': img_plus_ic_mask2.select('cosZ'),
																'cvalue': out_c });
		  
				return ee.Image(SCSc_output);
			
			img_SCSccorr = ee.Image([apply_SCSccorr(band) for band in bandList]).addBands(img_plus_ic.select('IC'));
		
			bandList_IC = ee.List([bandList, 'IC']).flatten();
			
			img_SCSccorr = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList);
  			
			return img_SCSccorr.unmask(img_plus_ic.select(bandList)) 			
			
		img = topoCorr_IC(img)
		img = topoCorr_SCSc(img).addBands(others )

		return img




def exportMap(self,img,studyArea,week):

	geom  = studyArea.getInfo();
	sd = str(self.env.startDate.getRelative('day','year').getInfo()).zfill(3);
	ed = str(self.env.endDate.getRelative('day','year').getInfo()).zfill(3);
	year = str(self.env.startDate.get('year').getInfo());
	regionName = self.env.regionName.replace(" ",'_') + "_"

	task_ordered= ee.batch.Export.image.toAsset(image=img, 
							  description = self.env.name + regionName + str(week).zfill(3) +'_'+ year + sd + ed, 
							  assetId= self.env.assetId + self.env.name + regionName + str(week).zfill(3)+'_'+ year + sd + ed,
							  region=geom['coordinates'], 
							  maxPixels=1e13,
							  crs=self.env.epsg,
							  scale=self.env.exportScale)

	task_ordered.start()
	
class Harmonize():
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 
	
	def harmonizeData(self,l8,s2):
		
		bands = ee.List(['blue','green','red','nir','swir1','swir2'])
		s2 = s2.select(bands)
		l8 = l8.select(bands)
		return ee.ImageCollection(l8.merge(s2))
		
	
	def medoidMosaic(self,collection):
		""" medoid composite with equal weight among indices """
	
		bands = ee.List(['blue','green','red','nir','swir1','swir2'])
		collection = collection.select(bands)
		median = ee.ImageCollection(collection).median()
        
		bandNumbers = ee.List.sequence(1,bands.length());
        
		def subtractmedian(img):
			diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
			return diff.reduce('sum').addBands(img);
        
		medoid = collection.map(subtractmedian)
  
		medoid = ee.ImageCollection(medoid).reduce(ee.Reducer.min(bands.length().add(1))).select(bandNumbers,bands)
  
		return medoid;

	def exportMap(self,img,studyArea,week):

		geom  = studyArea.bounds().getInfo();
		
		task_ordered= ee.batch.Export.image.toAsset(image=img.clip(studyArea.buffer(10000)), 
						  description = self.env.name + str(week), 
						  assetId= self.env.assetId + self.env.name ,
						  region=geom['coordinates'], 
						  maxPixels=1e13,
						  crs=self.env.epsg,
						  scale=self.env.exportScale)
	
		task_ordered.start() 
	

if __name__ == "__main__":        

	ee.Initialize()
		
	
	start = 0

	for i in range(0,1,1):
		#2018 starts at week 104
		startWeek = start+ i
		print "week", startWeek
	
		year = ee.Date("2018-01-01")
		startDay = 0 #(startWeek -1) *14
		endDay = 365 #(startWeek) *31 -1
		print startDay, endDay

		startDate = year.advance(startDay,'day') 
		endDate = year.advance(endDay,'day')

		regionName = 'SIERRA'
		studyArea =  ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
		studyArea = studyArea.filterMetadata('PROVINCIA','equals',regionName).geometry().bounds()
		
		studyArea =  ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
		studyArea = studyArea.geometry().bounds()
	
		landsat = Landsat().main(studyArea,startDate,endDate,startDay,endDay,startWeek,regionName)
		s2 = sentinel2().main(studyArea,startDate,endDate,startDay,endDay,startWeek,regionName)
	
		collection = Harmonize().harmonizeData(landsat,s2)
		#print(collection.size().getInfo())
		
		#medoid = Harmonize().medoidMosaic(collection)
		#medoid = medoid.set("system:time_start",startDate.millis())
		
		medoid = collection.median().set("system:time_start",startDate.millis())
		percentiles = collection.reduce(ee.Reducer.percentile([25,75]));
		inbands =  ['blue_p25', 'blue_p75', 'green_p25', 'green_p75', 'red_p25', 'red_p75', 'nir_p25', 'nir_p75', 'swir1_p25', 'swir1_p75', 'swir2_p25', 'swir2_p75']
		outbands = ['p25_blue', 'p75_blue', 'p25_green', 'p75_green', 'p25_red', 'p75_red', 'p25_nir', 'p75_nir', 'p25_swir1', 'p75_swir1', 'p25_swir2', 'p75_swir2']
		
		percentiles = percentiles.select(inbands,outbands)
		
		#print percentiles.bandNames().getInfo()
		nImages = ee.ImageCollection(collection).select([0]).count().rename('count')
		
		medoid = medoid.addBands(percentiles)
		medoid = medoid.multiply(10000).int16()
		medoid = medoid.addBands(nImages)
		Harmonize().exportMap(medoid,studyArea,startWeek)
		
		#print(medoid.bandNames().getInfo())


