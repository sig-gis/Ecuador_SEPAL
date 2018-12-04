# Sentinel-2 package

import ee
from Py6S import *
import math
import datetime
import os, sys
from utils import *
sys.path.append("/gee-atmcorr-S2/bin/")
from atmospheric import Atmospheric
import sun_angles
import view_angles
import time

class env(object):

	def __init__(self):
		"""Initialize the environment."""

		# Initialize the Earth Engine object, using the authentication credentials.
		ee.Initialize()
		
		self.dem =  ee.Image("JAXA/ALOS/AW3D30_V1_1").select(["AVE"])
		self.epsg = "EPSG:32717"	
		self.feature = 0	
		##########################################
		# variable for the getSentinel algorithm #
		##########################################

		self.metadataCloudCoverMax = 80;
		
		
		##########################################
		# variable for the shadowMask  algorithm #
		##########################################
	
		# zScoreThresh: Threshold for cloud shadow masking- lower number masks out 
		# less. Between -0.8 and -1.2 generally works well
		self.zScoreThresh = -1

		# shadowSumThresh: Sum of IR bands to include as shadows within TDOM and the 
		# shadow shift method (lower number masks out less)
		self.shadowSumThresh = 0.500;
		
		# contractPixels: The radius of the number of pixels to contract (negative buffer) clouds and cloud shadows by. Intended to eliminate smaller cloud 
		#    patches that are likely errors (1.5 results in a -1 pixel buffer)(0.5 results in a -0 pixel buffer)
		# (1.5 or 2.5 generally is sufficient)
		self.contractPixels = 1.5; 
		
		# dilatePixels: The radius of the number of pixels to dilate (buffer) clouds 
		# and cloud shadows by. Intended to include edges of clouds/cloud shadows 
		# that are often missed (1.5 results in a 1 pixel buffer)(0.5 results in a 0 pixel buffer)
		# (2.5 or 3.5 generally is sufficient)
		self.dilatePixels = 3.5;	
		
		
		##########################################
		# variable for cloudScore  algorithm     #
		##########################################	
		
		# 9. Cloud and cloud shadow masking parameters.
		# If cloudScoreTDOM is chosen
		# cloudScoreThresh: If using the cloudScoreTDOMShift method-Threshold for cloud 
		#    masking (lower number masks more clouds.  Between 10 and 30 generally works best)
		self.cloudScoreThresh = 20;
		
		# Percentile of cloud score to pull from time series to represent a minimum for 
		# the cloud score over time for a given pixel. Reduces commission errors over 
		# cool bright surfaces. Generally between 5 and 10 works well. 0 generally is a bit noisy	
		self.cloudScorePctl = 5; 	
		
		##########################################
		# variable for terrain  algorithm        #
		##########################################		
		
		self.terrainScale = 1000

		##########################################
		# Export variables		  		         #
		##########################################		

		self.assetId ="projects/Sacha/S2/S2_Biweekly/"
		self.name = "Sentinel2_SR_Biweek_V2_" 
		self.exportScale = 20
		
		##########################################
		# variable band selection  		         #
		##########################################		
		
		self.s2BandsIn = ee.List(['QA60','B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12','TDOMMask'])
		self.s2BandsOut = ee.List(['QA60','cb','blue','green','red','re1','re2','re3','nir','re4','waterVapor','cirrus','swir1','swir2','TDOMMask'])
		self.divideBands = ee.List(['blue','green','red','re1','re2','re3','nir','re4','cb','cirrus','swir1','swir2','waterVapor'])
		self.medianIncludeBands = ee.List(['blue','green','red','re1','re2','re3','nir','re4','cb','cirrus','swir1','swir2','waterVapor'])
		
		##########################################
		# enable / disable modules 		         #
		##########################################		
		self.calcSR = True     
		self.brdf = True
		self.QAcloudMask = True
		self.cloudMask = True
		self.shadowMask = True
		self.terrainCorrection = True


class functions():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
 
	    # get the environment
		self.env = env() 
	
	
	def main(self,studyArea,startDate,endDate,startDay,endDay,week):
		
		self.env.startDate = startDate
		self.env.endDate = endDate
		
		self.env.startDoy = startDay
		self.env.endDoy = endDay
		
		s2 = self.getSentinel2(startDate,endDate,studyArea);
		
		#print(s2.size().getInfo())
		if s2.size().getInfo() > 0:		
			
			print(self.env.startDate.getInfo())
			print(self.env.endDate.getInfo())

			print(ee.Image(s2.first()).get('system:time_start').getInfo())			

			s2 = s2.map(self.scaleS2)
			
			# masking the shadows
			print("Masking shadows..") 
			if self.env.shadowMask == True:
				s2 = self.maskShadows(s2,studyArea)
				
			self.collectionMeta = s2.getInfo()['features']
			#print(ee.Image(s2.first()).get('system:time_start').getInfo())

			print("scaling bands..")
			
			#print(ee.Image(s2.first()).get('system:time_start').getInfo())
			
			
			if self.env.QAcloudMask == True:
				print("use QA band for cloud Masking")
				s2 = s2.map(self.QAMaskCloud)
			#print(ee.Image(s2.first()).get('system:time_start').getInfo())
			

			if self.env.cloudMask == True:
				print("sentinel cloud score...")
				s2 = s2.map(self.sentinelCloudScore)
				s2 = self.cloudMasking(s2)
			#print(ee.Image(s2.first()).get('system:time_start').getInfo())
			
			# applying the atmospheric correction
			if self.env.calcSR  == True:
				s2 = s2.select(self.env.s2BandsOut,self.env.s2BandsIn)
				print("applying atmospheric correction")
				s2 = s2.map(self.TOAtoSR).select(self.env.s2BandsIn,self.env.s2BandsOut)
			#print(ee.Image(s2.first()).get('system:time_start').getInfo())

				
			if self.env.brdf == True:
				print("apply brdf correction..")
				s2 = s2.map(self.brdf)
					
			if self.env.terrainCorrection == True:
				print("apply terrain correction..")
				s2 = s2.map(self.getTopo)
				corrected = s2.filter(ee.Filter.gt("slope",20))
				notCorrected = s2.filter(ee.Filter.lt("slope",20))
				s2 = corrected.map(self.terrain).merge(notCorrected)			
			
			
			
			print("calculating medoid")
			img = self.medoidMosaic(s2)
						
			print("rescale")
			img = self.reScaleS2(img)
						
			print("set MetaData")
			img = self.setMetaData(img)
			
			print("exporting composite")
			self.exportMap(img,studyArea,week)


	def getSentinel2(self,start,end,studyArea):
	
		s2s = ee.ImageCollection('COPERNICUS/S2').filterDate(start,end) \
	                                             .filterBounds(studyArea) \
												 .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',self.env.metadataCloudCoverMax)) \
												 .filter(ee.Filter.lt('CLOUD_COVERAGE_ASSESSMENT',self.env.metadataCloudCoverMax))\
		
		return s2s


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

		allCollection = ee.ImageCollection('COPERNICUS/S2').filterBounds(studyArea)
	                                           
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



	def TOAtoSR(self,img):
		
		TDOMMask = img.select(['TDOMMask'])
		info = self.collectionMeta[self.env.feature]['properties']
		scene_date = datetime.datetime.utcfromtimestamp(info['system:time_start']/1000)# i.e. Python uses seconds, EE uses milliseconds
		solar_z = info['MEAN_SOLAR_ZENITH_ANGLE']
        
		geom = ee.Geometry(info['system:footprint']).centroid()	
		#geom = ee.Geometry.Point([info['centroid']['coordinates'][0],info['centroid']['coordinates'][1]])
		date = ee.Date.fromYMD(scene_date.year,scene_date.month,scene_date.day)
		
		h2o = Atmospheric.water(geom,date).getInfo()
		o3 = Atmospheric.ozone(geom,date).getInfo()
		aot = Atmospheric.aerosol(geom,date).getInfo()

		SRTM = ee.Image('CGIAR/SRTM90_V4')# Shuttle Radar Topography mission covers *most* of the Earth
		alt = SRTM.reduceRegion(reducer = ee.Reducer.mean(),geometry = geom).get('elevation').getInfo()
		
		if alt:
			km = alt/1000 # i.e. Py6S uses units of kilometers
		
		else:
			km = 0
		# Instantiate
		s = SixS()
		
		# Atmospheric constituents
		s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o,o3)
		s.aero_profile = AeroProfile.Continental
		s.aot550 = aot
		
		# Earth-Sun-satellite geometry
		s.geometry = Geometry.User()
		s.geometry.view_z = 0               # always NADIR (I think..)
		s.geometry.solar_z = solar_z        # solar zenith angle
		s.geometry.month = scene_date.month # month and day used for Earth-Sun distance
		s.geometry.day = scene_date.day     # month and day used for Earth-Sun distance
		s.altitudes.set_sensor_satellite_level()
		s.altitudes.set_target_custom_altitude(km)
		
		def spectralResponseFunction(bandname):
			"""
			Extract spectral response function for given band name
			"""
            
			bandSelect = {
				'B1':PredefinedWavelengths.S2A_MSI_01,
				'B2':PredefinedWavelengths.S2A_MSI_02,
				'B3':PredefinedWavelengths.S2A_MSI_03,
				'B4':PredefinedWavelengths.S2A_MSI_04,
				'B5':PredefinedWavelengths.S2A_MSI_05,
				'B6':PredefinedWavelengths.S2A_MSI_06,
				'B7':PredefinedWavelengths.S2A_MSI_07,
				'B8':PredefinedWavelengths.S2A_MSI_08,
				'B8A':PredefinedWavelengths.S2A_MSI_09,
				'B9':PredefinedWavelengths.S2A_MSI_10,
				'B10':PredefinedWavelengths.S2A_MSI_11,
				'B11':PredefinedWavelengths.S2A_MSI_12,
				'B12':PredefinedWavelengths.S2A_MSI_13}
								    
			return Wavelength(bandSelect[bandname])

		def toa_to_rad(bandname):
			"""
			Converts top of atmosphere reflectance to at-sensor radiance
			"""
			
			# solar exoatmospheric spectral irradiance
			ESUN = info['SOLAR_IRRADIANCE_'+bandname]
			solar_angle_correction = math.cos(math.radians(solar_z))
			
			# Earth-Sun distance (from day of year)
			doy = scene_date.timetuple().tm_yday
			d = 1 - 0.01672 * math.cos(0.9856 * (doy-4))# http://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year
			
			# conversion factor
			multiplier = ESUN*solar_angle_correction/(math.pi*d**2)
			
			# at-sensor radiance
			rad = img.select(bandname).multiply(multiplier)
			return rad
			
		
		def surface_reflectance(bandname):
			"""
			Calculate surface reflectance from at-sensor radiance given waveband name
			"""
			
			# run 6S for this waveband
			s.wavelength = spectralResponseFunction(bandname)
			s.run()
			
			# extract 6S outputs
			Edir = s.outputs.direct_solar_irradiance             #direct solar irradiance
			Edif = s.outputs.diffuse_solar_irradiance            #diffuse solar irradiance
			Lp   = s.outputs.atmospheric_intrinsic_radiance      #path radiance
			absorb  = s.outputs.trans['global_gas'].upward       #absorption transmissivity
			scatter = s.outputs.trans['total_scattering'].upward #scattering transmissivity
			tau2 = absorb*scatter                                #total transmissivity
			
			# radiance to surface reflectance
			rad = toa_to_rad(bandname)
			
			ref = rad.subtract(Lp).multiply(math.pi).divide(tau2*(Edir+Edif))
			
			return ref
		
		# all wavebands
		output = img.select('QA60')
		for band in ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12']:
			output = output.addBands(surface_reflectance(band))
			
		self.env.feature +=1
		
		return output.addBands(TDOMMask).copyProperties(img,['system:time_start','system:footprint','MEAN_SOLAR_ZENITH_ANGLE','MEAN_SOLAR_AZIMUTH_ANGLE'])


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


	def medoidMosaic(self,collection):
		""" medoid composite with equal weight among indices """

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
  
		return medoid.addBands(others);


	def medianMosaic(self,collection):
		
		""" median composite """ 
		median = collection.select(medianIncludeBands).median();
		othersBands = bandNames.removeAll(medianIncludeBands);
		others = collection.select(otherBands).mean();
    
		return median.addBands(others)


	def setMetaData(self,img):
		""" add metadata to image """
		
		img = ee.Image(img).set({'system:time_start':ee.Date(self.env.startDate).millis(), \
								 'startDOY':str(self.env.startDoy), \
								 'endDOY':str(self.env.endDoy), \
								 'useCloudScore':str(self.env.cloudMask), \
								 'useTDOM':str(self.env.shadowMask), \
								 'useQAmask':str(self.env.QAcloudMask), \
								 'useCloudProject':str(self.env.cloudMask), \
								 'terrain':str(self.env.terrainCorrection), \
								 'surfaceReflectance':str(self.env.calcSR), \
								 'cloudScoreThresh':str(self.env.cloudScoreThresh), \
								 'cloudScorePctl':str(self.env.cloudScorePctl), \
								 'zScoreThresh':str(self.env.zScoreThresh), \
								 'shadowSumThresh':str(self.env.shadowSumThresh), \
								 'contractPixels':str(self.env.contractPixels), \
								 'crs':self.env.epsg, \
								 'cloudFilter':str(self.env.metadataCloudCoverMax),\
								 'dilatePixels':str(self.env.dilatePixels)})

		return img

	def exportMap(self,img,studyArea,week):

		geom  = studyArea.bounds().getInfo();
		
		task_ordered= ee.batch.Export.image.toAsset(image=img.clip(studyArea.buffer(10000)), 
								  description = self.env.name + str(week), 
								  assetId= self.env.assetId + self.env.name + str(week).zfill(3),
								  region=geom['coordinates'], 
								  maxPixels=1e13,
								  crs=self.env.epsg,
								  scale=self.env.exportScale)
	
		task_ordered.start() 
		

if __name__ == "__main__":        

	ee.Initialize()

	studyArea = ee.FeatureCollection("users/apoortinga/countries/Ecuador_nxprovincias").geometry() #.bounds();

	
	start = 1
	
	for i in range(413,468,1):
		startWeek = start+ i
		print startWeek
	
		year = ee.Date("2000-01-01")
		startDay = (startWeek -1) *14
		endDay = (startWeek) *14 -1
		
		startDate = year.advance(startDay,"day")
		endDate = year.advance(endDay,"day")		
		studyArea = ee.FeatureCollection("users/apoortinga/countries/Ecuador_nxprovincias").geometry().bounds();
		
		functions().main(studyArea,startDate,endDate,startDay,endDay,startWeek)
