

import ee
import math 
from utils import *

class env(object):
	
	def __init__(self):
		"""Initialize the environment"""
		self.studyArea = ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
		self.polar = ["VV","VH"]
		self.direction = 'ASCENDING'; 
		self.startDate = ee.Date.fromYMD(2017,1,1)
		self.endDate = ee.Date.fromYMD(2017,12,31)
		
		self.ksize = 3
		self.enl = 7; 
		

class functions():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
	
	    # get the environment
		self.env = env() 	

	def main(self):
		collection = self.getSAR()
		
		def applySpeckleFilter(img):
			vv = img.select('VV');
			vv = self.speckleFilter(img).rename(['VV']);
			
			if len(self.env.polar) > 1:
				vh = self.speckleFilter(img).rename(['VH']);
				return ee.Image.cat(vv,vh).copyProperties(img,['system:time_start']);
			
			else:
				return ee.Image.cat(vv).copyProperties(img,['system:time_start']);
				
				# Filter speckle using Gamma Map algorithm
		collection = collection.map(applySpeckleFilter)
		
			
		if len(self.env.polar) > 1:
			collection = collection.map(self.addRatio)
		
		print ee.Image(collection.first()).bandNames().getInfo()
		
		
		
	def getSAR(self):
		""" get the Sentinel 1 data collection"""
		# Inventory of S1 dual polarization images.
		collection = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(self.env.studyArea)\
														   .filter(ee.Filter.eq('instrumentMode','IW'))\
														   .filter(ee.Filter.eq('orbitProperties_pass', self.env.direction))\
														   .filter(ee.Filter.eq('transmitterReceiverPolarisation', self.env.polar))\
														   .filter(ee.Filter.date(self.env.startDate, self.env.endDate));
		
		return collection

	def speckleFilter(self,image):
		""" apply the speckle filter """
		# Convert image from dB to natural values
		nat_img = self.toNatural(image);

		# Square kernel, ksize should be odd (typically 3, 5 or 7)
		weights = ee.List.repeat(ee.List.repeat(1,self.env.ksize),self.env.ksize);

		# ~~(ksize/2) does integer division in JavaScript
		kernel = ee.Kernel.fixed(self.env.ksize,self.env.ksize, weights, ~~(self.env.ksize/2), ~~(self.env.ksize/2), False);

		# Get mean and variance
		mean = nat_img.reduceNeighborhood(ee.Reducer.mean(), kernel);
		variance = nat_img.reduceNeighborhood(ee.Reducer.variance(), kernel);

		# "Pure speckle" threshold
		ci = variance.sqrt().divide(mean);# square root of inverse of enl

		# If ci <= cu, the kernel lies in a "pure speckle" area -> return simple mean
		cu = 1.0/math.sqrt(self.env.enl);

		# If cu < ci < cmax the kernel lies in the low textured speckle area
		# -> return the filtered value
		cmax = math.sqrt(2.0) * cu;

		alpha = ee.Image(1.0 + cu*cu).divide(ci.multiply(ci).subtract(cu*cu));
		b = alpha.subtract(self.env.enl + 1.0);
		d = mean.multiply(mean).multiply(b).multiply(b).add(alpha.multiply(mean).multiply(nat_img).multiply(4.0*self.env.enl));
		f = b.multiply(mean).add(d.sqrt()).divide(alpha.multiply(2.0));

		# If ci > cmax do not filter at all (i.e. we don't do anything, other then masking)

		# Compose a 3 band image with the mean filtered "pure speckle", 
		# the "low textured" filtered and the unfiltered portions
		out = ee.Image.cat(self.toDB(mean.updateMask(ci.lte(cu))),self.toDB(f.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax))),image.updateMask(ci.gte(cmax)));	
		
		return out.reduce(ee.Reducer.sum());

	def addRatio(self,img):
		vv = self.toNatural(img.select('VV')).rename('VV');
		vh = self.toNatural(img.select('VH')).rename('VH');
		ratio = vh.divide(vv).rename('ratio');
		return ee.Image.cat(vv,vh,ratio).copyProperties(img,['system:time_start']);


	def toNatural(self,img):
		"""Function to convert from dB to natural"""
		return ee.Image(10.0).pow(img.select(0).divide(10.0));
		
	def toDB(self,img):
		""" Function to convert from natural to dB """
		return ee.Image(img).log10().multiply(10.0);
		
if __name__ == "__main__":
	ee.Initialize()
	functions().main()


