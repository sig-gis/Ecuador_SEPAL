import ee
import math

ee.Initialize()
class addCovariates():
  def __init__(self):
    self.jrcImage = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
    self.elevation = ee.Image("USGS/SRTMGL1_003");
      
    self.ndCovariatesList = [['blue', 'green'],['blue', 'red'],['blue', 'nir'],['blue', 'swir1'],['blue', 'swir2'],
      ['green', 'red'],['green', 'nir'],['green', 'swir1'],['green', 'swir2'],['red', 'swir1'],['red', 'swir2'],
      ['nir', 'red'],['nir', 'swir1'],['nir', 'swir2'],['swir1', 'swir2']]


  rCovariatesList = [['swir1', 'nir'],['red', 'swir1']];

  def ComputeNDCovariatesList(self,season):
    l = [];
    for index in self.ndCovariatesList:
      currentIndex = self.ndCovariatesList.index(index)
      list_ = [season + self.ndCovariatesList[currentIndex][0], season + self.ndCovariatesList[currentIndex][1]];
      l.append(list_);

    return l;


  def addNDCovariates(self, season, image):
    l = self.ComputeNDCovariatesList(season);
    for index in l:
      currentIndex = l.index(index)
      image = image.addBands(image.normalizedDifference(index).rename(season + 'ND_'+ self.ndCovariatesList[currentIndex][0] + '_' + self.ndCovariatesList[currentIndex][1]));

    return image;


  def ComputeRCovariatesList(season):
    l = [];
    for index in rCovariatesList:
      list_ = [season  + rCovariatesList[index][0], season  + rCovariatesList[index][1]];
      l.append(list_);

    return l;


  def addRCovariates(season, image):
    l = ComputeRCovariatesList(season);
    for index in l:
      image = image.addBands(image.select(l[index][0]).divide(image.select(l[index][1]))
              .rename(season + '_R_' + rCovariatesList[index][0] + '_' + rCovariatesList[index][1]));

    return image;


  def addEVI(season, image):
    evi = image.expression('2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR' : image.select(season + '_nir'),
      'RED' : image.select(season + '_red'),
      'BLUE': image.select(season + '_blue')}).float();
    return image.addBands(evi.rename(season + '_EVI'));


  def addSAVI(season, image):
    #// Add Soil Adjust Vegetation Index (SAVI)
  	#// using L = 0.5;
  	savi = image.expression('(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
      'NIR': image.select(season + '_nir'),
      'RED': image.select(season + '_red')}).float();
  	return image.addBands(savi.rename(season + '_SAVI'));


  def addIBI(season, image):
    # // Add Index-Based Built-Up Index (IBI)
    ibiA = image.expression('2 * SWIR1 / (SWIR1 + NIR)', {
      'SWIR1': image.select(season + '_swir1'),
      'NIR'  : image.select(season + '_nir')
    }).rename(['IBI_A']);

    ibiB = image.expression('(NIR / (NIR + RED)) + (GREEN / (GREEN + SWIR1))', {
      'NIR'  : image.select(season + '_nir'),
      'RED'  : image.select(season + '_red'),
      'GREEN': image.select(season + '_green'),
      'SWIR1': image.select(season + '_swir1')
    }).rename(['IBI_B']);

    ibiAB = ibiA.addBands(ibiB);
    ibi = ibiAB.normalizedDifference(['IBI_A', 'IBI_B']);
    return image.addBands(ibi.rename([season + '_IBI']));


  # // Function to compute the Tasseled Cap transformation and return an image
  def getTassledCapComponents(season, image):

    coefficients = ee.Array([
      [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
      [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
      [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
      [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
      [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
      [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]]);
    bands = ee.List([season + '_blue', season + '_green', season + '_red', season + '_nir', season + '_swir1', season + '_swir2']);

    # // Make an Array Image, with a 1-D Array per pixel.
    arrayImage1D = image.select(bands).toArray();

    # // Make an Array Image with a 2-D Array per pixel, 6 x 1
    arrayImage2D = arrayImage1D.toArray(1);

    componentsImage = ee.Image(coefficients).matrixMultiply(arrayImage2D).arrayProject([0]).arrayFlatten([[season + '_brightness', season + '_greenness', season + '_wetness', season + '_fourth', season + '_fifth', season + '_sixth']]).float();
    # // Get a multi-band image with TC-named bands
    return image.addBands(componentsImage);


  # // Function to add Tasseled Cap angles and distances to an image. Assumes image has bands: 'brightness', 'greenness', and 'wetness'.
  def getTassledCapAngleAndDistance (season, image):

    brightness = image.select(season + '_brightness');
    greenness = image.select(season + '_greenness');
    wetness = image.select(season + '_wetness');

    # // Calculate tassled cap angles and distances
    tcAngleBG = brightness.atan2(greenness).divide(Math.PI).rename([season + '_tcAngleBG']);
    tcAngleGW = greenness.atan2(wetness).divide(Math.PI).rename([season + '_tcAngleGW']);
    tcAngleBW = brightness.atan2(wetness).divide(Math.PI).rename([season + '_tcAngleBW']);

    tcDistanceBG = brightness.hypot(greenness).rename([season + '_tcDistanceBG']);
    tcDistanceGW = greenness.hypot(wetness).rename([season + '_tcDistanceGW']);
    tcDistanceBW = brightness.hypot(wetness).rename([season + '_tcDistanceBW']);

    image = image.addBands(tcAngleBG).addBands(tcAngleGW).addBands(tcAngleBW).addBands(tcDistanceBG).addBands(tcDistanceGW).addBands(tcDistanceBW);

    return image;


  def computeTassledCap(season, image):
    image = getTassledCapComponents(season, image);
    image = getTassledCapAngleAndDistance(season, image);
    return image;


  def addTopography (self,image):

    # // Calculate slope, aspect and hillshade
    topo = ee.Algorithms.Terrain(self.elevation);

    # // From aspect (a), calculate eastness (sin a), northness (cos a)
    deg2rad = ee.Number(math.pi).divide(180);
    aspect = topo.select(['aspect']);
    aspect_rad = aspect.multiply(deg2rad);
    eastness = aspect_rad.sin().rename(['eastness']).float();
    northness = aspect_rad.cos().rename(['northness']).float();

    # // Add topography bands to image
    topo = topo.select(['elevation','slope','aspect']).addBands(eastness).addBands(northness);
    image = image.addBands(topo);
    return image;


  def addJRCDataset (self,image):
    # // Update the mask.
    jrcImage = self.jrcImage.unmask(0);

    image = image.addBands(jrcImage.select(['occurrence']).rename(['occurrence']));
    image = image.addBands(jrcImage.select(['change_abs']).rename(['change_abs']));
    image = image.addBands(jrcImage.select(['change_norm']).rename(['change_norm']));
    image = image.addBands(jrcImage.select(['seasonality']).rename(['seasonality']));
    image = image.addBands(jrcImage.select(['transition']).rename(['transition']));
    image = image.addBands(jrcImage.select(['max_extent']).rename(['max_extent']));

    return image;

  def addCovariates (self,season, image):
    image = self.addNDCovariates(season, image);
    # /*image = addRCovariates(season, image);
    # image = addEVI(season, image);
    # image = addSAVI(season, image);
    # image = addIBI(season, image);
    # image = computeTassledCap(season, image);*/
    return image;


  def addJRCAndTopo(self,image):
    image = self.addTopography(image);
    image = self.addJRCDataset(image);
    return image;








