import ee

class diffMedian():
    def __init__(self):
        self.exportPath = 'users/TEST/'
        self.epsg = "EPSG:32717"
        self.ecoregions = ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Buffered")
        
        self.nDayBuffer = 7*10
        self.diffCountMin = 3
        self.biweeklyIC = 'projects/Sacha/PreprocessedData/L8_Biweekly_V6'


    def smartJoin(self,primary,secondary,julianDiff):
        """Function for joining based on max julian difference. Assumes startJulian and endJulian are set."""
        #Create a time filter to define a match as overlapping timestamps.
        maxDiffFilter = ee.Filter.Or(
            ee.Filter.maxDifference(
                difference = julianDiff,
                leftField = 'startJulian',
                rightField = 'endJulian'),
            ee.Filter.maxDifference(
                difference = julianDiff,
                leftField = 'startJulian',
                rightField = 'endJulian'
            )
        )

        # Define the join.
        saveAllJoin = ee.Join.saveAll(
            matchesKey = 'matches',
            measureKey = 'dayDiff'
        )

        #Apply the join.
        joined = saveAllJoin.apply(primary, secondary, maxDiffFilter)

        return joined

    def mergeMany(self,img,secondaryProperty,sortProperty):
        """Function to get the many secondaries and choose the first non null value"""
        img = ee.Image(img)
        secondaries = img.get(secondaryProperty)
        secondaries = ee.ImageCollection.fromImages(secondaries).sort(sortProperty)
        secondaries1 = secondaries.filter(ee.Filter.calendarRange(preYearPrimary,preYearPrimary,'year'))
        secondaries2 = secondaries.filter(ee.Filter.calendarRange(preYearSecondary,preYearSecondary,'year'))

        secondary1Composite = ee.Image(self.weightedCombiner(secondaries1))
        secondary2Composite = ee.Image(self.weightedCombiner(secondaries2))
        secondariesMosaiced = ee.ImageCollection([secondary1Composite,secondary2Composite]).mosaic()

        return img.addBands(secondariesMosaiced)

    def weightedCombiner(self,matches):
        """Function to take a set of matches and create a weighted median composite. Assumes the dayDiff property is set."""

        #Find the unique dayDiffs.
        matchesHist = ee.Dictionary(matches.aggregate_histogram('dayDiff'))

        #Convert it back to a number.
        def convertkeys(n):
            return ee.Number.parse(n).float()

        keys = matchesHist.keys().map(convertkeys)

        #Find the min and max of the dayDiffs and min max 0-1 stretch. Then reverse it and add 1 sdo the repeated values are from 1-20.
        minKey = keys.reduce(ee.Reducer.min())
        maxKey = keys.reduce(ee.Reducer.max())

        def normedn(n):
            return (ee.Number(n).subtract(minKey)).divide(ee.Number(minKey).add(maxKey))

        def normedreverse(n):
            return ee.Number(n).multiply(-1).add(1).multiply(20).int16()

        normedkeys = keys.map(normedn)
        normed = normedkeys.map(normedreverse)

        #Zip them together
        zipped = keys.zip(normed)

        def keyWeight(kw):
            keyWeight = ee.List(kw)
            key = keyWeight.get(0)
            weight = keyWeight.get(1)

            #Get images for given dayDiff.
            imgs = matches.filter(ee.Filter.eq('dayDiff',ee.Number(key)))

            def keyweightrepeat(img):
                return ee.ImageCollection(ee.List.repeat(ee.Image(img),ee.Number(weight)))

            #Repeat the images based on the weight.
            rep = ee.ImageCollection(ee.FeatureCollection(imgs.map(keyweightrepeat)).flatten())

            return rep

        repeated = zipped.map(keyWeight)

        #Flatten and compute median
        out = ee.ImageCollection(ee.FeatureCollection(repeated).flatten()).median()

        return ee.Image(out)

    def setJulian(self,img):
        """Function for setting start and end julians based on system:time_start. Assumes a 14 day diff inclusive of the first day."""
        d = ee.Date(img.get('system:time_start'))
        startJulian = d.getRelative('day','year')
        endJulian = startJulian.add(13)
        return img.set({'startJulian':startJulian,'endJulian':endJulian})

    def simpleAddIndices(self,in_image):
        """Function for only adding common indices."""
        in_image = in_image.addBands(in_image.normalizedDifference(['nir','red']).select([0],['NDVI']))
        in_image = in_image.addBands(in_image.normalizedDifference(['nir','swir2']).select([0], ['NBR']))
        in_image = in_image.addBands(in_image.normalizedDifference(['nir','swir1']).select([0], ['NDMI']))
        in_image = in_image.addBands(in_image.normalizedDifference(['green','swir1']).select([0], ['NDSI']))
        return in_image

    def cReducer(self,img):
        m = img.mask().reduce(ee.Reducer.min()).focal_min(3.5)
        return img.updateMask(m)

    def joinedmerge(self,img):
        return ee.Image(self.mergeMany(img,'matches','dayDiff'))
    
    def joinedmerge2(self,l):
        def joinedmerge(img):
            return ee.Image(self.mergeMany(img,'matches','dayDiff'))
        l = l.map(joinedmerge)
        return ee.ImageCollection(l)
    
    #Find the t2-t1 difference for each time period.
    def joineddiff(self,img):
        t1T = img.select(['.*_2014'])
        t2T = img.select(['.*_2016'])
        return img.addBands(t2T.subtract(t1T).rename(self.bnsDiff))
    
    def addsuffix(self,l,suffix):
        def base(bn):
            return ee.String(bn).cat(suffix)
        return l.map(base)
    
    def exportMap(self,img,studyArea):
        img = img
        ed = str(postYear)
        sd = str(preYearPrimary)
        regionName = ProcessingRegion.replace(" ",'_') + "_"

        task_ordered= ee.batch.Export.image.toAsset(image=img.clip(studyArea), 
                                  description = regionName + '_Diff_Comp_rSA_2lst_' + sd + '_' + ed, 
                                  assetId = self.exportPath + regionName + '_Diff_Comp' + sd + '_' + ed,
                                  region = studyArea.bounds().getInfo()['coordinates'], 
                                  maxPixels = 1e13,
                                  crs = self.epsg,
                                  scale = 30)

        task_ordered.start()
        print('Export Started: ',self.exportPath + regionName + '_Diff_Comp' + sd + '_' + ed)
   
    def main(self,ProcessingRegion,postYear,preYearPrimary,preYearSecondary,exportImg=False):
        studyArea = self.ecoregions.filter(ee.Filter.eq("PROVINCIA", ProcessingRegion)).geometry().buffer(1000)

        c = ee.ImageCollection(self.biweeklyIC).filter(ee.Filter.eq('regionName',ProcessingRegion)).map(self.cReducer).map(self.simpleAddIndices)

        bns = ee.List(['blue','green','red','nir','swir1','swir2','NDVI','NBR','NDMI'])

        c = c.select(bns)

        #Append endings to band names
        bnsT1 = self.addsuffix(bns,'_2014')
        bnsT2 = self.addsuffix(bns,'_2016')
        self.bnsDiff = self.addsuffix(bns,'_2014_2016_diff')

        #Filter off the two years of data

        t1 = c.filter(ee.Filter.calendarRange(preYearPrimary,preYearSecondary,'year')).select(bns,bnsT1).map(self.setJulian)
        t2 = c.filter(ee.Filter.calendarRange(postYear,postYear,'year')).select(bns,bnsT2).map(self.setJulian)
        print(t2.first().bandNames().getInfo())

        joined = ee.ImageCollection(self.smartJoin(t2,t1,self.nDayBuffer))

        joined = joined.toList(500)#.map(self.joinedmerge)
        joined = self.joinedmerge2(joined)
        print(joined.first().bandNames().getInfo(),'join')

        diff = joined.map(self.joineddiff)

        diffMedian = diff.median()
        diffCount = diff.select(['.*_diff']).count().reduce(ee.Reducer.min())
        diffMedian = diffMedian.updateMask(diffCount.gte(self.diffCountMin))

        
        print(diffMedian.getInfo())
        if exportImg:
            self.exportMap(diffMedian,studyArea)

        return diffMedian

if __name__ == "__main__":
    ee.Initialize()

    ProcessingRegion = 'GALAPAGOS'
    
    postYear = 2016
    preYearPrimary = 2014
    preYearSecondary = 2013

    exportImg = True

    diffMedian().main(ProcessingRegion,postYear,preYearPrimary,preYearSecondary,exportImg)
