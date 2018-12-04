

class assemblage():

	def __init__(self):
		pass
		
	def createAssemblage(self,image,nodeStruct):

		names = image.bandNames().getInfo()
		classes = ['other']

		for name in names:
			classes.append(name.encode('ascii'))
		numbers = range(0,image.bandNames().length().getInfo()+1,1)
		
		classStruct = {}
		for i in numbers:
			classStruct[classes[i]] = {'number' : i}
		
		# The starting id, i.e. the first decision
		startId = 'key1';
		
		# The initial decision tree string (DO NOT CHANGE)
		DTstring = ['1) root 9999 9999 9999'];
		# Call the function to construct the decision tree string (DO NOT CHANGE)
		DTstring = "\n".join(self.decision(nodeStruct,classStruct,startId,1,DTstring))#.join("\n");

		classifier = ee.Classifier.decisionTree(DTstring)
		
		sd = 10
		nIter = 100
		nBands = ee.List.sequence(0,image.bandNames().length().getInfo()-1,1)
		
	
		def monteCarlo(i):
			def createRand(j):
				rand = ee.Image(image.select(0)).gt(-1).multiply(ee.Image.random(ee.Number(i).multiply(j))).subtract(0.5).multiply(2)
				random = rand.multiply(sd)
				return ee.Image(image.select(ee.Number(j))).add(random)
			img = ee.ImageCollection(nBands.map(createRand))
			img = self.collectionToImage(img)
			classified = img.classify(classifier);
			return classified
		
		
		iters = ee.List.sequence(1,nIter)
		assemblage = ee.ImageCollection(iters.map(monteCarlo))
				
		mode = ee.Image(assemblage.mode())
		
		def uncertainty(img):
			return img.eq(mode)
		
		prob = ee.Image(assemblage.map(uncertainty).sum()).rename('prob')
		
		#region =  ee.Geometry.Polygon([[103.876,18.552],[105.806,18.552],[105.806,19.999],[103.876,19.999],[103.876,18.552]])

		#task_ordered = ee.batch.Export.image.toAsset(image=ee.Image(prob), description='test', assetId='users/servirmekong/temp/outputAssemblage4',region=region['coordinates'], maxPixels=1e13,scale=300)
		#task_ordered.start() 

		#print(classified.getInfo())
		
		return mode, prob

	# Function to convert a dictionary of nodes into a decision tree string
	def decision(self,nodeStruct,classStruct,id1,node,DTstring):
		# Extract parameters
		lnode = 2*node; #// left child node number
		rnode = lnode + 1; #// right child node number
		dict1 = nodeStruct[id1]; #// current dictionary
		band = dict1['band']; #// decision band
		threshold = dict1['threshold']; #// decision threshold
		left = dict1['left']; #// left result (either 'terminal' or new id)
		right = dict1['right']; #// right result (either 'terminal' or new id)
		leftName = dict1['leftName']; #// left class name (if 'terminal')
		if right == 'terminal':
			rightName = dict1['rightName']; #// right class name (if 'terminal')
		
		leftLine = ''
		rightLine = '';
		leftNumber = 0;
		rightNumber = 0;

		# Add the left condition and right condition strings to the current decision 
		# tree string. If either condition is non-terminal, recursively call the function.
		if (left == 'terminal'): # left terminal condition
			leftNumber = str(classStruct[leftName]['number']);
			leftLine = str(lnode) + ') ' + str(band) + '>=' + str(threshold) + ' 9999 9999 ' + str(leftNumber) + ' *';
			DTstring.append(leftLine);
			if (right == 'terminal'): #// right terminal condition
				rightNumber = classStruct[rightName]['number'];
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 ' + str(rightNumber) + ' *';
				DTstring.append(rightLine);
				return DTstring;
			else: # right non-terminal condition
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 9999';
				DTstring.append(rightLine);
				return self.decision(nodeStruct,classStruct,right,rnode,DTstring);
    
		else: # left non-terminal condition
			leftLine = str(lnode) + ') ' + str(band) + '>=' + str(threshold) + ' 9999 9999 9999';
			DTstring.append(leftLine);	
			DTstring = self.decision(nodeStruct,classStruct,left,lnode,DTstring);
		
			if (right == 'terminal'): # // right terminal condition
				rightNumber = classStruct[rightName]['number'];
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 ' + str(rightNumber) + ' *';
				DTstring.append(rightLine);
				return DTstring;
			else:  # right non-terminal
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 9999';
				DTstring.append(rightLine);
				return self.decision(nodeStruct,classStruct,right,rnode,DTstring);
    
  
		return DTstring;

	def collectionToImage(self,collection):
		
		def iterate(img,prev):
			return ee.Image(prev).addBands(img)
		
		stack = ee.Image(collection.iterate(iterate,ee.Image(1)))
		 
		stack = stack.select(ee.List.sequence(1, stack.bandNames().size().subtract(1)));
		
		return stack;

		
import ee
ee.Initialize()


aquaculture = ee.Image(ee.ImageCollection("projects/servir-mekong/yearly_primitives_smoothed/aquaculture").first()).rename('aquaculture')
barren = ee.Image(ee.ImageCollection("projects/servir-mekong/yearly_primitives_smoothed/barren").first()).rename('barren')
cropland = ee.Image(ee.ImageCollection("projects/servir-mekong/yearly_primitives_smoothed/cropland").first()).rename('cropland')
deciduous = ee.Image(ee.ImageCollection("projects/servir-mekong/yearly_primitives_smoothed/deciduous").first()).rename('forest')
		
image = aquaculture.addBands(barren).addBands(cropland).addBands(deciduous)

nodeStruct = { 	'key1':  {'band': 'aquaculture','threshold': 50, 'left': 'terminal', 'leftName': 'aquaculture', 'right': 'key2'},
				'key2':  {'band': 'barren', 'threshold': 40, 'left': 'terminal', 'leftName': 'barren', 'right': 'key3'},
				'key3':  {'band': 'cropland', 'threshold': 60, 'left': 'terminal', 'leftName': 'cropland', 'right': 'key4'},
				'key4':  {'band': 'forest', 'threshold': 5, 'left': 'terminal', 'leftName': 'other', 'right': 'terminal', 'rightName': 'forest'}	};

m,p = assemblage().createAssemblage(image,nodeStruct)
print(m.getInfo())
print(p.getInfo())
