class Switcher(object):
	def paramSelect(self, argument):
		"""Dispatch method"""
		method_name = str(argument).replace(" ","_")
		# Get the method from 'self'. Default to a lambda.
		method = getattr(self, method_name, lambda: "nothing")

		return method()

	def AMAZONIA_NOROCCIDENTAL(self):
		cloudScoreThresh= 1
		cloudScorePctl= 8
		zScoreThresh= -0.9
		shadowSumThresh= 0.4
		contractPixels= 1.5
		dilatePixels= 3.25

		return cloudScoreThresh, cloudScorePctl, zScoreThresh, shadowSumThresh, contractPixels, dilatePixels

	def ANDES_DEL_NORTE(self):
		cloudScoreThresh= 11
		cloudScorePctl= 8
		zScoreThresh= -0.8
		shadowSumThresh= 0.15
		contractPixels= 1.25
		dilatePixels= 2.75

		return cloudScoreThresh, cloudScorePctl, zScoreThresh, shadowSumThresh, contractPixels, dilatePixels

	def CHOCO(self):
		cloudScoreThresh= 1
		cloudScorePctl= 8
		zScoreThresh= -0.9
		shadowSumThresh= 0.35
		contractPixels= 0.86
		dilatePixels= 2

		return cloudScoreThresh, cloudScorePctl, zScoreThresh, shadowSumThresh, contractPixels, dilatePixels

	def GALAPAGOS(self):
		cloudScoreThresh= 1
		cloudScorePctl= 8
		zScoreThresh= -1.11
		shadowSumThresh= 0.25
		contractPixels= 1.18
		dilatePixels= 2

		return cloudScoreThresh, cloudScorePctl, zScoreThresh, shadowSumThresh, contractPixels, dilatePixels

	def PACIFICO_ECUATORIAL(self):
		cloudScoreThresh= 1
		cloudScorePctl= 8
		zScoreThresh= -0.9
		shadowSumThresh= 0.35
		contractPixels= 0.86
		dilatePixels= 2

		return cloudScoreThresh, cloudScorePctl, zScoreThresh, shadowSumThresh, contractPixels, dilatePixels

	def SIERRA(self):
		cloudScoreThresh= 2
		cloudScorePctl= 11
		zScoreThresh= -0.8
		shadowSumThresh= 0.15
		contractPixels= 1.2
		dilatePixels= 2.5

		return cloudScoreThresh, cloudScorePctl, zScoreThresh, shadowSumThresh, contractPixels, dilatePixels
