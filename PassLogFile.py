import sys
import re
import numpy as np

from Distribution import CreateDistribution

LAND_PRIOR_STR = "\s+(VR_LS_BL|FabricBeta) - (.+)"
LAND_THRESHOLD_STR = "\s*RJ Local LandscapeBL:\s+True Threshold (.+)"



NODE_PRIOR_STR = "\s+VRNode - (.+)"
NODE_THRESHOLD_STR = "\s*RJ Local Node:\s+True Threshold(.+)"

HEADER_STR = "Iteration	Lh"


class PassLogFile():

	def __init__(self, FileName):

		self.FlieName = FileName 
		
		self.LandThreshold = 0
		self.LandPrior = CreateDistribution("uniform", [-100, 100])
		
		self.NodePrior = CreateDistribution("uniform", [0, 100])
		self.NodeThreshold = 0

		
		self.__pass_file(FileName)
##		self.__check_info()
		self.__post_process()

	def __post_process(self):

		lh1 = self.Lh[:-1]
		lh2 = self.Lh[1:]


		correl = np.corrcoef(lh1, lh2)

		self.LhR2 = correl[0][1]**2
		self.ESS = len(self.Lh) * (1.0-self.LhR2)
		
		return 

		self.Prior.Paramters = [0.0, 1.0]

		for x in np.arange(0.001, 2.5, 0.01):
			print(x, self.Prior.CDF(x), sep='\t')


		
	def __check_info(self):

		if hasattr(self, "LandThreshold") == False:
			raise Exception(f"Landscape Threshold could not be found in log flie {self.FlieName}")
			
		if hasattr(self, "LandPrior") == False:
			raise Exception(f"Landscape Prior could not be found in log flie {self.FlieName}")

		if hasattr(self, "NodeThreshold") == False:
			raise Exception(f"Node Threshold could not be found in log flie {self.FlieName}")

		if hasattr(self, "NodePrior") == False:
			raise Exception(f"Node Prior could not be found in log flie {self.FlieName}")

		if hasattr(self, "Lh") == False:
			raise Exception(f"Lh could not be found in log flie {self.FlieName}")


	def __pass_file(self, FileName):

		with open(FileName, "r") as FIn:
			for Line in FIn:

				Match = re.match(LAND_PRIOR_STR, Line)
				if Match:
					self.LandPrior = self.__pass_prior(Match[2])

				Match = re.match(LAND_THRESHOLD_STR, Line)
				if Match:
					self.LandThreshold = float(Match[1])

				Match = re.match(NODE_PRIOR_STR, Line)
				if Match:
					self.NodePrior = self.__pass_prior(Match[1])

				Match = re.match(NODE_THRESHOLD_STR, Line)
				if Match:
					self.NodeThreshold = float(Match[1])


				Match = re.match(HEADER_STR, Line)
				if Match:
					self.Lh = self.__pass_lh(FIn)
									
	def __pass_lh(self, FIn):

		Lh = []
		for Line in FIn:
			Line = Line.split("\t")

			Lh.append(float(Line[1]))

		return Lh

	def __pass_prior(self, PriorStr):
		Prior = {}
		PriorStr = PriorStr.rstrip()
		Line = PriorStr.split(" ")

		return CreateDistribution(Line[0], [float(x) for x in Line[1:]])
		