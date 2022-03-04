from .Tree import Tree
import sys, numpy, itertools, math



class GLSTree(Tree):
	def __init__(self, Root):
		Tree.__init__(self, Root)

	def GetDataVect(self):
		Ret = [0] * self.NoTaxa

		for Node in self:
			if Node.Tip == True:
				Ret[Node.TipID] = float(Node.Data[0])

		return Ret

	def CaclAlpha(self):

		P1 = 1.0 / self.InvV.sum()
		
		P2 = 0
		Index = 0
		for Row in self.InvV:
			RowSum = Row.sum()
			P2 += self.DataVect[Index] * RowSum 

			Index += 1

		Ret = P1 * P2

		return Ret

	def CaclZAlpha(self):
		return numpy.array([x - self.Alpha for x in self.DataVect])

	def CaclSigma2(self):

		Ret = self.ZAlpha.dot(self.InvV)
		Sig2 = Ret.dot(self.ZAlpha)
		
		Sig2 = Sig2 * (1.0 / self.NoTaxa)
		
		return Sig2

	def CaclKronecker(self):

		InvVSig = (1.0 / self.Sig2) * self.InvV

		return self.ZAlpha.dot(InvVSig)
		
	def CombineLh(self):

		Val = self.KroneckerVect.dot(self.ZAlpha) * -0.5

		Det = self.InvVDet + (self.NoTaxa * math.log(self.Sig2))		
		Det = -0.5 * Det 

		Ret = -self.NoTaxa / 2.0
		Ret = Ret * math.log(2*math.pi)

		return Ret + Det + Val

	def CalcLh(self):		

		# Find the number of ataxa
		self.NoTaxa = self.GetNoTaxa()

		# Convert the data to a vector
		self.DataVect = self.GetDataVect()

		# Get V
		self.V = self.TreeToV()

		# Invert V
		self.InvV = numpy.linalg.inv(self.V)
		
		# Get the log determinate 
		self.InvVDet = numpy.linalg.slogdet(self.V)
		self.InvVDet = self.InvVDet[1]

		# Calculate the mean 
		self.Alpha = self.CaclAlpha()
		
		# Calculate Z (Alpha - Data)
		self.ZAlpha = self.CaclZAlpha()

		# Calculate the variance  
		self.Sig2 = self.CaclSigma2()
		
		# Build the Kronecker vector
		self.KroneckerVect = self.CaclKronecker()
		
		# Combine the likelihood components
		self.Lh = self.CombineLh()

