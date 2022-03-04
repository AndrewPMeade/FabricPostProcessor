from .Tree import Tree

from .GenLib import NormVect
import numpy as np
import math

class MarkovTree(Tree):
	def __inti__(self, Root):
		Tree.__init__(self, Root)

	def CalcPMatrix(self, t, EigenVals, EigenVect, InvEigenVect):
		EigenVals = [math.exp(x*t) for x in EigenVals]
	
		P = (EigenVals * EigenVect)
		P = P.dot(InvEigenVect)
		return P

	def SetPMatrix(self, AMatrix):
		self.AMatrix = AMatrix

		(EigenVals, EigenVect) = np.linalg.eig(AMatrix)

		InvEigenVect = np.linalg.inv(EigenVect)

		for Node in self:
			if Node.Anc != None:
				Node.PMatrix = self.CalcPMatrix(Node.BranchLength, EigenVals, EigenVect, InvEigenVect)

	def GetNewSimState(self, PList):
		Point = np.random.random()
		CPoint = 0
		Index = 0
		for Val in PList:
			if Point >= CPoint and Point <= CPoint + Val:
				return Index
			CPoint += Val
			Index += 1
		return -1

	def BlankSim(self):

		for Node in self:
			Node.AncestralStates = []
	
	def PropSimStates(self, Node):

		AncState = Node.Anc.State
		Node.State = self.GetNewSimState(Node.PMatrix[AncState])

		for N in Node:
			self.PropSimStates(N)
	
	def SaveAncestralStates(self):
		
		for Node in self:
			Node.AncestralStates.append(Node.State)

	def SimData(self, NumberOfSites, AMatrix, RootStateFreq):

		self.BlankSim()
		self.SetPMatrix(AMatrix)

		RootStateFreq = NormVect(RootStateFreq);
	
		for i in range(NumberOfSites):
			self.Root.State = self.GetNewSimState(RootStateFreq)

			for N in self.Root:
				self.PropSimStates(N)

			self.SaveAncestralStates()