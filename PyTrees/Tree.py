import sys
import numpy
import itertools
import math
import copy

class Tree():

	def __init__(self, Root):
		self.NodeList = []
		self.Root = Root

		self.RecMakeNodeList(Root)

	def __iter__(self):
		return iter(self.NodeList)

	def	__getitem__(self, Key):
		return self.NodeList[Key]

	def __len__(self):
		return len(self.NodeList)

	def next(self):
		next(self.NodeList)


	def RecMakeNodeList(self, Node):
		self.NodeList.append(Node)

		if Node.Tip == True:
			return

		for N in Node.NodeList:
			self.RecMakeNodeList(N)

	def RecSetPartitions(self, Node):

		for N in Node:
			self.RecSetPartitions(N)

		if Node.Tip == True:
			Node.Partition = frozenset([Node.TipID])
			return

		TID = set()
		for N in Node:
			TID.update(N.Partition)

		Node.Partition = frozenset(TID)

	def SetPartitions(self):
		self.RecSetPartitions(self.Root)

	def GetNodeFromTaxaID(self, TaxaID):

		for Node in self:
			if Node.TipID == TaxaID:
				return Node

		return None

	def GetNodeFromTaxa(self, Taxa):

		for Node in self:
			if Node.Taxa == Taxa:
				return Node

		return None


	def GetMRCA(self, TaxaSet):
		for FID in TaxaSet:
			break
		
		Node = self.GetNodeFromTaxaID(FID)
		while Node != None:

			if len(Node.Partition & TaxaSet) == len(TaxaSet):
				return Node

			Node = Node.Anc

		return None

	def RecSetDistToRoot(self, Node):

		Node.DistToRoot = Node.Anc.DistToRoot + Node.BranchLength
		for N in Node:
			self.RecSetDistToRoot(N)
			
	def SetDistToRoot(self):

		for Node in self.Root:
			self.RecSetDistToRoot(Node)

	def GetMaxDistToRoot(self):
		Ret = 0

		for Node in self:
			if Node.DistToRoot > Ret:
				Ret = Node.DistToRoot
		
		return Ret

	def SetHeight(self):
		MaxDistToRoot = self.GetMaxDistToRoot()

		for Node in self:
			Node.Height = MaxDistToRoot - Node.DistToRoot
	
	def GetTreeLength(self):
		TreeLenght = 0.0
		
		for Node in self:
			TreeLenght += Node.BranchLength
		
		return TreeLenght

	def GetNoTaxa(self):
		Ret = 0

		for N in self:
			if N.Tip == True:
				Ret += 1

		return Ret

	def __RecTreeToV(self, V, Node, Part):

		Dist = Node.Anc.DistToRoot

		PDiff = Part.difference(Node.Partition)
		
		for (i,j) in itertools.product(PDiff, Node.Partition):
			V[i,j] = Dist
			V[j,i] = Dist

		if Node.Tip == True:
			V[Node.TipID, Node.TipID] = Node.DistToRoot
			return

		for N in Node:
			self.__RecTreeToV(V, N, Node.Partition)



	def TreeToV(self):
		"""Returns the variance covariance matrix for a tree"""	
		NoTaxa = self.GetNoTaxa()

		Ret = numpy.zeros(shape=(NoTaxa,NoTaxa))
		
		for N in self.Root:
			self.__RecTreeToV(Ret, N, self.Root.Partition)

		return Ret

	def FindValidRoot(self):

		Node = self.Root
		while Node.GetNoDecedent() == 1:
			Node = Node.NodeList[0]

		self.Root = Node
		self.Root.BranchLength = 0
		self.Root.Rate = 1.0
		
	def RecRebuildTree(self, Node):

		if Node.Tip == True:
			return

		if Node.GetNoDecedent() == 0:
			Anc = Node.Anc
			Anc.RemoveDecedent(Node)
			return

		if Node.GetNoDecedent() == 1:
			Anc = Node.Anc
			Dec = Node.NodeList[0]

			Dec.BranchLength += Node.BranchLength

			Anc.RemoveDecedent(Node)

			Anc.AddNode(Dec)
			
			self.RecRebuildTree(Anc)
			return


		for N in Node:
			self.RecRebuildTree(N)

	def RebuildTree(self):

		self.RecRebuildTree(self.Root)
		self.FindValidRoot()

		self.NodeList = list()
		self.RecMakeNodeList(self.Root)

	
	def DeleteTaxa(self, TaxaList):

		for N in self.NodeList:
			if N.Taxa in TaxaList:
				N.Anc.RemoveDecedent(N)

		self.RebuildTree()

	def RecSimBMData(self, Anc, Sig2, Node):

		SimData = numpy.random.normal(scale=Sig2 * math.sqrt(Node.BranchLength))
		Node.Change = SimData
				
		SimData += Anc

		Node.SimData = SimData

		for N in Node:
			self.RecSimBMData(SimData, Sig2, N)


	def SimBMData(self, Sig2):

		Sig2 = math.sqrt(Sig2)

		for N in self.Root:
			self.RecSimBMData(0.0, Sig2, N)

	def SetTreeData(self, DataDic):
		
		for N in self:
			if N.Tip == True:
				N.Data = DataDic[N.Taxa]

	def Clone(self):
		return copy.deepcopy(self)

	def RemoveData(self):

		for Node in self:
			if hasattr(Node, "Data"):
				del Node.Data

	def LinagesAtHeight(self, Height):

		Ret = 0
		for Node in self:
			if Height >= Node.Height and Height < Node.Height + Node.BranchLength:
				Ret += 1

		return Ret

	def LTT(self, StepSize):

		Ret = []
		Height = 0;
		while Height <= self.Root.Height:
			Linages = self.LinagesAtHeight(Height)

			Ret.append((Height, Linages))

			Height += StepSize

		return Ret



