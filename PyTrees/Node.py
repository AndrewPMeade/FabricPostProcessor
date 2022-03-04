import sys
import hashlib
import copy

class Node:
	def __init__(self):
		self.NodeList = list()
		self.BranchLength = 0
		self.Anc = None
		self.Tip = False
		self.Taxa = ""
		self.TipID = -1
		#Partition frozen set of taxa id
		self.Partition = None
		self.DistToRoot = 0
		self.Height = -1
		self.Rate = 1
		self.Comment = None

	def __iter__(self):
		return iter(self.NodeList)

	def __len__(self):
		return len(self.NodeList)

	def next(self):
		next(self.NodeList)
	
	def __hash__(self):
		return hash(self.Partition)
	
	def __eq__(self, Other):
		if isinstance(Other, Node):
			return Other.__hash__() == self.__hash__()
		return False

	def __ne__(self, Other):
		return not self.__eq__(Other)

	def AddNode(self, Dec):
		Dec.Anc = self
		self.NodeList.append(Dec)

	def RecNodeToTipList(self, List):
		if self.Tip == True:
			List.append(self)
		
		for N in self:
			N.RecNodeToTipList(List)
	
	def NodeToTipList(self):
		"""return a list of tips from a node"""
		Ret = []

		self.RecNodeToTipList(Ret)

		return Ret

	def RemoveDecedent(self, Decedent):
		self.NodeList.remove(Decedent)

	def GetNoDecedent(self):
		return len(self.NodeList)
	
	def GetNodeHash(self, Trees):
		Hash = hashlib.blake2b()

		PList = list(self.Partition)
		PList.sort()

		for TID in PList:
			Taxa = Trees.IDTaxaMap[TID]
			Hash.update(Taxa.encode('utf-8'))

		return Hash.hexdigest()

	def RecPrintTaxa(self, sep='\t'):
		if self.Tip == True:
			print(self.Taxa, end=sep)
			return

		for N in self:
			N.RecPrintTaxa(sep)

	def RecScaleNode(self, Scale):

		self.BranchLength = self.BranchLength * Scale
		for N in self:
			N.RecScaleNode(Scale)

	def ScaleNode(self, Scale, ScaleSelf=False):

		if ScaleSelf == True:
			self.RecScaleNode(Scale)
		else:
			for N in self:
				N.RecScaleNode(Scale)
	
	def RecGetTaxaList(self, TList):

		if self.Tip == True:
			TList.append(self.Taxa)

		for Node in self:
			Node.RecGetTaxaList(TList)

	def GetTaxaList(self):
		TList = []

		self.RecGetTaxaList(TList)

		return TList

	def Clone(self):
		return copy.deepcopy(self)

	