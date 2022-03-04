import sys
import re
import copy

from .TreeFormat import TreeFormat
from .Node import Node
from .StringToTree import StringToTree
from .Tree import Tree
from .PassNexus import PassNexus

class Trees:
	def __init__(self, FName, Format=TreeFormat.Nexus, TType=Tree, NType=Node):
		""" Load a sample of trees in nexus for phylip format.
			Takes:
			A file name names
			Options 
				Format: TreeFormat.Nexus (default) or TreeFormat.Phylip
				TType: Tree to use, must be inherited from type Tree (default)
				NType: Node to use, must be inherited from type Node (default) 
	    """
		# input file name
		self.FName = FName

		# list of tree
		self.TreeList = list()

		# hash to map taxa Name to ID
		self.TaxaIDMap = {}

		# has to map ID to taxa Names
		self.IDTaxaMap = {}
		
		if Format == TreeFormat.Phylip:
			self.LoadPhylip(TType, NType)

		if Format == TreeFormat.Nexus:
			self.LoadNexus(TType, NType)

		self.InputFormat = Format

		self.PostProcesses()

		self.Data = {}
		

	def __iter__(self):
		return iter(self.TreeList)

	def	__getitem__(self, Key):
		return self.TreeList[Key]

	def next(self) -> Tree:
		next(self.TreeList)

	def __len__(self):
		return len(self.TreeList)

	def LoadPhylip(self, TType, NType):
		
		FIn = open(self.FName, "r")

		for Line in FIn:
			Line = Line.rstrip()
			Root = StringToTree(Line, NType)
			self.TreeList.append(TType(Root))

		FIn.close()
		self.SetTreeTaxaID()


	def LoadNexus(self, TType, NType):

		Nexus = PassNexus(self.FName, NType)
		for Root in Nexus.RootList:
			self.TreeList.append(TType(Root))
		Nexus = None

	def MakeTaxaNameList(self):
		""" Return a list of taxa names """
		TSet = set()

		for T in self.TreeList:
			for N in T.NodeList:
				if N.Tip == True:
					TSet.add(N.Taxa)
		return list(TSet)
	

	def BuildTaxaIDMaps(self):
		""" 
		Build the dictionaries that map taxa id to names and reverse. 
		For nexus files the taxa numbers are reset, starting from 0. 
		"""
		TaxaList = self.MakeTaxaNameList()
		TaxaList.sort()

		self.TaxaIDMap = dict()
		self.IDTaxaMap = dict()

		for Taxa in TaxaList:
			self.TaxaIDMap[Taxa] = TaxaList.index(Taxa)
			self.IDTaxaMap[TaxaList.index(Taxa)] = Taxa

	def SetTreeTaxaID(self):
		self.BuildTaxaIDMaps()

		for Tree in self.TreeList:
			for Node in Tree:
				if Node.Tip == True:
					Node.TipID = self.TaxaIDMap[Node.Taxa]

	def SetPartitions(self):
		for Tree in self.TreeList:
			Tree.SetPartitions()
	
	def SetDistToRoot(self):
		for Tree in self.TreeList:
			Tree.SetDistToRoot()

	def SetHeight(self):
		for Tree in self.TreeList:
			Tree.SetHeight()


	def PostProcesses(self):
		self.SetTreeTaxaID()
		self.SetPartitions()
		self.SetDistToRoot()
		self.SetHeight()
		
	def TaxaListToIDSet(self, TaxaList):

		Ret = []

		for TaxaName in TaxaList:
			Ret.append(self.TaxaIDMap[TaxaName])

		return set(Ret)

	# given a list of taxa, find the MRCA in every tree
	def GetMRCA(self, TaxaList):
		
		Ret = []
		MRCASet = self.TaxaListToIDSet(TaxaList)

		for Tree in self.TreeList:
			Ret.append(Tree.GetMRCA(MRCASet))

		return Ret

	def GetTreeMRCA(self, TreeNo, TaxaList):
		MRCASet = self.TaxaListToIDSet(TaxaList)
		Tree = self.TreeList[TreeNo]

		return Tree.GetMRCA(MRCASet)
	
	def NodeToTaxaList(self, Node):

		Ret = []

		for PID in Node.Partition:
			Ret.append(self.IDTaxaMap[PID])

		return Ret

	def PrintNode(self, Node):
		
		TList = self.NodeToTaxaList(Node)
		for T in TList:
			print(T, end='\t')


	def SaveBL(self, Node, FOut, SaveRates):

		print(":%8.8f" % Node.BranchLength, end='', file=FOut)

		if SaveRates == False:
			return
		
		print("@%f" % Node.Rate, end='', file=FOut)

	def SaveNexusTree(self, Node, FOut, SaveRates):
		if Node.Tip  == True:
			print("%d" % Node.TipID, end='', file=FOut)
			self.SaveBL(Node, FOut, SaveRates)
			return

		print("(", end='', file=FOut)

		for N in Node.NodeList:
			self.SaveNexusTree(N, FOut, SaveRates)
			if N != Node.NodeList[-1]:
				print(",", end='', file=FOut)

		if Node.Anc != None:
			print(")", end='', file=FOut)
			self.SaveBL(Node, FOut, SaveRates)
		else:
			print(")", end='', file=FOut)

	def SaveNexusHeader(self, OutFile):
		
		print("#NEXUS", file=OutFile)
		print("begin trees;", file=OutFile)
		print("\ttranslate", file=OutFile)

		i = 0
		Size = len(self.IDTaxaMap.items())
		for (K, V) in sorted (self.IDTaxaMap.items()):
			if i < Size-1:
				print("\t\t%d %s," % (K, V), file=OutFile)
			else:
				print("\t\t%d %s;" % (K, V), file=OutFile)

			i += 1
			
	def SaveNexusFooter(self, OutFile):
		print("end;\n", file=OutFile)

	def SaveNexus(self, FName, SaveRates=False):

		FOut = open(FName, "w")

		self.SaveNexusHeader(FOut)

		i = 0
		for Tree in self:
			print("\t\ttree Tree_%08d = " % i, end='', file=FOut)
			self.SaveNexusTree(Tree.Root, FOut, SaveRates)
			print(";", file=FOut)
			i += 1
		print("end;", file=FOut)	
		FOut.close()
	
	def CheckData(self, FName):

		for (Key, Val) in self.TaxaIDMap.items():
			if Key not in self.Data:
				raise Exception("Taxa (%s) not found in Data file %s" % (Key, FName))

	def SetTreeData(self):

		for Tree in self:
			Tree.SetTreeData(self.Data)

	def LoadData(self, FName):
		""" Load data into Data Dictionary, file format should be BayesTraits """
		FIn = open(FName, "r")

		for Line in FIn:
			Line = Line.rstrip()
			Line = re.split("[\t ]+", Line)
			
			if Line[0] in self.TaxaIDMap:
				self.Data[Line[0]] = Line[1:]

		self.CheckData(FName)

		self.SetTreeData()

		FIn.close()

	def FreqDicToList(self, FreqDic):

		Max = 0
		for (K,V) in FreqDic.items():
			print(K, Max)
			if K > Max:
				Max = K

		Ret = [0] * (Max+1)
		for (K,V) in FreqDic.items():
			Ret[K] = V

		return Ret


	def PolytomyDistribution(self):
		""" Return A list of the frequency of the polytomy  """
		
		PolyDist = {}
		for Tree in self.TreeList:
			for N in Tree:
				if N.Tip == False:
					NoP = len(N.NodeList)
					if NoP in PolyDist:
						PolyDist[NoP] = PolyDist[NoP] + 1
					else:
						PolyDist[NoP] = 1
		
		return self.FreqDicToList(PolyDist)

	
	def DeleteTaxa(self, TaxaList):
		"""Delete a set of taxa from the tree, takes a list of taxa names."""
		for Tree in self:
			Tree.DeleteTaxa(TaxaList)

		self.PostProcesses()
						

	def CollapseNode(self, NewName, TaxaList):
		""" Convert an internal node to a tip with a given name. 
		The node specified must be monophyletic """
		TaxaIDSet = self.TaxaListToIDSet(TaxaList)

		for Tree in self:
			MRCA = Tree.GetMRCA(TaxaIDSet)
			
			if len(MRCA.Partition) != len(TaxaIDSet):
				 raise Exception("Collapsed node must be monophyletic")

			MRCA.Tip = True
			MRCA.Taxa = NewName

			Tree.RebuildTree()

		self.PostProcesses()


	def GetNoTrees(self):
		""" return the number of trees """
		return len(self.TreeList)

	def Clone(self):
		return copy.deepcopy(self)