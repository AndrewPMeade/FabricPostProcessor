import sys
import numpy
import timeit

class PartFreq():
	def __init__(self, Node):
		self.BLList = []
		self.TaxaSet = Node.Partition
		
	def AddNode(self, Node):
		self.BLList.append(Node.BranchLength)

	def GetSize(self):
		return len(self.BLList)

	def SetBLInfo(self):

		self.Mean = numpy.mean(self.BLList)
		self.SD = numpy.std(self.BLList)
		self.Median = numpy.median(self.BLList)

	def SetProb(self, Prob):
		self.Prob = Prob 

	def Print(self, Trees):

		print(self.Prob, self.Mean, self.SD, self.Median, len(self.BLList), sep='\t', end='\t')

		Taxa = []
		for ID in self.TaxaSet:
			Taxa.append(Trees.IDTaxaMap[ID])
		print(*Taxa, sep=',', end='\t')

		print()

class ConsensusesTrees():

	def __init__(self, Trees, CutPoint=0.0):
				
		self.PartitionList = self.BuildPartitionListNew(Trees)
		self.ProcPartList(CutPoint, Trees.GetNoTrees())
#		self.PrintPartFreq(Trees)

	def BuildPartitionListNew(self, Trees):

		PartDic = {}

		for Tree in Trees:
			for Node in Tree:
				if Node.Partition in PartDic.keys():
					Part = PartDic[Node.Partition]
				else:
					Part = PartFreq(Node)
					PartDic[Node.Partition] = Part

				Part.AddNode(Node)

		return PartDic.values()
	
	def ProcPartList(self, CutPoint, NoTrees):

		for PFreq in self.PartitionList:
			PFreq.SetProb(PFreq.GetSize() / NoTrees)

		NList = [PFreq for PFreq in self.PartitionList if PFreq.Prob >= CutPoint]

		for PFreq in NList:
			PFreq.SetBLInfo()

		self.PartitionList = NList


	def PrintPartFreq(self, Trees):

		Header = ["Prob", "Mean BL", "SD BL", "Median BL", "Freq", "Taxa"]
		print(*Header, sep='\t')

		for PFreq in self.PartitionList:
			PFreq.Print(Trees)


#def CompConTreeTable(CTreeTableA, CTreeTableB):

