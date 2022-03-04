import sys
import PyTrees
import numpy
import scipy.stats
import math

def Chunk(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def LoadDataFile(FName):
	FIn = open(FName, "r")

	DataFile = FIn.readlines()

	DataFile = list(map(str.rstrip, DataFile))

	FIn.close()

	return DataFile


def LoadTaxa(Data):

	No = int(Data[0])
	MapData = Data[1:No+1]
	IDTaxaMap = {}

	for Line in MapData:
		Line = Line.split("\t")
		IDTaxaMap[int(Line[0])] = Line[1]

	return IDTaxaMap, Data[No+1:]

def MakeNode(NodeLine, MapData):

	Node = {}
	NodeLine = NodeLine.split("\t")

	Node["NodeID"] = int(NodeLine[0])
	Node["BranchLength"] = float(NodeLine[1])

	TaxaID = list(map(int, NodeLine[3:]))

	Node["Taxa"] = set([MapData[x] for x in TaxaID])

	return Node

def LoadNodes(Data, MapData):

	NodeList = []

	No = int(Data[0])
	NodeData = Data[1:No+1]

	for Line in NodeData:
		NodeList.append(MakeNode(Line, MapData))

	return NodeList, Data[No+1:]

def PassSampleNode(List):
	Ret = {}

	Ret["NodeID"] = int(List[0])
	Ret["Scalar"] = float(List[1])
	Ret["CrateIT"] = int(List[2])
	Ret["Type"] = List[3]

	return Ret

def PassSample(Line):

	Ret = {}

	Line = Line.split("\t")
	
	Ret["It"] = int(Line[0])
	Ret["Lh"] = float(Line[1])
	Ret["Lh+Prior"] =  float(Line[2])
	Ret["Alpha"] = float(Line[4])
	Ret["Sig2"] = float(Line[5])

	NodeList = list(Chunk(Line[7:], 4))
	Sample = []

	for N in NodeList:
		Sample.append(PassSampleNode(N))

	Ret["Sample"] = Sample

	return Ret

def LoadPosterior(DataFile):

	Sample = []
	for Line in DataFile:
		Sample.append(PassSample(Line))

	return Sample

#def GetNodeFromTaxaList(Trees, TaxaList):

def CreateNodeHash(NodeList, Trees):

	Ret = {}

	for Node in NodeList:
				
		NodeList = Trees.GetMRCA(Node["Taxa"])
		Ret[Node["NodeID"]] = NodeList[0]

	return Ret

def LoadData(FName):

	Data =  {}

	Data["FName"] = FName

	DataFile = LoadDataFile(FName)

	MapData, DataFile = LoadTaxa(DataFile)

	NodeList, DataFile = LoadNodes(DataFile, MapData)

	Sample = LoadPosterior(DataFile[1:])


	Data["IDTaxaMap"] = MapData
	Data["NodeList"] = NodeList
	Data["Sample"] = Sample
	Data["SampleSize"] = len(Sample)
	
	return Data

def InitTrees(Tree):

	for Node in Tree:
		Node.Beta = 0.0
		Node.BetaList = []

		Node.Scalar = 1.0
		Node.ScalarList = []

def ReSetTree(Tree):

	for Node in Tree:
		Node.Beta = 0.0
		Node.Scalar = 1.0



def GetNode(ID, NodeList):
		
	Ret = [N for N in NodeList if N["NodeID"] == ID]

	assert len(Ret) == 1

	return Ret[0]

def RecPropNode(Node, Scalar):
	Node.Scalar = Node.Scalar * Scalar

	for N in Node:
		RecPropNode(N, Scalar)

def SetNodeScalar(Node, Scalar, PropNodeScalar):
	if PropNodeScalar == False:
		Node.Scalar = Node.Scalar * Scalar
		return

	RecPropNode(Node, Scalar)

def RunSampleLine(Trees, Tree, Sample, Data, PropNodeScalar):

	for SNode in Sample:

		Node = Data["NodeHash"][SNode["NodeID"]]

		if SNode["Type"] == "Node":
			SetNodeScalar(Node, SNode["Scalar"], PropNodeScalar)

		if SNode["Type"] == "LS_Beta":
			Node.Beta = SNode["Scalar"]

		if SNode["Type"] == "Branch":
			Node.Scalar = Node.Scalar * SNode["Scalar"]


def SaveSample(Tree):

	for N in Tree:
		N.BetaList.append(N.Beta)
		N.ScalarList.append(N.Scalar)


def MapData(Trees, Data, PropNodeScalar):

	Tree = Trees.TreeList[0]

	InitTrees(Tree)

	for Sample in Data["Sample"]:
		ReSetTree(Tree)

		RunSampleLine(Trees, Tree, Sample["Sample"], Data, PropNodeScalar)

		SaveSample(Tree)

def ReSampleData(Data, SampleSize):
	SampleF = int(math.floor(len(Data["Sample"]) / SampleSize))

	Pos = 0

	Sample = Data["Sample"]

	ReSample = []

	while Pos < len(Sample):

		if Pos % SampleF == 0:
			ReSample.append(Sample[Pos])

		Pos += 1

	Data["Sample"] = ReSample[:SampleSize]

	Data["SampleSize"] = len(Data["Sample"])



def PassVarRates(TreeFN, VarRatesFN, PropNodeScalar=False, SampleSize=-1):

	Data = LoadData(VarRatesFN)

	if SampleSize != -1:
		ReSampleData(Data, SampleSize)

	Trees = PyTrees.Trees(TreeFN)

	Data["NodeHash"] = CreateNodeHash(Data['NodeList'], Trees)

	MapData(Trees, Data, PropNodeScalar)

	return Trees, Data

def PrintNodeInfo(Tree, Node):

	ID = Tree.NodeList.index(Node)
		
	print(ID, len(Node.Partition), Node.BranchLength,  sep='\t', end='\t')
	
	ScaleList = [x for x in Node.ScalarList if x != 1]
	NoScaled = len(ScaleList)
	print(NoScaled, NoScaled / len(Node.ScalarList), sep='\t', end='\t')

	Mean = numpy.mean(Node.ScalarList)
	
	Median = numpy.median(Node.ScalarList)
	SD = numpy.std(Node.ScalarList)
		
	print(Mean, Median, SD, sep='\t',end='\t')
	Node.RecPrintTaxa(sep=',')



def OutputVarRates(Trees, Data):
	Header = ["NodeID", "No Taxa", "BL", "No Scaled", "Pct Scaled", "Mean Scalar",  "Median Scalar", "SD Scalar", "Taxa"]

	print(*Header, sep='\t')

	Tree = Trees[0]

	for Node in Tree:
		PrintNodeInfo(Tree, Node)
		print()
		

if __name__ == "__main__":
	
	Trees, Data = PassVarRates(sys.argv[1], sys.argv[2], PropNodeScalar=True)
	OutputVarRates(Trees, Data)
