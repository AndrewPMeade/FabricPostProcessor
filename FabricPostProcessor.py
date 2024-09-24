import sys
import PyTrees
import numpy
import hashlib
import math
import random
import pandas as pd

from PassVarRates import PassVarRates
from PassLogFile import PassLogFile
from Distribution import CreateDistribution
from Utilities import RecGetTaxa, GetMD5Sum

JOIN_TAXA_NAMES = False
Z_SIGNIFICANCE  = 2.0
Z_OVERFLOW_VAL = 100.0

BETA_Z_PLACE_HOLDER = "BETA_Z_PLACE_HOLDER"
SCALAR_Z_PLACE_HOLDER = "SCALAR_Z_PLACE_HOLDER"

BETA_SIG_PLACE_HOLDER = "BETA_SIG_PLACE_HOLDER"
SCALAR_SIG_PLACE_HOLDER = "SCALAR_SIG_PLACE_HOLDER"


def PostProcess(Trees):
	Tree = Trees[0]

	for Node in Tree:
		NonZero = [x for x in Node.BetaList if x != 0]
		Node.BetaPP = len(NonZero) / len(Node.BetaList)
		Node.BetaListBin = [0 if x == 0 else 1 for x in Node.BetaList]


def CaclSDP(ObsP, ESS):

	if ObsP == 1:
		return 0

	Ret = ObsP * (1.0 - ObsP)
	
	return math.sqrt(Ret / ESS)

def SignificanceError(Value, Distibution, Node):

	print()
	print(f"Cannot calculate valid probably for value {Value} using {Distibution} on node, ", end='\t')
	TList = Node.GetTaxaList()
	TList.sort()

	print(*TList, sep=',')
	print()
	print("One cause of this can be, the x value is 0, if this happens it may be the all the taxa in the node have the same trait value.\n")
	print()
	sys.exit(1)


def GetParameterSignificance(Value, Distibution, NoEqZero, Thershold, N, CorrelR2, ESS, Node):
	
	if NoEqZero == 1.0:
		return 0, 0, 0, 0, 0

	try:
		P = math.exp(math.log(Distibution.PDF(math.fabs(Value))) + Thershold)
	except ValueError:
		SignificanceError(Value, Distibution, Node)


	ObsP = 1.0 - NoEqZero
	SDObsP = CaclSDP(ObsP, ESS)

	if SDObsP != 0:
		Z = (ObsP - P) / SDObsP
		Sig = 1 if Z >= Z_SIGNIFICANCE else 0
	else:
		Z = Z_OVERFLOW_VAL
		Sig = 1

	return P, ObsP, SDObsP, Z, Sig

def MakeSigInfo(P, ObsP, SDObsP, Z, Sig):
	SigInfo = {}
	SigInfo["P"] = P
	SigInfo["ObsP"] = ObsP
	SigInfo["SDObsP"] = SDObsP
	SigInfo["Z"] = Z
	SigInfo["Sig"] = Sig

	return SigInfo

def OutputBeta(Node, LogFileInfo, NData):

	BetaList = Node.BetaList

	Mean = numpy.average(BetaList)
	SD = numpy.std(BetaList)
	N = len(BetaList)

	NonZeroBeta = [x for x in BetaList if x != 0]

	MeanNZ = 0
	SDNZ = 0
	if len(NonZeroBeta) > 0:
		MeanNZ = numpy.average(NonZeroBeta)
		SDNZ = numpy.std(NonZeroBeta)
	
	NoGZero = len([x for x in BetaList if x > 0]) / N
	NoLZero = len([x for x in BetaList if x < 0]) / N
	NoEZero = len([x for x in BetaList if x == 0]) / N

	NData.extend([Mean, SD, MeanNZ, SDNZ, len(BetaList), LogFileInfo.ESS, NoLZero, NoEZero, NoGZero])

	P, ObsP, SDObsP, Z, Sig = GetParameterSignificance(MeanNZ, LogFileInfo.LandPrior, NoEZero, LogFileInfo.LandThreshold, N, LogFileInfo.LhR2, LogFileInfo.ESS, Node)

	Node.BetaSigInfo = MakeSigInfo(P, ObsP, SDObsP, Z, Sig)

	NData.extend([P, ObsP, SDObsP, BETA_Z_PLACE_HOLDER, BETA_SIG_PLACE_HOLDER])


def OutputSclars(Node, LogFileInfo, NData):
	
	ScalarList = Node.ScalarList
	Mean = numpy.average(ScalarList)
	SD = numpy.std(ScalarList)
	N = len(ScalarList)

	NonOneScalars = [x for x in ScalarList if x != 1] 
	MeanNO = 1
	SDNO = 0
	MedianNO = 1
	if len(NonOneScalars) > 0:
		MeanNO = numpy.average(NonOneScalars)
		SDNO = numpy.std(NonOneScalars)
		MedianNO = numpy.median(NonOneScalars)

	PctGOne = len([x for x in ScalarList if x > 1]) / N
	PctLOne = len([x for x in ScalarList if x < 1]) / N
	PctEOne = len([x for x in ScalarList if x == 1])/ N

	NData.extend([Mean, SD, MeanNO, MedianNO, SDNO, PctLOne, PctEOne, PctGOne])

	P, ObsP, SDObsP, Z, Sig = GetParameterSignificance(MeanNO, LogFileInfo.NodePrior, PctEOne, LogFileInfo.NodeThreshold, N, LogFileInfo.LhR2, LogFileInfo.ESS, Node)

	Node.ScalarSigInfo = MakeSigInfo(P, ObsP, SDObsP, Z, Sig)

	NData.extend([P, ObsP, SDObsP, SCALAR_Z_PLACE_HOLDER, SCALAR_SIG_PLACE_HOLDER])


def BetaNodeCoOccurrence(Node):
	
	ScalarList = Node.ScalarList
	BetaList = Node.BetaList

	Numerator = 0
	Both = []
	for NodeS, Beta in zip(ScalarList, BetaList):
		if NodeS != 1.0 or Beta != 0:
			Both.append((NodeS, Beta))

		if NodeS != 1.0 and Beta != 0:
			Numerator += 1

	Denominator = len(Both)

	if Denominator == 0:
		return 0, 0
	
	return Numerator / Denominator, Denominator

def OutputCoOccurrence(Node, NData):

	Freq, Denominator = BetaNodeCoOccurrence(Node)

	NData.extend([Freq, Denominator])

	return Freq, Denominator

def ThenBeta(Pos, Node):

	for N in Node:
		if N.BetaList[Pos] != 0.0:
			return True

	return False

def ThenNode(Pos, Node):

	for N in Node:
		if N.ScalarList[Pos] != 1.0:
			return True

	return False

def PropNodeThenBeta(Node):

	Ret = 0

	for Index, Val in enumerate(Node.ScalarList):
		if Val != 1.0:
			if ThenBeta(Index, Node) == True:
				Ret += 1

	return Ret / len(Node.ScalarList)

def PropBetaThenNode(Node):
	Ret = 0

	for Index, Val in enumerate(Node.BetaList):
		if Val != 0.0:
			if ThenNode(Index, Node) == True:
				Ret += 1

	return Ret / len(Node.BetaList)

def GetPropDecBeta(Node):
	Ret = 0
	for Index, Val in enumerate(Node.ScalarList):
		if ThenBeta(Index, Node) == True:
			Ret += 1

	return Ret / len(Node.ScalarList)

def RecGetNodeSumBL(Node):
	
	Ret = Node.BranchLength

	for N in Node:
		Ret += RecGetNodeSumBL(N)

	return Ret

def GetNodeSumBL(Node):

	Ret = 0

	for N in Node:
		Ret += RecGetNodeSumBL(Node)

	return Ret



def GetHeader():
	Ret = []
	Ret.extend(["ID", "Branch Length", "Height",  "Height Mid-Point", "Sum BL", "No Descendants"])
	Ret.extend(["No Concurrent Linages Start", "No Concurrent Linages Mid-point", "No Concurrent Linages End"])

	Ret.extend(["Mean (Beta * BL)","SD (Beta * BL)", "Mean (Beta * BL) NZ", "SD (Beta * BL) NZ", "Sample Size", "ESS", "P < 0", "P == 0", "P > 0"])
	Ret.extend(["P (Beta * BL)","Observed P (Beta * BL)","SD Observed P(Beta * BL)","Z (Beta * BL)","Sig (Beta * BL)"])

	Ret.extend(["Mean Scalar","SD Scalar","Mean Non 1 Scalar", "Median Non 1", "SD Non 1 Scalar", "P < 1", "P == 1", "P > 1"])
	Ret.extend(["P Scalar","Observed P Scalar","SD Observed P Scalar", "Z Scalar","Sig Scalar"])

	Ret.extend(["Prop beta followed by a variance scalar", "No Beta or Node"])
	Ret.extend(["Prop Variance scalar followed by a beta ", "Prop Beta on Dec Nodes"])
	Ret.extend(["Beta Scalar Trade off"])
	Ret.extend(["Md5 Sum","No Taxa"])

	return Ret

def SetTaxaHeader(Header, Data):

	MaxNoTaxa = max([len(x) for x in Data]) - len(Header)
	
	TaxaH = []
	for i in range(MaxNoTaxa):
		TaxaH.append(f"Taxa-{i+1}")

	Header.extend(TaxaH)

def PadData(Data):

	Ret = []

	MaxNoTaxa = max([len(x) for x in Data])

	for Row in Data:
		PadNo = MaxNoTaxa - len(Row)

		Row.extend("" * PadNo)

		Ret.append(Row)

	return Ret

def SaveData(Data, VarRatesFN):

	Header = GetHeader()
	if JOIN_TAXA_NAMES == False:
		SetTaxaHeader(Header, Data)
		Data = PadData(Data)
	else:
		Header.append("Taxa")

	df = pd.DataFrame(Data, columns = Header)
	df.set_index("ID", inplace=True)

	df.to_csv(f"{VarRatesFN}.csv")
	
def SetRowBetaZScalarZ(Row, BetaZ, SclarZ):
	Row[Row.index(BETA_Z_PLACE_HOLDER)] = BetaZ
	Row[Row.index(SCALAR_Z_PLACE_HOLDER)] = SclarZ

	Row[Row.index(BETA_SIG_PLACE_HOLDER)] = 1 if BetaZ >= Z_SIGNIFICANCE else 0
	Row[Row.index(SCALAR_SIG_PLACE_HOLDER)] = 1 if SclarZ >= Z_SIGNIFICANCE else 0


def AssingBetaNode(Row, Node, ZObsChancePairing, TradeOff):
	if ZObsChancePairing >= -2:
		SetRowBetaZScalarZ(Row, Node.BetaSigInfo["Z"], Node.ScalarSigInfo["Z"])

		if ZObsChancePairing >= 2:
			Row.append("Retain Both - positively correlated")
		else:
			Row.append("Retain Both - indeterminate")

		return

	if TradeOff >= 2:
		Row.append("Retain Beta")
		SetRowBetaZScalarZ(Row, Node.BetaSigInfo["Z"], 0.0)
		return
		
	if TradeOff <= -2:
		Row.append("Retain Node")
		SetRowBetaZScalarZ(Row, 0.0, Node.ScalarSigInfo["Z"])
		return

	SetRowBetaZScalarZ(Row, Node.BetaSigInfo["Z"], Node.ScalarSigInfo["Z"])
	Row.append("Retain Both - indeterminate")

def SetNodeBetaTraidOff(Node, Row, ESS, PropBetaVarScalar):
	
	if Node.BetaSigInfo["Sig"] == 0 or Node.ScalarSigInfo["Sig"] == 0:
		SetRowBetaZScalarZ(Row, Node.BetaSigInfo["Z"], Node.ScalarSigInfo["Z"])
		Row.append("No Change")
		return 

	BSInfo = Node.BetaSigInfo
	SSInfo = Node.ScalarSigInfo
		
	ChancePairing = BSInfo["ObsP"]  * SSInfo["ObsP"]

	if ChancePairing == 1.0:
		SetRowBetaZScalarZ(Row, Node.BetaSigInfo["Z"], Node.ScalarSigInfo["Z"])
		Row.append("Retain Both - both present 100% ")
		return 
		
	A = (ChancePairing + PropBetaVarScalar) / 2
	A = A * (1.0 - A)
	A = A / ESS
	ZObsChancePairing = (PropBetaVarScalar - ChancePairing) / math.sqrt(A)


	A = (BSInfo["ObsP"] + SSInfo["ObsP"]) / 2.0 
	A = A * (1.0 - A)
	A = math.sqrt(A / ESS)
	TradeOff = (BSInfo["ObsP"] - SSInfo["ObsP"]) / A

	AssingBetaNode(Row, Node, ZObsChancePairing, TradeOff)
	

def CreateData(Trees, LogFileInfo):

	Tree = Trees.TreeList[0]

#	PrintHeader()

	Data = []

	for Index, Node in enumerate(Tree):

		NData = []

		NonZero = [x for x in Node.BetaList if x != 0]
		Node.BetaPP = len(NonZero) / len(Node.BetaList)
		Node.BetaListBin = [0 if x == 0 else 1 for x in Node.BetaList]
				
		NData = [Index, Node.BranchLength, Node.Height, Node.Height + (Node.BranchLength * .5), GetNodeSumBL(Node), len(Node.NodeList)]

		CLin = Tree.LinagesAtHeight(Node.Height) 
		CLinMidPoint = Tree.LinagesAtHeight(Node.Height + (Node.BranchLength * .5)) 
		CLinEnd = Tree.LinagesAtHeight(Node.Height + Node.BranchLength) 

		NData.extend([CLin, CLinMidPoint, CLinEnd])

		OutputBeta(Node, LogFileInfo, NData)

		OutputSclars(Node, LogFileInfo, NData)

		PropBetaVarScalar, _ = OutputCoOccurrence(Node, NData)
		
		NodeTBeta = PropNodeThenBeta(Node)
		PropOnDecBeta = GetPropDecBeta(Node)

		NData.extend([NodeTBeta, PropOnDecBeta])
		
		SetNodeBetaTraidOff(Node, NData, LogFileInfo.ESS, PropBetaVarScalar)

		TList = Node.GetTaxaList()
		TList.sort()
		
		NData.extend([GetMD5Sum(TList), len(TList)])

		if JOIN_TAXA_NAMES == False:
			NData.extend(TList)
		else:
			NData.append('|'.join(TList))

		
		Data.append(NData)

	return Data

def SetTipHeight(Tree):

	for Node in Tree:
		if Node.Tip == True:
			Node.Height = 0

def main():
	
	if len(sys.argv) != 4:
		print(f"{sys.argv[0]} takes a tree file, a VarRates file and a log file.")
		sys.exit(0)

	(Trees, Data) = PassVarRates(sys.argv[1], sys.argv[2])
	LogFileInfo = PassLogFile(sys.argv[3])
	
	#SetTipHeight(Trees[0])


	Data = CreateData(Trees, LogFileInfo)

	SaveData(Data, sys.argv[2])



# pyinstaller --onefile --hidden-import pandas._libs FabricPostProcessor.py

# OS X 
# pyinstaller --onefile --hidden-import pandas._libs --target-architecture x86_64 FabricPostProcessor.py

if __name__ == "__main__":
	main()