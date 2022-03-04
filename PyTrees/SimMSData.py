import sys
import math
import random

def CreatePMat(t, NOS):
	Temp = (t * NOS) / (NOS - 1)
	Temp = math.exp(-Temp)
	InvNOS = 1/NOS;


	Miss = InvNOS - (InvNOS * Temp);
	Hit = (((NOS-1) / NOS) * Temp) + InvNOS;

	P = [[Miss for i in range(NOS)] for j in range(NOS)]
	for i in range(NOS):
		P[i][i] = Hit

	return P

def SetNodePMatrix(Node, StateSet):

	Node.PDic = {}

	for No in StateSet:
		Node.PDic[No] = CreatePMat(Node.BranchLength, No)

def SetPMatrix(Tree, NoStateList):

	UStateList = list(set(NoStateList))

	for Node in Tree:
		if Node != Tree.Root:
			SetNodePMatrix(Node, UStateList)

def GetNewState(PList, Point):

	CPoint = 0
	Index = 0
	for Val in PList:
		if Point >= CPoint and Point <= CPoint + Val:
			return Index
		CPoint += Val
		Index += 1

	print("Prop errror.\n")
	sys.exit(1)
	return -1
			

def SimNode(Node, NOS):
	PMat =  Node.PDic[NOS]
	AnsState = Node.Anc.State
	Point = random.random()

	Node.State = GetNewState(PMat[AnsState], Point)

	for N in Node:
		SimNode(N, NOS)

	if Node.Tip == True:
		Node.Data.append(Node.State)

def RunData(Tree, NoStateList):

	Root = Tree.Root
	for NOS in NoStateList:
		Root.State = random.randint(0, NOS-1)
		for Node in Root:
			SimNode(Node, NOS)

def BlankTipData(Tree):

	for Node in Tree:
		if Node.Tip == True:
			Node.Data = []

def MakeTaxaBin(Node, NOSList):
	
	BinList = []
	Index = 0
	for State in Node.Data:
		BinData = [0] * NOSList[Index]
		BinData[State] = 1
		BinList.extend(BinData)
		Index += 1

	Node.Data = BinList

def MakeDataBin(Tree, NOSList):
	for Node in Tree:
		if Node.Tip == True:
			MakeTaxaBin(Node, NOSList)


def SimMSData(Tree, NoStateList, MakeBinary=False):
	BlankTipData(Tree)
	SetPMatrix(Tree, NoStateList)
	RunData(Tree, NoStateList)

	if MakeBinary == True:
		MakeDataBin(Tree, NoStateList)



