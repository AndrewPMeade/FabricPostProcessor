import sys
import ply.lex as lex
from .Node import Node

#tokens = ("lParen", "rParen", "Comma", "Length", "Colon", "SimiColon", "Name", "Rate", "Comment")
tokens = ("Comment", "lParen", "rParen", "Comma", "Length", "Colon", "SimiColon",  "Name", "Rate")
#tokens = ("lParen", "rParen", "Comma", "Length", "Colon", "SimiColon", "Name")

t_Comment = "(\[.+?\])"
t_lParen = r'\('
t_rParen = r'\)'
t_Comma = ","
t_Colon = ":"
t_SimiColon = ";"

t_Name = "(\w|\_|\ |\'|-|[0-9]|[.])+"

def t_error(t):
    print("Illegal character '%s'" % t.value[0])
    t.lexer.skip(1)

def t_Length(t):
	"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
	t.value = float(t.value)
	return t

def t_Rate(t):
	"@(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
	Str = str(t.value)
	t.value = float(Str.replace("@", ""))
	return t

#def t_Comment(t):	
#	"\[.+?\]"
#	ComStr = str(t.value)
#	print(ComStr)
#	sys.exit(0)

def GetComment(Node, TokeList):

	if len(TokeList) == 0:
		return

	if TokeList[0].type != "Comment":
		return

	Comment = TokeList.pop(0)
	Node.Comment = Comment.value


def GetBL(N, TokeList):

	GetComment(N, TokeList)

	if TokeList[0].type == "Colon":
		TokeList.pop(0)
		Len = TokeList.pop(0)
		N.BranchLength = Len.value

	GetComment(N, TokeList)

	if TokeList[0].type == "Rate":
		Rate = TokeList.pop(0)
#		Node.Rate = Rate.value
		N.Rate = Rate.value

	GetComment(N, TokeList)


def GeTip(N, TokeList):

	GetComment(N, TokeList)
	Toke = TokeList.pop(0)
	GetComment(N, TokeList)
	N.Taxa = Toke.value
	N.Tip = True
	GetBL(N, TokeList)



def RecBuildTree(Node, TokeList, NType):
		
	GetComment(Node, TokeList)

	while True:
		TokeList.pop(0)

		NNode = NType()
		Node.AddNode(NNode)
		
		GetComment(NNode, TokeList)

		Toke = TokeList[0]
		if Toke.type == "lParen":
			RecBuildTree(NNode, TokeList, NType)
		else:
			GeTip(NNode, TokeList)

		if TokeList[0].type != "Comma":
			break
	TokeList.pop(0)
		
	GetBL(Node, TokeList)
	
	GetComment(Node, TokeList)


def BuildTree(TokeList, NType):
	
	Root = NType()

	RecBuildTree(Root, TokeList, NType)

	return Root

def PrintTaxa(Node):
	if Node.Tip == True:
		print(Node.Taxa, Node.BranchLength)
		return

	for N in Node.NodeList:
		print("Node:\t", N.BranchLength)
		PrintTaxa(N)
		
def PassTree(TreeStr, NType):
	Lexer = lex.lex()
	
	Lexer.input(TreeStr)
	
	TokeList = list(Lexer)
		
	Root = BuildTree(TokeList, NType)
	
	return Root

def StringToTree(TreeStr, NType=Node):
	return PassTree(TreeStr, NType)



