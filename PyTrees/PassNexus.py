import sys, re

from .StringToTree import StringToTree

class PassNexus:

	def __init__(self, FName, NType):
		self.RootList = []
		self.TaxaLookup = {}

		self.PassFile(FName, NType)
			
	def AddTaxa(self, TaxaMatch):

		Num = int(TaxaMatch.group(1))
		Taxa = TaxaMatch.group(2)
		
		Taxa = Taxa.replace(",", "")
		Taxa = Taxa.replace(";", "")
				
		self.TaxaLookup[Num] = Taxa

	def MapTaxa(self, Node):

		if Node.Tip == True:
			Node.TipID = int(Node.Taxa)
			Node.Taxa = self.TaxaLookup[Node.Taxa]


		for N in Node.NodeList:
			self.MapTaxa(N)

	def AddTree(self, TressStr, NType):

		Root = StringToTree(TressStr, NType)
		self.MapTaxa(Root)
		self.RootList.append(Root)

	def GetNexusTreeStr(self, Str):
		
		TreeRe = re.compile("(\(.+;)")

		TreeMatch = re.search(TreeRe , Str)

		if not TreeMatch:
			raise Exception("cannot find valid tree in %s" % Str) 

		return TreeMatch[0]
		
	def PassFile(self, FName, NType):
	
			Start = False
		
			BeginTreesRE = re.compile("\s*begin\s+trees", re.IGNORECASE)
			EndTreesRE = re.compile("\s*(end|endblock)", re.IGNORECASE)
			TaxaRE = re.compile("\s*(-?[0-9]+)\s+(.+)")
			TreeRE = re.compile("^\s*((u|r)*(tree))\s+", re.IGNORECASE)

			FIn = open(FName, "r")
		
			for Line in FIn:

				Line = Line.rstrip()
#				Line = re.sub("\[.+?\]", '', Line)

				BeginMatch = BeginTreesRE.match(Line)
				if BeginMatch:
					Start = True

				EndMatch = EndTreesRE.match(Line)
				if EndMatch:
					Start = False
					
				if Start == True:
					TaxaMatch = TaxaRE.match(Line)
					if TaxaMatch:
						self.AddTaxa(TaxaMatch)
					
					TreeMatch = TreeRE.match(Line)

					if TreeMatch:
						TreeStr = self.GetNexusTreeStr(Line)
						self.AddTree(TreeStr, NType)

			FIn.close()


