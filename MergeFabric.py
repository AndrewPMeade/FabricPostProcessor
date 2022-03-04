import sys
import pandas as pd

def LoadTables(FileList):
	Tables = [pd.read_csv(FileName) for FileName in FileList]
		
	MasterSet = set(Tables[0]["Md5 Sum"]) 

	Index = 1
	for Table in Tables[1:]:

		Diff = MasterSet ^ set(Table["Md5 Sum"]) 

		if Diff:
			print(f"Differences btween checksums for files {FileList[0]} and {FileList[Index]}")
			print(*Diff)
			sys.exit(1)

		Index += 1

	[df.set_index(df["Md5 Sum"], inplace=True) for df in Tables]
	return Tables

def CreatRow(Md5Sum, Tables):

	

	BetaSig = 0
	BetaZ = []

	NodeSig = 0
	NodeZ = []

	for df in Tables:
		row = df.loc[Md5Sum]
	
		if row["Sig (Beta * BL)"] == 1:
			BetaSig += 1

		BetaZ.append(float(row["Z (Beta * BL)"]))

		if row["Sig Scalar"] == 1:
			NodeSig += 1

		NodeZ.append(float(row["Z Scalar"]))
				
	return [Md5Sum, BetaSig, *BetaZ, NodeSig, *NodeZ]

def CreateHeader(FileList):
	BetaName = [f"Z (Beta * BL) - {FileName}" for FileName in FileList]
	NodeName = [f"Z Scalar - {FileName}" for FileName in FileList]

	return ["Md5 Sum", "No Sig Beta", *BetaName, "No Sig Nodes", *NodeName]

def BuildCombinedDF(Tables, FileList):

	
	Data = []
	for Md5Sum in Tables[0]["Md5 Sum"]:
		Data.append(CreatRow(Md5Sum, Tables))

	return Data

def Output(Data, FileList):

	Header = CreateHeader(FileList)	

	DF = pd.DataFrame(Data, columns=Header)
	DF.set_index("Md5 Sum", inplace=True)
	DF.to_csv(sys.stdout,line_terminator='\n')


#pyinstaller --onefile --hidden-import pandas._libs MergeFabric.py

sys.stdout = open("out.csv", "w")

FileList = sys.argv[1:]

Tables = LoadTables(FileList)

Data = BuildCombinedDF(Tables, FileList)
#for row in Data:
#	print(row)
Output(Data, FileList)