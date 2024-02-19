import sys
import pandas as pd


def CreateHeader(FileList):
	BetaName = [f"Z (Beta * BL) - {FileName}" for FileName in FileList]
	NodeName = [f"Z Scalar - {FileName}" for FileName in FileList]

	return ["Md5 Sum", "No Sig Beta", *BetaName, "No Sig Nodes", *NodeName]


def FastMerge(FileList):
	df = pd.read_csv(FileList[0], low_memory=False)

	checksum = df["Md5 Sum"].copy()
	sig_beta = df["Sig (Beta * BL)"].copy()
	sig_node = df["Sig Scalar"].copy()

	beta_z = [df["Z (Beta * BL)"].copy()]
	node_z = [df["Z Scalar"].copy()]

	for fname in FileList[1:]:
		df_run = pd.read_csv(fname, low_memory=False)

		if df["Md5 Sum"].equals(df_run["Md5 Sum"]) == False:
			print(f"Md5 Sum differs in {fname}")
			sys.exit(1)
		
		sig_beta = sig_beta.add(df_run["Sig (Beta * BL)"])
		sig_node = sig_node.add(df_run["Sig Scalar"])

		beta_z.append(df_run["Z (Beta * BL)"].copy())
		node_z.append(df_run["Z Scalar"].copy())

	header = CreateHeader(FileList)
	output = pd.DataFrame([checksum, sig_beta, *beta_z, sig_node, *node_z])
	output = output.transpose()

	output = output.set_axis(header, axis=1)
	output.set_index("Md5 Sum", inplace=True)

	output.to_csv("out.csv")


#pyinstaller --onefile --hidden-import pandas._libs MergeFabric.py
# OS X 
# pyinstaller --onefile --hidden-import pandas._libs --target-architecture x86_64 MergeFabric.py

FileList = sys.argv[1:]

FastMerge(FileList)
