import sys
import pandas as pd

def sum_sig(sig_cols):

	sum_col = sig_cols[0]

	for col in sig_cols[1:]:
		sum_col = sum_col + col

	return sum_col

def append_cols(master, cols, file_names):

	cols = reversed(cols)
	file_names = reversed(file_names)
	for col, file_name in zip(cols, file_names):
		col_name = f"{col.name} - {file_name}"
		master.insert(loc=0, column=col_name, value=col)

files = sys.argv[1:]

sig_beta = []
sig_node = []
beta_z = []
node_z = []
beta_non0 = []
node_non1 = []

for file_name in files:
	df = pd.read_csv(file_name, low_memory=False)
	
	sig_beta.append(df["Sig (Beta * BL)"])
	sig_node.append(df["Sig Scalar"])
	beta_z.append(df["Z (Beta * BL)"])
	node_z.append(df["Z Scalar"])
	beta_non0.append(df["Mean (Beta * BL) NZ"])
	node_non1.append(df["Mean Non 1 Scalar"])
	 

master = pd.read_csv(files[0], low_memory=False)
master.set_index("ID", inplace=True)

master.insert(loc=0, column="", value="")

append_cols(master, node_z, files)
master.insert(loc=0, column=" ", value="")

append_cols(master, beta_z, files)
master.insert(loc=0, column="  ", value="")

append_cols(master, node_non1, files)
master.insert(loc=0, column="   ", value="")

append_cols(master, beta_non0, files)
master.insert(loc=0, column="    ", value="")

sum_sig_beta = sum_sig(sig_beta)
sum_sig_nodes = sum_sig(sig_node)

master.insert(loc=0, column="No Sig Nodes", value=sum_sig_nodes)
master.insert(loc=0, column="No Sig Betas", value=sum_sig_beta)

master.to_csv("master.csv")
