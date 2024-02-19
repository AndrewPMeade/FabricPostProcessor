import hashlib

def RecGetTaxa(Node, TList):

	if Node.Tip == True:
		TList.append(Node.Taxa)

	for N in Node:
		RecGetTaxa(N, TList)
	
def GetMD5Sum(TaxaList):

	TaxaList.sort()

	md5 = hashlib.md5()

	for Taxa in TaxaList:
		md5.update(Taxa.encode('utf-8'))

	return md5.hexdigest()
