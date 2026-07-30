import scipy.io as sio
import numpy as np
import random as rand
from Bio import Entrez


Entrez.email = "xmzno1@stu.xjtu.edu.cn"

def test():
	network_array = np.zeros((670,40230),dtype=float)

	for i in range(1,2):
		filename = "yrvm_exp_"+str(i)+".mat"
		print('Integrating '+filename+'...')
		mat_content = sio.loadmat(filename)
		network_array[670*(i-1):len(mat_content['yrvm_exp'])*i,:]=mat_content['yrvm_exp']
	
	sio.savemat('gene_net.mat',{'yrvm_exp':network_array})
	print('done')

def matToNparray():
    network_array = np.zeros((40230,40230),dtype=float)
    #geneSeriesNo = 0
    for i in range(1, 62):
        filename = "yrvm_exp_"+str(i)+".mat"
        mat_content = sio.loadmat(filename)
        mat_array = mat_content['yrvm_exp']
        for j in range(len(mat_array)):
            network_array[670*(i-1) + j] = mat_array[j]
            print("\rGene No. ", 670*(i-1) + j, "added.", end='')
    print("\nwritting to file...")
    np.save('rvm_net.npy', network_array, allow_pickle=False, fix_imports=False)

def matToText():
	network_array = np.zeros((670,40230),dtype=float)

	for i in range(1,2):
		filename = "yrvm_exp_"+str(i)+".mat"
		print('Opening '+filename+'...')
		mat_content = sio.loadmat(filename)
		network_array=mat_content['yrvm_exp']
		for rn,row in enumerate(network_array):
			geneidx = 670*(i-1)+rn
			with open("txtmat/"+str(geneidx),'w') as w:
				for cn,col in enumerate(row):
					w.write(str(rn)+'\t'+str(cn)+'\t'+str(col)+'\n')
	print('done')
	
def getRandGenes(onconum=50,randnum=450,filesufix='-1'):
	oncogene = rand.sample(range(0,699),onconum)
	randgene = rand.sample(range(0,40230),randnum)
	oncogeneinfo = []
	with open('oncogenes.txt') as r:
		for line in r.readlines():
			oncogeneinfo.append(line.strip().split())
	mat_content = sio.loadmat('gene_index.mat')
	gene_array = np.zeros((40230),dtype=int)
	gene_array = mat_content['gene_index']
	with open('txtnet/testgene'+str(onconum)+'-'+str(randnum)+filesufix+'.txt','w') as w:
		for gi in oncogene:
			w.write(oncogeneinfo[gi][1]+'\t'+oncogeneinfo[gi][2]+'\n')
		for gi in randgene:
			w.write(str(gene_array[gi])+'\t'+str(gi)+'\n')
	print('done.')#gene duplicates may happen

def addRandGenes(orifile='aml/933124/symbols',totalnum=600,fileprefix='aml/933124/testgene-',filesufix='-1'):
	oncogeneinfo = []
	with open(orifile+'.txt') as r:
		for line in r.readlines():
			oncogeneinfo.append(line.strip().split())
	onconum=len(oncogeneinfo)
	randnum=totalnum-onconum
	randgene = rand.sample(range(0,40230),randnum)
	mat_content = sio.loadmat('gene_index.mat')
	gene_array = np.zeros((40230),dtype=int)
	gene_array = mat_content['gene_index']
	with open(fileprefix+str(onconum)+'-'+str(randnum)+filesufix+'.txt','w') as w:
		for gi in oncogeneinfo:
			w.write(gi[0]+'\t'+gi[1]+'\n')
		for gi in randgene:
			w.write(str(gene_array[gi])+'\t'+str(gi)+'\n')
	print('done.')#gene duplicates may happen

def getAmlRandGenes(filedir,totalnum=1000):
	basedir='geneinfo/aml_by_rand_pp/0/'
	fdir=basedir+filedir
	oncofile=fdir+'/symbols'
	oncogeneinfo = []
	with open(oncofile) as r:
		for line in r.readlines():
			oncogeneinfo.append(line.strip().split())
	randnum=totalnum-len(oncogeneinfo)
	randgene = rand.sample(range(0,40230),randnum)
	mat_content = sio.loadmat('gene_index.mat')
	gene_array = np.zeros((40230),dtype=int)
	gene_array = mat_content['gene_index']
	with open(fdir+'/testgene-'+str(totalnum-randnum)+'-'+str(totalnum)+'.txt','w') as w:
		for gi in oncogeneinfo:
			w.write(gi[0]+'\t'+gi[1]+'\n')
		for gi in randgene:
			w.write(str(gene_array[gi])+'\t'+str(gi)+'\n')
	print('done.')#gene duplicates may happen

def getGeneSymbInfo(genesymbol,searchBy='name'):
	genename=''
	entrezid=''
	rvmidx=''
	ensembleid=''
	idx=[]
	print('searching gene: '+genesymbol)
	with open('matnet/rvm_entrez_idx.txt') as r:
		for line in r.readlines():
			idx.append(line.strip().split()[0])
	idmap=[]
	with open('matnet/identifier_mappings.txt') as r:
		for line in r.readlines():
			idmap.append(line.strip().split('\t'))
	namefound=False
	if(searchBy=='name'):
		entrezidfound=False
		for el in idmap:
			if((not namefound) and el[2]=='Entrez Gene ID'):
				ensembleid=el[0]
				entrezid=el[1]
			if(el[1]==genesymbol and el[2]=='Gene Name'):
				genename=el[1]
				namefound=True
				if (ensembleid==el[0]):
					entrezidfound=True
					break
				else:
					ensembleid=el[0]
					entrezid=''
					continue
			if(el[0]==ensembleid and el[2]=='Entrez Gene ID' and namefound):
				entrezid=el[1]
				entrezidfound=True
		if entrezidfound:
			for eidx,eid in enumerate(idx):
				if entrezid==eid:
					rvmidx=eidx
		if not namefound:
			genename=genesymbol
			ensembleid=''
			entrezid=''
		if entrezid=='':
			while True:
				try:
					handle = Entrez.esearch(db='gene',term=genename+'[gene] AND human[orgn]')
					record = Entrez.read(handle)
				except:
					print('Error happening!')
					continue
				break
			if(int(record['Count'])>=1):
				entrezid=record['IdList'][0]
				for eidx,eid in enumerate(idx):
					if entrezid==eid:
						rvmidx=eidx
		return [genename,entrezid,rvmidx,ensembleid]
	elif(searchBy=='rvmidx'):
		entrezidfound=False
		for eidx,eid in enumerate(idx):
				if genesymbol==eidx:
					entrezid=eid
					rvmidx=eidx
					entrezidfound=True
		if entrezidfound:
			for el in idmap:
				if(el[1]==entrezid and el[2]=='Entrez Gene ID'):
					ensembleid=el[0]
				if(el[0]==ensembleid and el[2]=='Gene Name'):
					genename=el[1]
		return [genename,entrezid,rvmidx,ensembleid]

def findGeneByPos(chr,start,end='',orien=''):
	geneannot=[]
	with open('matnet/gene_result.txt') as r:
		for line in r.readlines():
			geneannot.append(line.split('\t'))
	#Entrez id:col[2];	Gene symbol:col[5];	Chromosome:col[10];	Start:col[12];	End:col[13];	Orientation[14]
	results=[]
	for gene in geneannot:
		foundstart=False
		foundend=False
		overlap=False
		if(gene[12]!='' and gene[13]!=''):
			if (gene[10]==chr and int(start) >= int(gene[12]) and int(start) <=int(gene[13])):
				if(orien!=''):
					if((orien=='+' or orien=='plus') and gene[14]=='plus'):
						foundstart=True
					elif((orien=='-' or orien=='minus') and gene[14]=='minus'):
						foundstart=True
				else:
					foundstart=True
			if(gene[10]==chr and end!='' and int(end) >= int(gene[12]) and int(end) <=int(gene[13])):
				if(orien!=''):
					if((orien=='+' or orien=='plus') and gene[14]=='plus'):
						foundend=True
					elif((orien=='-' or orien=='minus') and gene[14]=='minus'):
						foundend=True
				else:
					foundend=True
		if(foundstart or foundend):
			results.append([gene[2],gene[5],gene[10],gene[12],gene[13],gene[14]])
			if(not(foundstart and foundend) and end!=''):
				overlap=True
			results[-1].append(overlap)
	return results

def chkIsOnco(file='aml/933124/symbols.txt', oncofile='oncogenes.txt'):
	genelist=[]
	oncolist=[]
	with open(file) as r:
		for line in r.readlines():
			genelist.append(line.strip().split())
	with open(oncofile) as r:
		for line in r.readlines():
			oncolist.append(line.strip().split())
	for gene in genelist:
		gene.append('0')
		for onco in oncolist:
			if gene[2]==onco[1]:
				gene[-1]='1'
				break
	with open(file, 'w') as w:
		w.write('\n'.join('\t'.join(el for el in gene) for gene in genelist))

'''Need to connect to NCBI database
def searchEnsgFromEntrezDbByGeneSymbol(symb):
    entrezid = ''
    while True:
        try:
            handle = Entrez.esearch(db='gene',term=symb+'[gene] AND human[orgn]')
            record = Entrez.read(handle)
        except:
            print('Error happening!')
            continue
        break
    if(int(record['Count'])>=1):
        entrezid = record['IdList'][0]
    while True:
        try:
            stream = Entrez.efetch(db="gene", id=entrezid, retmode='xml')
            res = Entrez.read(stream)
        except:
            print('Error happening!')
            continue
        break
    for row in res[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_db']:
        if row['Dbtag_db'] == 'Ensembl':
            return row['Dbtag_tag']['Object-id']['Object-id_str']
    return False

def searchEnsgFromEntrezDbByEntrezId(eid):
    while True:
        try:
            stream = Entrez.efetch(db="gene", id=eid, retmode='xml')
            res = Entrez.read(stream)
        except:
            print('Error happening!')
            continue
        break
    for row in res[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_db']:
        if row['Dbtag_db'] == 'Ensembl':
            return row['Dbtag_tag']['Object-id']['Object-id_str']
    return False
'''

def getMappingList():
    mapping_list = []
    with open("identifier_mappings.txt") as file:
        for line in file:        
            l=line.strip().split('\t')
            mapping_list.append(l)
    #names_header = mapping_list.pop(0)
    return mapping_list

def getGeneEntrezidBySymb(symbol, mapping_list):  
    for idx, item in enumerate(mapping_list):
        if item[2].startswith('Gene Name') and item[1] == symbol:
            return mapping_list[idx - 1][1]
    return 'False'

def getGeneEnsgBySymb(symb, mapping_list):
    ensglist = []
    for idx, item in enumerate(mapping_list):
        if item[2].startswith('Ensembl Gene ID') and item[1].startswith('ENSG'):
            ensglist.append(mapping_list[idx])
        if item[2].startswith('Gene Name') and item[1] == symb:
            return ensglist[-1][1]
    return 'False'

def getGeneSymbByEntrezid(entzid, mapping_list):
    for idx, item in enumerate(mapping_list):
        if item[2].startswith('Entrez Gene ID') and item[1] == entzid:
            return mapping_list[idx + 1][1]
    return 'False'

def getGeneEnsgByEntrezid(symb, mapping_list):
    ensglist = []
    for idx, item in enumerate(mapping_list):
        if item[2].startswith('Ensembl Gene ID') and item[1].startswith('ENSG'):
            ensglist.append(mapping_list[idx])
        if item[2].startswith('Entrez Gene ID') and item[1] == symb:
            return ensglist[-1][1]
    return 'False'

def getComnetEnsgMappingList():
    idx=[]
    with open('com_net_ensg_nodes.txt') as r:
        for line in r.readlines():
            idx.append(line.strip().split()[0])
    return idx

def getGeneComnetidxByEnsg(ensg, comnetMappingList):
    for eidx,eid in enumerate(comnetMappingList):
        if ensg==eid:
            comnetidx=eidx
            return comnetidx
    return 'False'

def getRvmIdxSymbolMappingList():
    idx=[]
    with open('matnet/rvm_entrez_symb_idx.tsv') as r:
        for line in r.readlines():
            idx.append(line.strip().split('\t'))
    return idx

def getGeneRvmidxBySymbol(symbol, rvmIdxSymbolMappingList):
    for eidx,eid in enumerate(rvmIdxSymbolMappingList):
        if symbol==eid[1]:
            rvmidx=eidx
            return rvmidx
    return 'False'

def getRvmIdxMappingList():
    idx=[]
    with open('matnet/rvm_entrez_ensg_idx.tsv') as r:
        for line in r.readlines():
            idx.append(line.strip().split('\t'))
    return idx

def getGeneRvmidxByEntrezid(entrezid, rvmIdxMappingList):
    for eidx,eid in enumerate(rvmIdxMappingList):
        if entrezid==eid[0]:
            rvmidx=eidx
            return rvmidx
    return 'False'

def getGeneRvmidxByEnsg(ensg, rvmIdxMappingList):
    for eidx,eid in enumerate(rvmIdxMappingList):
        if ensg==eid[1]:
            rvmidx=eidx
            return rvmidx
    return 'False'