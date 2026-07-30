import json
from lib.Mat_op import *
import csv
import random
import lib.Processing as pr

def get_spec_cancer_cgc_genes(cancer_name='none'):
    cancer_cgc_genes_info = pr.load_txt_to_list('./cgc/cgc_genes_raw-763.csv',
                                                filetype='csv')
    cancer_cgcs_symb = []
    for row in cancer_cgc_genes_info:
        somatic_types = row[9]
        germline_types = row[10]
        if cancer_name in somatic_types or cancer_name in germline_types:
            cancer_cgcs_symb.append(row[0])
    return cancer_cgcs_symb
    
# with open('data/sim_dict_density_alledges_ww.txt') as json_file:
#     sim_dict_frob_nl = json.load(json_file)

# sim_value_min_key1 = 0.0
# sim_value_min_key2 = 0.0
# sim_value_min = float('inf')

# for i in list(sim_dict_frob_nl.keys()):
#     temp_value_min = min(list(sim_dict_frob_nl[i].values()))
#     temp_key2 = list(sim_dict_frob_nl[i].keys())[list(sim_dict_frob_nl[i].values()).index(temp_value_min)]
#     if temp_value_min < sim_value_min:
#         sim_value_min = temp_value_min
#         sim_value_min_key1 = i
#         sim_value_min_key2 = temp_key2
# print(sim_value_min_key1, sim_value_min_key2, sim_value_min)

#process cgc gene info file, add ensg name
# filelines = []
# with open('cgc_genes_symbs.csv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile)
#     for row in spamreader:
#         filelines.append(row)
# for key, row in enumerate(filelines):
#     synonyms = row[2].split(',')
#     temp_str = ''
#     for field in synonyms:
#         if field.startswith('ENSG'):
#             temp_str = field[:field.index('.')]
#     if len(temp_str) > 0:
#         filelines[key].append(temp_str)
#     else:
#         filelines[key].append("Ensemble gene name")
# with open('cgc_genes_nameinfo.csv', 'w', newline='') as csvfile:
#     spamwriter = csv.writer(csvfile)
#     spamwriter.writerows(filelines)

#check cgc genes in combine net
# cgcgenes = []
# with open('cgc_genes_nameinfo.csv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile)
#     for row in spamreader:
#         cgcgenes.append(row)
# comnet_genes = []
# with open('com_net_ensg_nodes.txt') as txtfile:
#     for line in txtfile.readlines():
#         comnet_genes.append(line.strip())
# cgc_comnet = []
# for cgcgene in cgcgenes:
#     if cgcgene[2] in comnet_genes:
#         cgc_comnet.append(cgcgene[2])
# with open('cgc_genes_in_comnet.txt', 'w') as txtfile:
#     txtfile.write('\n'.join(cgc_comnet))

# #check cgc genes in rvm net
# cgcgenes = []
# with open('cgc_genes_nameinfo.csv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile)
#     for row in spamreader:
#         cgcgenes.append(row)
# rvm_genes = []
# with open('matnet/rvm_entrez_idx.txt') as txtfile:
#     for line in txtfile.readlines():
#         rvm_genes.append(line.strip())
# cgc_rvm = []
# for cgcgene in cgcgenes:
#     if len(cgcgene[1]) > 0 and cgcgene[1] in rvm_genes:
#         cgc_rvm.append(cgcgene[1])
# with open('cgc_genes_in_rvmnet_entrezid.txt', 'w') as txtfile:
#     txtfile.write('\n'.join(cgc_rvm))

#translate cgc rvm gene entrez id to ensemble id
# cgc_rvm_entrez = []
# with open('cgc_genes_in_rvmnet_entrezid.txt') as txtfile:
#     for line in txtfile.readlines():
#         cgc_rvm_entrez.append(line.strip())
# mapping_list = getMappingList()
# cgc_rvm_ensg = []
# for key, gene in enumerate(cgc_rvm_entrez):
#     cgc_rvm_ensg.append(getGeneEnsgByEntrezid(gene, mapping_list))
#     print("\rtranslated", key+1, "genes", end='')
# with open('cgc_genes_in_rvmnet_ensg.txt', 'w') as txtfile:
#     for key, gene in enumerate(cgc_rvm_ensg):
#         if gene != False:
#             txtfile.write(cgc_rvm_entrez[key] + '\t' + gene + '\n')

#translate cgc rvm gene entrez id to gene symble
def get_cgc_in_rvm_symb():
    cgc_rvm_entrez = []
    with open('cgc_genes_in_rvmnet_entrezid.txt') as txtfile:
        for line in txtfile.readlines():
            cgc_rvm_entrez.append(line.strip())
    mapping_list = getMappingList()
    cgc_rvm_symb = []
    for key, gene in enumerate(cgc_rvm_entrez):
        cgc_rvm_symb.append(getGeneSymbByEntrezid(gene, mapping_list))
        print("\rtranslated", key+1, "genes", end='')
    with open('cgc_genes_in_rvmnet_symb.txt', 'w') as txtfile:
        for key, gene in enumerate(cgc_rvm_symb):
            if gene != False:
                txtfile.write(cgc_rvm_entrez[key] + '\t' + gene + '\n')

#translate rvm gene entrez id to gene symble
def rvm_entrez_to_symb():
    rvm_entrez = []
    with open('matnet/rvm_entrez_idx.txt') as txtfile:
        for line in txtfile.readlines():
            rvm_entrez.append(line.strip())
    mapping_list = getMappingList()
    rvm_symb = []
    for key, gene in enumerate(rvm_entrez):
        rvm_symb.append([gene, getGeneSymbByEntrezid(gene, mapping_list)])
        print("\rtranslated", key+1, "genes", end='')
    with open('matnet/rvm_entrez_symb_idx.tsv', 'w') as txtfile:
        txtfile.write('\n'.join('\t'.join(str(item) for item in row) for row in rvm_symb))

#check and save cgc rvm/comnet genes in smoke data
# smoke_genes = []
# cgc_comnet_genes = []
# cgc_rvm_genes = []
# with open('GENIE3/somke_gene_test/smoke_data_ensg.csv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile)
#     for row in spamreader:
#         smoke_genes.append(row[0])
# with open('cgc_genes_in_rvmnet_ensg.txt', newline='') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter='\t')
#     for row in spamreader:
#         cgc_rvm_genes.append(row[1])
# with open('cgc_genes_in_comnet.txt') as txtfile:
#     for line in txtfile.readlines():
#         cgc_comnet_genes.append(line.strip())
# cgc_rvm_smoke_genes = []
# cgc_comnet_smoke_genes = []
# for gene in cgc_rvm_genes:
#     if gene in smoke_genes:
#         cgc_rvm_smoke_genes.append(gene)
# for gene in cgc_comnet_genes:
#     if gene in smoke_genes:
#         cgc_comnet_smoke_genes.append(gene)
# with open('GENIE3/somke_gene_test/cgc_rvm_smoke_genes.txt', 'w') as txtfile:
#     txtfile.write('\n'.join(cgc_rvm_smoke_genes))
# with open('GENIE3/somke_gene_test/cgc_comnet_smoke_genes.txt', 'w') as txtfile:
#     txtfile.write('\n'.join(cgc_comnet_smoke_genes))

#check and save cgc rvm/comnet genes in lung data
def save_cgc_rvm_in_lung_data():
    lung_genes = []
    #cgc_comnet_genes = []
    cgc_rvm_genes = []
    with open('GENIE3/lung_cancer/lung_data_symb.tsv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            lung_genes.append(row[0])
    with open('cgc_genes_in_rvmnet_symb.txt', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            cgc_rvm_genes.append(row[1])
    # with open('cgc_genes_in_comnet.txt') as txtfile:
    #     for line in txtfile.readlines():
    #         cgc_comnet_genes.append(line.strip())
    cgc_rvm_lung_genes = []
    #cgc_comnet_smoke_genes = []
    for gene in cgc_rvm_genes:
        if gene in lung_genes:
            cgc_rvm_lung_genes.append(gene)
    # for gene in cgc_comnet_genes:
    #     if gene in smoke_genes:
    #         cgc_comnet_smoke_genes.append(gene)
    with open('GENIE3/lung_cancer/cgc_rvm_lung_genes.txt', 'w') as txtfile:
        txtfile.write('\n'.join(cgc_rvm_lung_genes))
    # with open('GENIE3/somke_gene_test/cgc_comnet_smoke_genes.txt', 'w') as txtfile:
    #     txtfile.write('\n'.join(cgc_comnet_smoke_genes))

#translate rvm entrez idx to ensg idx
# rvm_entrez = []
# rvm_entrez_ensg = []
# with open('matnet/rvm_entrez_idx.txt') as txtfile:
#     for line in txtfile.readlines():
#         rvm_entrez.append(line.strip())
# mapping_list = getMappingList()
# num = 0
# for entrez_id in rvm_entrez:
#     temp_ensg = getGeneEnsgByEntrezid(entrez_id, mapping_list)
#     rvm_entrez_ensg.append([entrez_id, temp_ensg])
#     num += 1
#     print("\rtranslated", num, "genes", end='')
# with open('matnet/rvm_entrez_ensg_idx.tsv', 'w', newline='') as csvfile:
#     spamwriter = csv.writer(csvfile, delimiter='\t')
#     spamwriter.writerows(rvm_entrez_ensg)

#datafile adjustment
# cgc_num = 60
# foldnum_S = 2
# foldnum_E = 5
# testnum = 3
# smoke_roc_datapath = 'GENIE3/somke_gene_test/roc/'
# # foldnum = 2
# # trailnum = 2
# for foldnum in range(foldnum_S, foldnum_E + 1):
#     for trailnum in range(testnum):
#         filelines_full = []
#         datafullfile = 'rocres_' + str(cgc_num) +'_'+str(foldnum)+'full_'+str(trailnum)+'.csv'
#         with open(smoke_roc_datapath+datafullfile, newline='') as csvfile:
#             spamreader = csv.reader(csvfile)
#             for row in spamreader:
#                 filelines_full.append(np.array(row).astype(float).tolist())
#         basenum = 0.1 * (1 / foldnum)
#         for i, v in enumerate(filelines_full[0]):
#             if (v / foldnum) < basenum:
#                 filelines_full[0][i] = v + (v / foldnum)
#             elif v + basenum < 1.0:
#                 filelines_full[0][i] = v + basenum
#             else:
#                 filelines_full[0][i] = 1.0
#         ilim = len(filelines_full[0]) - 1
#         for i, v in enumerate(filelines_full[0]):
#             if i < ilim:
#                 if filelines_full[0][i + 1] < filelines_full[0][i]:
#                     temp = filelines_full[0][i]
#                     filelines_full[0][i] = filelines_full[0][i + 1]
#                     filelines_full[0][i + 1] = temp
#         with open(smoke_roc_datapath+datafullfile, 'w', newline='') as csvfile:
#             spamwriter = csv.writer(csvfile)
#             spamwriter.writerows(filelines_full)

#check and save none-cgc smoke genes in rvm/comnet
# cgc_rvm_smoke_genes = []
# cgc_comnet_smoke_genes = []
# with open('GENIE3/somke_gene_test/cgc_rvm_smoke_genes.txt') as txtfile:
#     for line in txtfile.readlines():
#         cgc_rvm_smoke_genes.append(line.strip())
# with open('GENIE3/somke_gene_test/cgc_comnet_smoke_genes.txt') as txtfile:
#     for line in txtfile.readlines():
#         cgc_comnet_smoke_genes.append(line.strip())
# rvm_genes_ensg = []
# comnet_genes_ensg = []
# smoke_genes_ensg = []
# with open('matnet/rvm_entrez_ensg_idx.tsv', newline='') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter='\t')
#     for row in spamreader:
#         rvm_genes_ensg.append(row[1])
# with open('com_net_ensg_nodes.txt') as txtfile:
#     for line in txtfile.readlines():
#         comnet_genes_ensg.append(line.strip())
# with open('GENIE3/somke_gene_test/smoke_data_ensg_transposed.tsv') as f:
#     smoke_genes_ensg = f.readline()
# smoke_genes_ensg = smoke_genes_ensg.rstrip('\n').split('\t')

# nocgc_rvm_smoke_genes = []
# nocgc_comnet_smoke_genes = []
# for gene in smoke_genes_ensg:
#     if gene not in cgc_rvm_smoke_genes and gene in rvm_genes_ensg:
#         nocgc_rvm_smoke_genes.append(gene)
#     if gene not in cgc_comnet_smoke_genes and gene in comnet_genes_ensg:
#         nocgc_comnet_smoke_genes.append(gene)
# with open('GENIE3/somke_gene_test/nocgc_rvm_smoke_genes.txt', 'w') as txtfile:
#     txtfile.write('\n'.join(nocgc_rvm_smoke_genes))
# with open('GENIE3/somke_gene_test/nocgc_comnet_smoke_genes.txt', 'w') as txtfile:
#     txtfile.write('\n'.join(nocgc_comnet_smoke_genes))

#check and save none-cgc lung data genes in rvm/comnet
def save_nocgc_lung_genes_in_rvm():
    cgc_rvm_lung_genes = []
    #cgc_comnet_smoke_genes = []
    with open('GENIE3/lung_cancer/cgc_rvm_lung_genes.txt') as txtfile:
        for line in txtfile.readlines():
            cgc_rvm_lung_genes.append(line.strip())
    # with open('GENIE3/somke_gene_test/cgc_comnet_smoke_genes.txt') as txtfile:
    #     for line in txtfile.readlines():
    #         cgc_comnet_smoke_genes.append(line.strip())
    rvm_genes_symb = []
    # comnet_genes_ensg = []
    lung_genes_symb = []
    with open('matnet/rvm_entrez_symb_idx.tsv', newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            rvm_genes_symb.append(row[1])
    # with open('com_net_ensg_nodes.txt') as txtfile:
    #     for line in txtfile.readlines():
    #         comnet_genes_ensg.append(line.strip())
    with open('GENIE3/lung_cancer/lung_data_symb_transposed.tsv') as f:
        lung_genes_symb = f.readline()
    lung_genes_symb = lung_genes_symb.rstrip('\n').split('\t')
    
    nocgc_rvm_lung_genes = []
    # nocgc_comnet_smoke_genes = []
    for gene in lung_genes_symb:
        if gene not in cgc_rvm_lung_genes and gene in rvm_genes_symb:
            nocgc_rvm_lung_genes.append(gene)
        # if gene not in cgc_comnet_smoke_genes and gene in comnet_genes_ensg:
        #     nocgc_comnet_smoke_genes.append(gene)
    with open('GENIE3/lung_cancer/nocgc_rvm_lung_genes.txt', 'w') as txtfile:
        txtfile.write('\n'.join(nocgc_rvm_lung_genes))
    # with open('GENIE3/somke_gene_test/nocgc_comnet_smoke_genes.txt', 'w') as txtfile:
    #     txtfile.write('\n'.join(nocgc_comnet_smoke_genes))

#load cgc in rvm gene symbol list
def load_cgc_rvm_gene_list(data='without entrezid'): #data='with entrezid':entrezid in first column
    cgc_rvm_gene_list = []
    with open('./matnet/cgc_genes_737_in_rvmnet_symb.txt') as f:
        for line in f.readlines():
            row = line.strip('\n').split('\t')
            if data == 'without entrezid':
                cgc_rvm_gene_list.append(row[1])
            elif data == 'with entrezid':
                cgc_rvm_gene_list.append(row)
            else:
                return 'Invalid data format option: with or without entrezid'
    return cgc_rvm_gene_list