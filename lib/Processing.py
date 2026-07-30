import csv
import random as rand
import copy
import numpy as np
import scipy.io as sio
import lib.Mat_op as mat

def load_txt_to_list(filepath, filetype = 'tsv'):
    res = []
    if filetype == 'csv':
        with open(filepath, newline='') as csvfile:
            filereader = csv.reader(csvfile)
            for row in filereader:
                if len(row) == 1:
                    res.append(row[0])
                else:
                    res.append(row)
    else:
        with open(filepath) as f:
            for line in f.readlines():
                row = line.strip('\n').split('\t')
                if len(row) == 1:
                    res.append(row[0])
                else:
                    res.append(row)
    return res

def write_list_to_file(inputlist, outpath, outtype='tsv'):
    file_out = open(outpath, 'w', newline='')
    if outtype == 'tsv':
        spamwriter = csv.writer(file_out, delimiter='\t')
        spamwriter.writerows(inputlist)
    elif outtype == 'csv':
        spamwriter = csv.writer(file_out)
        spamwriter.writerows(inputlist)
    else:
        print("Input file type not supported")
        return False

def transpose_tsv_file(filepath, outpath, filetype='tsv', outtype='tsv'):
    file_in = open(filepath, newline='')
    if filetype == 'csv':
        a = zip(*csv.reader(file_in))
    elif filetype == 'tsv':
        a = zip(*csv.reader(file_in, delimiter='\t'))
    else:
        print("Input file type not supported")
        return False
    write_list_to_file(a, outpath, outtype)

def sort_dict_by_value(data_d, reverse = True):
    sorted_res_by_value = sorted(data_d.items(), key=lambda x:x[1], reverse=reverse)
    converted_dict = dict(sorted_res_by_value)
    return converted_dict

#get gene expression data
def get_gene_exp_data(ori_data_file, ori_data_file_type, ori_gene_name_file,
                      cgc_gene_file, cgc_num_portion = 0.3, nocgc_ratio = 2.0,
                      outfilepath='./', outfilesufix=''):
    #global cgc_num
    #cgc_num = 50
    data_ori_list = load_txt_to_list(ori_data_file, filetype=ori_data_file_type)
    gene_ori_list = load_txt_to_list(ori_gene_name_file)
    cgc_genes = load_txt_to_list(cgc_gene_file)
    nocgc_genes = copy.deepcopy(gene_ori_list)
    for gene in cgc_genes:
        nocgc_genes.remove(gene)
    cgc_num = int(len(cgc_genes) * cgc_num_portion)
    nocgc_num = int(cgc_num * nocgc_ratio)
    rand.seed()
    genes_sampled = rand.sample(cgc_genes, cgc_num)
    genes_sampled.extend(rand.sample(nocgc_genes, nocgc_num))
    data_sampled = []
    for gene in genes_sampled:
        data_sampled.append(data_ori_list[gene_ori_list.index(gene)])
    outfile = outfilepath + 'data_sampled_' + outfilesufix + str(cgc_num)+'_'+str(nocgc_num)+ '.tsv'
    write_list_to_file(data_sampled, outfile)
    outfileT = outfilepath + 'data_sampled_transposed_' + outfilesufix + str(cgc_num)+'_'+str(nocgc_num)+ '.tsv'
    transpose_tsv_file(outfile, outfileT)
    return outfileT, cgc_num, nocgc_num

def loadData(file, delimiter='\t'):
    data_full = np.loadtxt(file, delimiter=delimiter, skiprows=1)
    f = open(file)
    gene_names = f.readline()
    f.close()
    gene_names = gene_names.rstrip('\n').split('\t')
    return data_full, gene_names

def getFoldData(data_full, foldnum=5):
    data_full_len = len(data_full)
    folded_len = data_full_len - int(data_full_len/foldnum)
    data_folded = np.zeros((folded_len, data_full.shape[1]))
    arr = np.arange(data_full_len)
    rng = np.random.default_rng()
    rng.shuffle(arr)
    for i in range(folded_len):
        data_folded[i] = data_full[arr[i]]
    return data_folded

def getUndirectedGraphAM(d_graph_adjmatrix):
    g_num = len(d_graph_adjmatrix)
    upside_adjmatrix = np.zeros((g_num, g_num))
    downside_adjmatrix = np.zeros((g_num, g_num))
    d_graph_adjmatrix_T = np.transpose(d_graph_adjmatrix)
    for i in range(g_num):
        upside_adjmatrix[i][:i] = d_graph_adjmatrix_T[i][:i]
        upside_adjmatrix[i][i:] = d_graph_adjmatrix[i][i:]
        downside_adjmatrix[i][:i] = d_graph_adjmatrix[i][:i]
        downside_adjmatrix[i][i:] = d_graph_adjmatrix_T[i][i:]
    return upside_adjmatrix, downside_adjmatrix

def combineUndirectedGraphAM(upside, downside):
    g_num = len(upside)
    d_graph_adjmatrix = np.zeros((g_num, g_num))
    for i in range(g_num):
        d_graph_adjmatrix[i][i:] = upside[i][i:]
        d_graph_adjmatrix[i][:i] = downside[i][:i]
    return d_graph_adjmatrix

def crossval_match_rocauc(data_pred, data_match, cgc_num, gene_names):
    gene_num = len(data_pred)
    if gene_num != len(data_match):
        return "Can't match!"
    m_node_cen_dict = dict()
    p_node_cen_dict = dict()
    for i in range(gene_num):
        m_node_cen_dict[gene_names[i]] = data_match[i].sum() + data_match[:,i].sum()
        p_node_cen_dict[gene_names[i]] = data_pred[i].sum() + data_pred[:,i].sum()
    m_node_cen_dict_sorted = sort_dict_by_value(m_node_cen_dict)
    p_node_cen_dict_sorted = sort_dict_by_value(p_node_cen_dict)
    pred_genes_sorted = list(p_node_cen_dict_sorted.keys())
    match_genes_sorted = list(m_node_cen_dict_sorted.keys())
    tprs = np.zeros(gene_num)
    fprs = np.zeros(gene_num)
    for i in range(1, gene_num):
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        for idx, gene in enumerate(pred_genes_sorted):
            if idx < i:
                #if gene in match_genes_sorted[:int(0.5*gene_num)]:
                if gene in gene_names[:cgc_num]:
                    tp += 1
                else:
                    fp += 1
            else:
                #if gene in match_genes_sorted[int(0.5*gene_num):]:
                if gene in gene_names[cgc_num:]:
                    tn += 1
                else:
                    fn += 1
        tprs[i] = tp / (tp + fn)
        fprs[i] = fp / (fp + tn)
        # tprs[i] = tp / len(match_genes_sorted[:i])
        # fprs[i] = fp / len(match_genes_sorted[i:])
    roc = np.vstack((fprs, tprs))
    auc = 0.0
    for i in range(roc.shape[1]):
        if i != 0:
            temp = roc[1][i] * (roc[0][i] - roc[0][i - 1])
        else:
            temp = roc[1][i] * roc[0][i]
        auc += temp
    return roc, auc

def cgc_pred_rocauc_by_nodelist(pred_node_list, cgc_list):
    gene_num = len(pred_node_list)
    tprs = np.zeros(gene_num)
    fprs = np.zeros(gene_num)
    for i in range(1, gene_num):
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        for idx, gene in enumerate(pred_node_list):
            if idx < i:
                if gene in cgc_list:
                    tp += 1
                else:
                    fp += 1
            else:
                if gene not in cgc_list:
                    tn += 1
                else:
                    fn += 1
        tprs[i] = tp / (tp + fn)
        fprs[i] = fp / (fp + tn)
    roc = np.vstack((fprs, tprs))
    auc = 0.0
    for i in range(roc.shape[1]):
        if i != 0:
            temp = roc[1][i] * (roc[0][i] - roc[0][i - 1])
        else:
            temp = roc[1][i] * roc[0][i]
        auc += temp
    return roc, auc

def crossval_match_f1(data_pred, data_match, cgc_num, gene_names, grate=2.0):#2.15
    gene_num = len(data_pred)
    if gene_num != len(data_match):
        return "Can't match!"
    m_node_cen_dict = dict()
    p_node_cen_dict = dict()
    for i in range(gene_num):
        m_node_cen_dict[gene_names[i]] = data_match[i].sum() + data_match[:,i].sum()
        p_node_cen_dict[gene_names[i]] = data_pred[i].sum() + data_pred[:,i].sum()
    m_node_cen_dict_sorted = sort_dict_by_value(m_node_cen_dict)
    p_node_cen_dict_sorted = sort_dict_by_value(p_node_cen_dict)
    pred_genes_sorted = list(p_node_cen_dict_sorted.keys())
    match_genes_sorted = list(m_node_cen_dict_sorted.keys())
    threshold = cgc_num / gene_num
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for idx, gene in enumerate(pred_genes_sorted):
        if idx < int(threshold * gene_num * grate):
            if gene in match_genes_sorted[:int(threshold * gene_num * grate)]:
                tp += 1
            else:
                fp += 1
        else:
            if gene in match_genes_sorted[int(threshold * gene_num * grate):]:
                tn += 1
            else:
                fn += 1
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    beta = 1.0
    f1 = (1+beta**2) * (precision*recall) / (precision*(beta**2)+recall)
    print("f1 = ", str(f1))
    return precision, recall, f1

def getSromGeneData(gene_names, geneNameType='Ensg', quickmode=True):
    gene_num = len(gene_names)
    rdata = np.zeros((gene_num, gene_num))
    rvmIdxSymbolMappingList = mat.getRvmIdxSymbolMappingList()
    rvmIdxMappingList = mat.getRvmIdxMappingList()
    rvmidxs = []
    if geneNameType == "Symbol":
        for gene_name in gene_names:
            rvmidxs.append(mat.getGeneRvmidxBySymbol(gene_name, rvmIdxSymbolMappingList))
    elif geneNameType == "Ensg":
        for gene_name in gene_names:
            rvmidxs.append(mat.getGeneRvmidxByEnsg(gene_name, rvmIdxMappingList))
    elif geneNameType == "Entrezid":
        for gene_name in gene_names:
            rvmidxs.append(mat.getGeneRvmidxByEntrezid(gene_name, rvmIdxMappingList))
            print("\r", gene_name, end='')
    if quickmode:
        netmat = np.load('matnet/rvm_data.npy')
        for i, rvmidx in enumerate(rvmidxs):
            if rvmidx == 'False':
                continue
            for j, edge_idx in enumerate(rvmidxs):
                if edge_idx == 'False':
                    continue
                rdata[i,j] = netmat[rvmidx,edge_idx]
            print("\rprepared", i+1, "//", gene_num, "genes data", end='')
    else:
        for i, rvmidx in enumerate(rvmidxs):
            if rvmidx == 'False':
                continue
            filename = "matnet/mat/yrvm_exp_"+str(int(rvmidx/670 + 1))+".mat"
            mat_content = sio.loadmat(filename)
            for j, edge_idx in enumerate(rvmidxs):
                if edge_idx == 'False':
                    continue
                rdata[i][j] = mat_content['yrvm_exp'][rvmidx%670][edge_idx]
            print("\rprepared", i+1, "//", gene_num, "genes data", end='')
    return rdata

def getRvmGeneData(gene_names, geneNameType='Ensg', quickmode=True):
    gene_num = len(gene_names)
    rdata = np.zeros((gene_num, gene_num))
    rvmIdxSymbolMappingList = mat.getRvmIdxSymbolMappingList()
    rvmIdxMappingList = mat.getRvmIdxMappingList()
    rvmidxs = []
    if geneNameType == "Symbol":
        for gene_name in gene_names:
            rvmidxs.append(mat.getGeneRvmidxBySymbol(gene_name, rvmIdxSymbolMappingList))
    elif geneNameType == "Ensg":
        for gene_name in gene_names:
            rvmidxs.append(mat.getGeneRvmidxByEnsg(gene_name, rvmIdxMappingList))
    elif geneNameType == "Entrezid":
        for gene_name in gene_names:
            rvmidxs.append(mat.getGeneRvmidxByEntrezid(gene_name, rvmIdxMappingList))
            print("\r", gene_name, end='')
    if quickmode:
        netmat = np.load('matnet/rvm_data.npy')
        for i, rvmidx in enumerate(rvmidxs):
            if rvmidx == 'False':
                continue
            for j, edge_idx in enumerate(rvmidxs):
                if edge_idx == 'False':
                    continue
                rdata[i,j] = netmat[rvmidx,edge_idx]
            print("\rprepared", i+1, "//", gene_num, "genes data", end='')
    else:
        for i, rvmidx in enumerate(rvmidxs):
            if rvmidx == 'False':
                continue
            filename = "matnet/mat/yrvm_exp_"+str(int(rvmidx/670 + 1))+".mat"
            mat_content = sio.loadmat(filename)
            for j, edge_idx in enumerate(rvmidxs):
                if edge_idx == 'False':
                    continue
                rdata[i][j] = mat_content['yrvm_exp'][rvmidx%670][edge_idx]
            print("\rprepared", i+1, "//", gene_num, "genes data", end='')
    return rdata