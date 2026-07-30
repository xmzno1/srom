import lib.Processing as pr
import lib.MI_betwn_nodes as mi
import lib.Node_edges_pruning as ep
import lib.Norm_betwn_nodes as nbn
import lib.GENIE3 as g3
import numpy as np
#import matplotlib.pyplot as plt

def get_edgelist_by_value(value_arr, edge_value):
    indices = np.where(value_arr == edge_value)
    reslist = []
    for i in range(len(indices[0])):
        reslist.append([indices[0][i], indices[1][i]])
    return reslist

def test(pred, target, threshold=0.5):
    pred1d = pred.reshape((len(pred)*len(pred)))
    target1d = target.reshape((len(target)*len(target)))
    ''' sort edge values in descending order '''
    pred1d_sort = np.sort(pred1d)[::-1]
    pred1d_sort_nozero = np.delete(pred1d_sort, np.where(pred1d_sort == 0.0))
    target1d_sort = np.sort(target1d)[::-1]
    target1d_sort_nozero = np.delete(target1d_sort, np.where(target1d_sort == 0.0))
    
    pred1d_nodup = np.unique(pred1d_sort_nozero)
    target1d_nodup = np.unique(target1d_sort_nozero)
    threshold_num = int(threshold*len(target1d_nodup))
    
    target_edges = []
    for edge_value in target1d_nodup[:threshold_num]:
        tmp = get_edgelist_by_value(target, edge_value)
        target_edges.extend(tmp)
    ''' calculate tprs & fprs '''
    value_num = len(pred1d_nodup)
    tprs = np.zeros(value_num - 1)
    fprs = np.zeros(value_num - 1)
    precision = np.zeros(value_num - 1)
    for i in range(1, value_num - 1):
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        for key, value in enumerate(pred1d_nodup):
            edges = get_edgelist_by_value(pred, value)
            for edge in edges:
                if edge in target_edges:
                    if key < i:
                        tp += 1
                    else:
                        fn += 1
                else:
                    if key < i:
                        fp += 1
                    else:
                        tn += 1
        tprs[i] = tp / (tp + fn) #recall
        fprs[i] = fp / (fp + tn)
        precision[i] = tp / (tp + fp)
    return fprs,tprs,precision

def srom_process(ori_grn):
    rdata = np.load('gene_rel.npy')
    vim_fold_u, vim_fold_d = pr.getUndirectedGraphAM(ori_grn)
    merge_param_index_array_u = mi.calMergeParam(rdata, vim_fold_u, 
                                                 savepath=None)
    merge_param_index_array_d = mi.calMergeParam(rdata, vim_fold_d, 
                                                 savepath=None)
    
    com_array_pruned_u, per_array_pruned_u = ep.edgePruning(rdata, vim_fold_u, merge_param_index_array_u, save_com_path=None, save_per_path=None)
    com_array_pruned_d, per_array_pruned_d = ep.edgePruning(rdata, vim_fold_d, merge_param_index_array_d, save_com_path=None, save_per_path=None)
    
    merged_array_u, merge_edge_weights_u = nbn.calMergedArray(com_array_pruned_u, per_array_pruned_u, savePath=None)
    merged_array_d, merge_edge_weights_d = nbn.calMergedArray(com_array_pruned_d, per_array_pruned_d, savePath=None)
    
    merged_net = pr.combineUndirectedGraphAM(merged_array_u, 
                                             merged_array_d)
    return merged_net

'''-----------------------ROC & AUC test----------------------------------'''
def test_10gene_crossval_edgematch(run_count = 30, outfileT='expression_data.txt', 
                                   nthreads=1, threshold=0.5):
    res_p3 = dict()
    res_p3['fprs'] = []
    res_p3['tprs'] = []
    res_p3['precision'] = []
    res_p5 = dict()
    res_p5['fprs'] = []
    res_p5['tprs'] = []
    res_p5['precision'] = []
    res_p7 = dict()
    res_p7['fprs'] = []
    res_p7['tprs'] = []
    res_p7['precision'] = []
    
    for i in range(run_count):
        print('Run count: ', str(i), 'start.')
        data_full, gene_names = pr.loadData(outfileT)
        np.random.shuffle(data_full)
        row_num = data_full.shape[0]
        
        data_p3 = data_full[:int(row_num*0.25)]
        data_p5 = data_full[:int(row_num*0.5)]
        data_p7 = data_full[:int(row_num*0.75)]
        
        target_vim = g3.GENIE3(data_full, nthreads=nthreads)
        
        vim_p3 = g3.GENIE3(data_p3, nthreads=nthreads)
        vim_p5 = g3.GENIE3(data_p5, nthreads=nthreads)
        vim_p7 = g3.GENIE3(data_p7, nthreads=nthreads)
        
        ''' merging process '''
        pred_p3 = srom_process(vim_p3)
        pred_p5 = srom_process(vim_p5)
        pred_p7 = srom_process(vim_p7)
        
        ''' crossvalidation '''
        fprs,tprs,precision = test(pred_p3, target_vim, threshold=threshold)
        res_p3['fprs'].append(fprs)
        res_p3['tprs'].append(tprs)
        res_p3['precision'].append(precision)
        
        fprs,tprs,precision = test(pred_p5, target_vim, threshold=threshold)
        res_p5['fprs'].append(fprs)
        res_p5['tprs'].append(tprs)
        res_p5['precision'].append(precision)
        
        fprs,tprs,precision = test(pred_p7, target_vim, threshold=threshold)
        res_p7['fprs'].append(fprs)
        res_p7['tprs'].append(tprs)
        res_p7['precision'].append(precision)
    
    return res_p3, res_p5, res_p7
        

if __name__ == '__main__': #this line is required for multi-process computing
    run_count = 3
    nthreads = 1
    threshold = 0.5
    res_p3, res_p5, res_p7 = test_10gene_crossval_edgematch(run_count = run_count, 
                                                            outfileT='expression_data.txt', 
                                                            nthreads=nthreads, 
                                                            threshold=threshold)
    aucs = np.zeros((3, run_count))
    for count in range(run_count):
        for key, res in enumerate([res_p3, res_p5, res_p7]):
            roc = np.vstack((res['fprs'][count], res['tprs'][count]))
            auc = 0.0
            for i in range(roc.shape[1]):
                if i != 0:
                    temp = roc[1][i] * (roc[0][i] - roc[0][i - 1])
                else:
                    temp = roc[1][i] * roc[0][i]
                auc += temp
            aucs[key][count] = auc
    
    print('\n=== SROM GRN Edges Prediction CV Results: ===')
    print('Test trail\t|\tPred on 25% data vs Truth\t|\tPred on 50% data vs Truth\t|\tPred on 75% data vs Truth')
    for count in range(run_count):
        print(f'Test #{count} \t | \t {aucs[0][count]:.2f} \t | \t  {aucs[1][count]:.2f} \t | \t {aucs[2][count]:.2f}')
