import numpy as np
from sklearn.metrics import mutual_info_score
import copy
#import json
#import sys

def get_discrete_values(continuous_value_array):
    discrete_value_list = []
    for i in range(len(continuous_value_array)):
        if continuous_value_array[i] > 0:
            discrete_value_list.append(1)
        else:
            discrete_value_list.append(0)
    return np.array(discrete_value_list)

def normalize_minmax(arr, t_min, t_max):
    arr_1d = arr.flatten()
    norm_arr = []
    diff = t_max - t_min
    diff_arr = max(arr_1d) - min(arr_1d)
    for i in arr_1d:
        if diff_arr == 0:
            temp = t_max
        else:
            temp = (((i - min(arr_1d))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return np.array(norm_arr).reshape(arr.shape)

# C_W = 1.0
# P_W = 1.0
# # com_array = np.load('com_array.npy')
# # per_array = np.load('per_array.npy')
# com_array = np.load('GENIE3/rvmdata.npy')
# per_array = np.load('GENIE3/vim_fold_d.npy')
def calMergeParam(com_array, per_array, savepath='GENIE3/merge_param_index_array_d.npy'):
    node_num = len(com_array)
    merge_param_index = []
    merge_param_index_array = np.zeros((node_num, 2))
    #merge_param = []
    
    for node_idx in range(node_num):
        #sys.stdout.write('\033[2K\033[1G')
        com_node = com_array[node_idx]
        per_node = per_array[node_idx]
        discrete_com_node = get_discrete_values(com_node)
        discrete_per_node = get_discrete_values(per_node)
        #mi_com_per = np.zeros((10, 10))
        #mi_merge_com = np.zeros((10, 10))
        #mi_merge_per = np.zeros((10, 10))
        #mi_total = np.zeros((10, 10))
        #com_hist, com_bin_edge = np.histogram(com_node)
        #per_hist, per_bin_edge = np.histogram(per_node)
        
        # check edge number
        com_nozero_count = 0
        per_nozero_count = 0
        com_node_edges = []
        per_node_edges = []
        for i in range(node_num):
            if com_node[i] > 0:
                com_nozero_count += 1
                com_node_edges.append(com_node[i])
            if per_node[i] > 0:
                per_nozero_count += 1
                per_node_edges.append(per_node[i])
        com_prune_flag = True
        per_prune_flag = True
        prune_lvl = 0.001
        if com_nozero_count <= node_num * prune_lvl:
            com_prune_flag = False
        if per_nozero_count <= node_num * prune_lvl:
            per_prune_flag = False
        pruned_com_node = copy.deepcopy(com_node)
        pruned_per_node = copy.deepcopy(per_node)
        
        com_node_edges.sort()
        per_node_edges.sort()
        com_bin_length = int(com_nozero_count / 10)
        per_bin_length = int(per_nozero_count / 10)
        com_bin_edge = [0]
        per_bin_edge = [0]
        for i in range(1, 10):
            if com_prune_flag:
                com_bin_edge.append(com_node_edges[com_bin_length * i])
            else:
                com_bin_edge.append(0)
            if per_prune_flag:
                per_bin_edge.append(per_node_edges[per_bin_length * i])
            else:
                per_bin_edge.append(0)  
        
        # if com_nozero_count == 0 or per_nozero_count == 0:
        #     merge_param_index.append([0, 0])
        # else:
        if com_nozero_count != 0 and per_nozero_count != 0:
            mi_max = 0.0
            for i in range(10):
                for c_idx in range(node_num):
                    if pruned_com_node[c_idx] < com_bin_edge[i] and com_prune_flag == True:
                        pruned_com_node[c_idx] = 0
                for j in range(10):
                    for p_idx in range(node_num):
                        if pruned_per_node[p_idx] < per_bin_edge[j] and per_prune_flag == True:
                            pruned_per_node[p_idx] = 0
                
                    temp_discrete_pruned_com_node = get_discrete_values(pruned_com_node)
                    temp_discrete_pruned_per_node = get_discrete_values(pruned_per_node)
                    #temp_merge_node = (C_W * temp_discrete_pruned_com_node) + (P_W * temp_discrete_pruned_per_node)
                    #temp_discrete_merge_node = get_discrete_values(temp_merge_node)
                    
                    temp_mi_com_per = mutual_info_score(temp_discrete_pruned_com_node, temp_discrete_pruned_per_node)
                    if temp_mi_com_per > mi_max:
                        mi_max = temp_mi_com_per
                        merge_param_index_array[node_idx][0] = i
                        merge_param_index_array[node_idx][1] = j
                    #temp_mi_merge_com = mutual_info_score(temp_discrete_merge_node, discrete_com_node)
                    #temp_mi_merge_per = mutual_info_score(temp_discrete_merge_node, discrete_per_node)
                    #mi_com_per[i][j] = temp_mi_com_per
                    # mi_merge_com[i][j] = temp_mi_merge_com
                    # mi_merge_per[i][j] = temp_mi_merge_per
                    #mi_total[i][j] = temp_mi_per_com + temp_mi_merge_com + temp_mi_merge_per
                # mi_per_com[i] = temp_mi_per_com
                # mi_merge_com[i] = temp_mi_merge_com
                # mi_merge_per[i] = temp_mi_merge_per
                # mi_total[i] = temp_mi_per_com + temp_mi_merge_com + temp_mi_merge_per
            
            # normalized_mi_com_per = normalize_minmax(mi_com_per, 1, 10)
            # com_cut_idx, per_cut_idx = np.where(np.isclose(normalized_mi_com_per, 10.0))
            # merge_param_index.append([com_cut_idx[0], per_cut_idx[0]])
        #print("\033[2K\033[1G")
        print("\rcalculated", node_idx+1, "//", node_num, "nodes", end='')
    
    #merge_param_index_array = np.array(merge_param_index)
    if savepath is not None:
        np.save(savepath, merge_param_index_array, allow_pickle=False, fix_imports=False)
    return merge_param_index_array
# normalized_mi_merge_com = normalize_minmax(mi_merge_com, 1, 10)
# normalized_mi_merge_per = normalize_minmax(mi_merge_per, 1, 10)
# for i in range(10):
#     for j in range(10):
#         mi_total[i][j] = normalized_mi_per_com[i][j] + normalized_mi_merge_com[i][j] + normalized_mi_merge_per[i][j]

# for i in range(10):
#     mi_total[i] = normalized_mi_per_com[i] + normalized_mi_merge_com[i] + normalized_mi_merge_per[i]

#merge_node = (C_W * com_node) + (P_W * per_node)

#mi_per_com = mutual_info_score(per_node, com_node)
#mi_merge_com = mutual_info_score(merge_node, com_node)
#mi_merge_per = mutual_info_score(merge_node, per_node)

