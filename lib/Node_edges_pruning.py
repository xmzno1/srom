import numpy as np
#import copy



# merge_param_index_array = np.load('GENIE3/merge_param_index_array_d.npy')
# # merge_param = []
# # for i in range(len(merge_param_index_array)):
# #     merge_param.append([merge_param_index_array[i][0] * 0.1, merge_param_index_array[i][1] * 0.1])

# C_W = 1.0
# P_W = 1.0
# com_array = np.load('GENIE3/rvmdata.npy')
# per_array = np.load('GENIE3/vim_fold_d.npy')
def edgePruning(com_array, per_array, merge_param_index_array, save_com_path='GENIE3/com_array_pruned_d.npy', save_per_path='GENIE3/per_array_pruned_d.npy'):
    node_num = len(com_array)
    com_array_pruned = np.array([])
    per_array_pruned = np.array([])
    com_pruned_list = []
    per_pruned_list = []
    for node_idx in range(node_num):
        #sys.stdout.write('\033[2K\033[1G')
        com_node = com_array[node_idx]
        per_node = per_array[node_idx]
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
        
        
        for c_idx in range(node_num):
            if com_node[c_idx] < com_bin_edge[int(merge_param_index_array[node_idx][0])] and com_prune_flag == True:
                com_node[c_idx] = 0
    
        for p_idx in range(node_num):
            if per_node[p_idx] < com_bin_edge[int(merge_param_index_array[node_idx][1])] and per_prune_flag == True:
                per_node[p_idx] = 0
                
        # if node_idx == 0:
        #     com_array_pruned = com_node
        #     per_array_pruned = per_node
        # else:
        #     com_array_pruned = np.vstack((com_array_pruned, com_node))
        #     per_array_pruned = np.vstack((per_array_pruned, per_node))
        
        # increase running speed using python list, may cause large memory usage
        com_pruned_list.append(com_node.tolist())
        per_pruned_list.append(per_node.tolist())
        print("\rcalculated", node_idx+1, "//", node_num, "nodes", end='')
    
    com_array_pruned = np.array(com_pruned_list)
    per_array_pruned = np.array(per_pruned_list)
    if save_com_path is not None:
        np.save(save_com_path, com_array_pruned, allow_pickle=False, fix_imports=False)
    if save_per_path is not None:
        np.save(save_per_path, per_array_pruned, allow_pickle=False, fix_imports=False)
    return com_array_pruned,per_array_pruned