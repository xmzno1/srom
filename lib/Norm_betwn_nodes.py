import numpy as np

#weight_array = np.arange(0.1, 4.1, 0.1)
#C_W = 1.0
#P_W = 1.0

# com_array_pruned = np.load('GENIE3/com_array_pruned_d.npy')
# per_array_pruned = np.load('GENIE3/per_array_pruned_d.npy')
def calMergedArray(com_array_pruned, per_array_pruned, savePath='GENIE3/merged_net_adj_array_d.npy'):
    merge_edge_weights = []
    node_num = len(com_array_pruned)
    merged_array = np.zeros((node_num, node_num))
    # com_merge_weights_array = np.ones(node_num)
    # per_merge_weights_array = np.ones(node_num)
    
    for node_idx in range(node_num):
        com_node = com_array_pruned[node_idx]
        per_node = per_array_pruned[node_idx]
        
        norm_com_node = np.linalg.norm(com_node, ord = 2)
        norm_per_node = np.linalg.norm(per_node, ord = 2)
        if norm_com_node != 0.0 and norm_per_node != 0.0:
            if norm_com_node > norm_per_node:
                merge_edge_weights.append([1.0, norm_com_node / norm_per_node])
            elif norm_per_node > norm_com_node:
                merge_edge_weights.append([norm_per_node / norm_com_node, 1.0])
        else:
            merge_edge_weights.append([1.0, 1.0])
        merged_array[node_idx] = (com_node * merge_edge_weights[node_idx][0]) + (per_node * merge_edge_weights[node_idx][1])
        # if np.linalg.norm(com_node, ord = 2) != 0.0 and np.linalg.norm(per_node, ord = 2) != 0.0:
        #     diff_min = float('inf')
        #     temp_com_node_merge_weight = 0.0
        #     temp_per_node_merge_weight = 0.0
        #     for i in range(len(weight_array)):
        #         temp_com_node_merge_weight = weight_array[i]
        #         for j in range(len(weight_array)):
        #             temp_per_node_merge_weight = weight_array[j]
        #             temp_com_node_norm = np.linalg.norm(com_node * temp_com_node_merge_weight, ord = 2)
        #             temp_per_node_norm = np.linalg.norm(per_node * temp_per_node_merge_weight, ord = 2)
        #             if abs(temp_com_node_norm - temp_per_node_norm) < diff_min:
        #                 diff_min = abs(temp_com_node_norm - temp_per_node_norm)
        #                 com_merge_weights_array[node_idx] = temp_com_node_merge_weight
        #                 per_merge_weights_array[node_idx] = temp_per_node_merge_weight
        print("\rcalculated", node_idx+1, "//", node_num, "nodes", end='')
    
    if savePath is not None:
        np.save(savePath, merged_array, allow_pickle=False, fix_imports=False)
    return merged_array, merge_edge_weights
#np.save('com_merge_weights_array.npy', com_merge_weights_array, allow_pickle=False, fix_imports=False)
#np.save('per_merge_weights_array.npy', per_merge_weights_array, allow_pickle=False, fix_imports=False)