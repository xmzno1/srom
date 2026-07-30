import numpy as np
import networkx as nx

def create_weighted_graph_from_adjacency(node_names, adjacency_matrix, weight_threshold=0.5, create_directed=False):
    """
    根据节点名称列表和邻接矩阵创建带权重的networkx图，并过滤低权重边
    
    参数:
        node_names: list，节点名称列表，长度应与邻接矩阵维度一致
        adjacency_matrix: numpy.array，邻接矩阵，应为方阵
        weight_threshold: float，权重阈值，低于此值的边将被过滤掉，默认为0.5
        create_directed: bool，是否创建有向图，默认为False（创建无向图）
        
    返回:
        networkx.Graph 或 networkx.DiGraph: 带权重的图对象
        
    异常:
        ValueError: 当输入参数不匹配时抛出
    """
    
    # 参数验证
    if not isinstance(node_names, list):
        raise ValueError("node_names 必须是列表类型")
    
    if not isinstance(adjacency_matrix, np.ndarray):
        raise ValueError("adjacency_matrix 必须是numpy.array类型")
    
    if len(node_names) != adjacency_matrix.shape[0] or len(node_names) != adjacency_matrix.shape[1]:
        raise ValueError(f"节点数量({len(node_names)})与邻接矩阵维度({adjacency_matrix.shape})不匹配")
    
    # 创建图对象
    if create_directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    
    # 添加节点
    G.add_nodes_from(node_names)
    
    # 添加边（根据权重阈值过滤）
    n = len(node_names)
    edges_added = 0
    
    for i in range(n):
        for j in range(n):
            # 对于无向图，只处理上三角矩阵（避免重复边）
            if not create_directed and i >= j:
                continue
                
            weight = adjacency_matrix[i, j]
            
            # 过滤权重低于阈值的边（包括权重为0的边）
            if weight >= weight_threshold:
                G.add_edge(node_names[i], node_names[j], weight=weight)
                edges_added += 1
    
    print(f"图创建完成:")
    print(f"  - 节点数量: {G.number_of_nodes()}")
    print(f"  - 边数量: {G.number_of_edges()}")
    print(f"  - 添加的边数量: {edges_added}")
    print(f"  - 权重阈值: {weight_threshold}")
    print(f"  - 图类型: {'有向图' if create_directed else '无向图'}")
    
    return G

def set_uniform_edge_weights(G, weight_value=0.5, weight_attribute='weight'):
    """
    为图中的所有边设置统一的权重值
    
    参数:
        G: networkx.Graph 或 networkx.DiGraph，输入的图
        weight_value: float，要设置的权重值，默认为0.5
        weight_attribute: str，权重属性的名称，默认为'weight'
        
    返回:
        networkx.Graph: 更新权重后的图（修改原图，也返回）
    """
    
    if not isinstance(G, (nx.Graph, nx.DiGraph)):
        raise ValueError("输入必须是networkx.Graph或networkx.DiGraph类型")
    
    # 获取图中的所有边
    edges = list(G.edges())
    
    if not edges:
        print("警告: 图中没有边")
        return G
    
    # 为每条边设置权重
    for u, v in edges:
        # 如果图有多个平行边（MultiGraph），需要特殊处理
        if G.is_multigraph():
            # 获取节点u和v之间的所有边键
            edge_keys = G[u][v].keys()
            for key in edge_keys:
                G[u][v][key][weight_attribute] = weight_value
        else:
            G[u][v][weight_attribute] = weight_value
    
    print(f"已为 {len(edges)} 条边设置权重为 {weight_value}")
    print(f"权重属性名称: '{weight_attribute}'")
    
    return G

def get_adjacency_with_uniform_weights(G, weight_value=0.5, weight_attribute='weight'):
    """
    在不修改原图的情况下，返回添加统一权重的新图的邻接矩阵
    
    参数:
        G: networkx.Graph 或 networkx.DiGraph，输入的图
        weight_value: float，要设置的权重值，默认为0.5
        weight_attribute: str，权重属性的名称，默认为'weight'
        
    返回:
        tuple: (邻接矩阵numpy.array, 节点列表)
    """
    
    if not isinstance(G, (nx.Graph, nx.DiGraph)):
        raise ValueError("输入必须是networkx.Graph或networkx.DiGraph类型")
    
    # 获取节点列表（保持原图节点顺序）
    nodes = list(G.nodes())
    n = len(nodes)
    
    # 创建节点到索引的映射
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    
    # 初始化邻接矩阵（全零）
    adjacency_matrix = np.zeros((n, n), dtype=float)
    
    # 遍历原图的边，应用统一权重
    for u, v, data in G.edges(data=True):
        i = node_to_idx[u]
        j = node_to_idx[v]
        
        # 如果是无向图，设置对称位置
        if not G.is_directed():
            adjacency_matrix[i][j] = weight_value
            adjacency_matrix[j][i] = weight_value
        else:
            adjacency_matrix[i][j] = weight_value
    
    print(f"邻接矩阵信息:")
    print(f"  - 维度: {adjacency_matrix.shape}")
    print(f"  - 统一权重值: {weight_value}")
    print(f"  - 非零元素数量: {np.count_nonzero(adjacency_matrix)}")
    print(f"  - 矩阵密度: {np.count_nonzero(adjacency_matrix) / (n*n):.4f}")
    
    return adjacency_matrix, nodes

# def create_weighted_graph_from_adjacency_with_metrics(node_names, adjacency_matrix, weight_threshold=0.5, 
#                                                        create_directed=False, add_self_loops=False):
#     """
#     增强版：根据节点名称列表和邻接矩阵创建带权重图，并提供更多统计信息
    
#     参数:
#         node_names: list，节点名称列表
#         adjacency_matrix: numpy.array，邻接矩阵
#         weight_threshold: float，权重阈值
#         create_directed: bool，是否创建有向图
#         add_self_loops: bool，是否添加自环边（对角线元素），默认为False
        
#     返回:
#         tuple: (图对象, 统计信息字典)
#     """
    
#     # 参数验证
#     if len(node_names) != adjacency_matrix.shape[0]:
#         raise ValueError("节点数量与邻接矩阵维度不匹配")
    
#     # 创建图对象
#     G = nx.DiGraph() if create_directed else nx.Graph()
    
#     # 添加节点
#     G.add_nodes_from(node_names)
    
#     # 统计信息
#     stats = {
#         'total_possible_edges': 0,
#         'edges_added': 0,
#         'edges_filtered': 0,
#         'min_weight': float('inf'),
#         'max_weight': float('-inf'),
#         'avg_weight': 0,
#         'weight_sum': 0
#     }
    
#     n = len(node_names)
    
#     # 添加边（根据权重阈值过滤）
#     for i in range(n):
#         for j in range(n):
#             # 对于无向图，只处理上三角矩阵（避免重复边）
#             if not create_directed and i > j:
#                 continue
            
#             # 是否跳过自环
#             if not add_self_loops and i == j:
#                 continue
            
#             weight = adjacency_matrix[i, j]
            
#             # 更新权重统计
#             if weight != 0:  # 只统计非零权重
#                 stats['min_weight'] = min(stats['min_weight'], weight)
#                 stats['max_weight'] = max(stats['max_weight'], weight)
#                 stats['weight_sum'] += weight
#                 stats['total_possible_edges'] += 1
            
#             # 检查是否添加边
#             if weight >= weight_threshold:
#                 G.add_edge(node_names[i], node_names[j], weight=weight)
#                 stats['edges_added'] += 1
#             else:
#                 stats['edges_filtered'] += 1
    
#     # 计算平均权重
#     if stats['total_possible_edges'] > 0:
#         stats['avg_weight'] = stats['weight_sum'] / stats['total_possible_edges']
    
#     # 如果最小权重还是无穷大，说明没有边
#     if stats['min_weight'] == float('inf'):
#         stats['min_weight'] = 0
    
#     if stats['max_weight'] == float('-inf'):
#         stats['max_weight'] = 0
    
#     return G, stats


# # 使用示例
# if __name__ == "__main__":
#     # 示例数据
#     node_names = ["A", "B", "C", "D", "E"]
    
#     # 创建一个示例邻接矩阵（对称矩阵，表示无向图）
#     adjacency_matrix = np.array([
#         [0.0, 0.8, 0.3, 0.9, 0.1],
#         [0.8, 0.0, 0.7, 0.2, 0.6],
#         [0.3, 0.7, 0.0, 0.4, 0.8],
#         [0.9, 0.2, 0.4, 0.0, 0.5],
#         [0.1, 0.6, 0.8, 0.5, 0.0]
#     ])
    
#     # 创建一个不对称矩阵，表示有向图
#     directed_adj_matrix = np.array([
#         [0.0, 0.8, 0.0, 0.0, 0.0],
#         [0.1, 0.0, 0.7, 0.0, 0.0],
#         [0.0, 0.2, 0.0, 0.6, 0.0],
#         [0.0, 0.0, 0.3, 0.0, 0.9],
#         [0.5, 0.0, 0.0, 0.1, 0.0]
#     ])
    
#     print("示例1: 创建无向图")
#     print("-" * 40)
#     graph1 = create_weighted_graph_from_adjacency(
#         node_names=node_names,
#         adjacency_matrix=adjacency_matrix,
#         weight_threshold=0.5,
#         create_directed=False
#     )
    
#     # 输出图的边和权重
#     print("\n图的边和权重:")
#     for edge in graph1.edges(data=True):
#         print(f"  {edge[0]} -- {edge[1]}: 权重 = {edge[2]['weight']:.2f}")
    
#     print("\n" + "=" * 60)
#     print("示例2: 创建有向图")
#     print("-" * 40)
#     graph2 = create_weighted_graph_from_adjacency(
#         node_names=node_names,
#         adjacency_matrix=directed_adj_matrix,
#         weight_threshold=0.3,
#         create_directed=True
#     )
    
#     # 输出图的边和权重
#     print("\n图的边和权重:")
#     for edge in graph2.edges(data=True):
#         print(f"  {edge[0]} -> {edge[1]}: 权重 = {edge[2]['weight']:.2f}")
    
#     print("\n" + "=" * 60)
#     print("示例3: 使用增强版函数")
#     print("-" * 40)
#     graph3, stats = create_weighted_graph_from_adjacency_with_metrics(
#         node_names=node_names,
#         adjacency_matrix=adjacency_matrix,
#         weight_threshold=0.5,
#         create_directed=False,
#         add_self_loops=False
#     )
    
#     print(f"\n统计信息:")
#     print(f"  总可能边数: {stats['total_possible_edges']}")
#     print(f"  添加的边数: {stats['edges_added']}")
#     print(f"  过滤的边数: {stats['edges_filtered']}")
#     print(f"  最小权重: {stats['min_weight']:.2f}")
#     print(f"  最大权重: {stats['max_weight']:.2f}")
#     print(f"  平均权重: {stats['avg_weight']:.2f}")
    
#     # 图的基本分析
#     print(f"\n图分析:")
#     print(f"  节点数量: {graph3.number_of_nodes()}")
#     print(f"  边数量: {graph3.number_of_edges()}")
#     print(f"  是否连通: {nx.is_connected(graph3) if not graph3.is_directed() else 'N/A（有向图）'}")
