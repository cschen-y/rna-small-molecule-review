import networkx as nx

def get_node_topology(adj_matrix):
    """
    计算给定邻接矩阵的节点度和接近中心性。

    参数：
        adj_matrix (numpy.ndarray): 图的邻接矩阵。

    返回：
        list: 每个节点的属性列表 [度, 接近中心性]。
    """
    
    G = nx.from_numpy_matrix(adj_matrix)

    
    DG = dict(G.degree())

    
    CL = nx.closeness_centrality(G)  

    
    all_properties = []
    for node in G.nodes():
        node_properties = [DG[node], CL[node]]
        all_properties.append(node_properties)

    return all_properties
