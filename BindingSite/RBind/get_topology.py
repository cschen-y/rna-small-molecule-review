import numpy as np
import networkx as nx

def get_node_topology(adj_matrix):
    """
    计算给定邻接矩阵的图的节点拓扑属性。

    参数：
        adj_matrix (numpy.ndarray): 图的邻接矩阵。

    返回：
        list: 每个节点的属性列表 [度, 接近中心性, 邻居平均度, 中介中心性, 离心率]。
    """
    
    G = nx.from_numpy_array(adj_matrix)

    
    DG = dict(G.degree())

    
    NC = {
        node: np.mean([DG[neighbor] for neighbor in G.neighbors(node)]) if DG[node] > 0 else 0
        for node in G.nodes()
    }

    
    BC = nx.betweenness_centrality(G)  
    CL = nx.closeness_centrality(G)  

    
    try:
        EC = nx.eccentricity(G)  
    except nx.NetworkXError:
        
        EC = {}
        for component in nx.connected_components(G):
            subgraph = G.subgraph(component)
            sub_eccentricity = nx.eccentricity(subgraph)
            EC.update(sub_eccentricity)

    
    all_properties = []
    for node in G.nodes():
        node_properties = [DG[node], CL[node], NC[node], BC[node], EC[node]]
        all_properties.append(node_properties)

    return all_properties