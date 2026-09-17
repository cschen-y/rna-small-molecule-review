from RGCN.utils import available_pdbids
from RGCN.utils import graph_from_pdbid

pdbids = available_pdbids()
graphs = [graph_from_pdbid(p) for p in pdbids]
for i,G in enumerate(graphs):
    try:
        binding_nodes = [(n, d['binding_ion']) for n, d in G.nodes(data=True) if d['binding_ion'] is not None]
    except KeyError:
        print(pdbids[i])

















from RGCN.prepare_data.main import cif_to_graph























