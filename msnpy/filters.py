
from networkx import isolates, get_node_attributes

def filter_ms1_by_injection_time(G, max_injection_time, label=" ms "):

    h = list(G.edges())[0][0]
    if label in h:
        if G.node[h]["injectiontime"] > max_injection_time:
            return False
        else:
            return True
    else:
        raise ValueError("Root node in tree is not a full scan")


def filter_by_replicates(G, min_replicates):

    def max_mslevel(G):
        return max(get_node_attributes(G, "mslevel").values())

    if isinstance(min_replicates, int):
        min_replicates = max_mslevel(G) * [min_replicates]
    elif not isinstance(min_replicates, list) or len(min_replicates) < max_mslevel(G):
        raise ValueError("Provide a singler integer or a list of integers to specify the minimum number of replciates for each ms level.")

    for n in G.nodes(data=True):
        if len(n[1]["scanids"]) < min_replicates[n[1]["mslevel"]-1]:
            G.node[n[0]]["scanids"] = []

    for e in G.edges():
        if len(G.node[e[0]]["scanids"]) < min_replicates[G.node[e[0]]["mslevel"]-1]:
            G.node[e[0]]["scanids"] = []
            G.node[e[1]]["scanids"] = []

    for n in isolates(G):
        G.node[n[0]]["scanids"] = []

    return G
