from networkx import isolates, get_node_attributes
import networkx as nx


def validate_injection_time_ms1(G: nx.OrderedDiGraph, max_injection_time: float, label: str = " ms "):

    """

    :param G:
    :param max_injection_time:
    :param label:
    :return:
    :rtype: nx.OrderedDiGraph
    """

    h = list(nx.topological_sort(G))[0]
    if len(G.node[h]["ioninjectiontimes"]) > 1:
        raise ValueError("Validation only valid for trees no replication for MS1")
    if label in h:
        # should only have a single scanid
        if G.node[h]["ioninjectiontimes"][0] > max_injection_time:
            return False
        else:
            return True
    else:
        raise ValueError("Root node (parent) in tree is not a full scan")


def filter_by_isolation(G: nx.OrderedDiGraph):

    """
    Flag 'isolated' scan events.

    :param G:
    :return:
    :rtype: nx.OrderedDiGraph
    """

    for n in isolates(G):
        G.node[n[0]]["flag"] = False
    # isolated edges / dependent scan events
    for e in G.edges(data=True):
        if len(G.node[e[0]]["scanids"]) == 0 and len(G.node[e[1]]["scanids"]) > 0:
            G.node[e[0]]["flag"] = False
    return G


def filter_by_replicates(G: nx.OrderedDiGraph, min_replicates: int):

    """

    :param G:
    :param min_replicates:
    :return:
    :rtype: nx.OrderedDiGraph
    """

    def max_mslevel(G_: nx.OrderedDiGraph):
        return max(get_node_attributes(G_, "mslevel").values())

    if isinstance(min_replicates, int):
        min_replicates = max_mslevel(G) * [min_replicates]
    elif not isinstance(min_replicates, list) or len(min_replicates) < max_mslevel(G):
        raise ValueError("Provide a singler integer or a list of integers to specify the minimum number of replciates for each ms level.")

    for n in G.nodes(data=True):
        if len(n[1]["scanids"]) < min_replicates[n[1]["mslevel"]-1]:
            G.node[n[0]]["flag"] = False

    for e in G.edges():
        if len(G.node[e[0]]["scanids"]) < min_replicates[G.node[e[0]]["mslevel"]-1]:
            G.node[e[0]]["flag"] = False
            G.node[e[1]]["flag"] = False
    return G
