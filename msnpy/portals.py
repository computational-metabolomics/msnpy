
from typing import Sequence
import networkx as nx
from networkx.readwrite import json_graph


def save_trees(trees: Sequence[nx.classes.ordered.OrderedDiGraph], filename: str, format: str = "json"):

    """

    :param trees:
    :param filename:
    :param format:
    :return:
    """

    with open(filename, "w") as out:
        for t in trees:
            tc = t.copy()  # required to not update Graph
            for i, e in enumerate(tc.edges()):
                tc[e[0]][e[1]]["order"] = i
            for i, n in enumerate(tc.nodes()):
                tc.node[n]["order"] = i
            if format == "json":
                out.write(str(nx.readwrite.json_graph.node_link_data(tc)) + "\n")
            elif format == "gml":
                for i, n in enumerate(tc.edges()):
                    if "mf" in tc[n[0]][n[1]]:
                        tc[n[0]][n[1]]["mf"] = str(tc[n[0]][n[1]]["mf"])
                for i, n in enumerate(tc.nodes()):
                    for k in ["scanids", "ioninjectiontimes", "mf", "coltype", "template"]:
                        if k in tc.node[n]:
                            tc.node[n][k] = str(tc.node[n][k])
                for line in nx.readwrite.generate_gml(tc):
                    out.write((line + "\n"))

            else:
                raise ValueError("Incorrect format - json or gml")


def load_trees(filename: str, format: str = "json"):

    """

    :param filename:
    :param format:
    :return:
    """

    def sort_graph(G_):
        G_sort = nx.OrderedDiGraph()
        G_sort.graph["id"] = G_.graph["id"]
        G_sort.add_nodes_from(sorted(G_.nodes(data=True), key=lambda x: x[1]['order']))
        G_sort.add_edges_from(sorted(G_.edges(data=True), key=lambda x: x[2]['order']))
        return G_sort

    def remove_attr(G_, atr):
        for n_ in G_.nodes():
            del G_.node[n_][atr]
        for e in G_.edges():
            del G_[e[0]][e[1]][atr]
        return G_

    with open(filename, "r") as inp:
        graphs = []
        if format == "json":
            for line in inp.readlines():
                G = json_graph.node_link_graph(eval(line))
                graphs.append(remove_attr(sort_graph(G), "order"))
            return graphs
        elif format == "gml":
            for gml_str in inp.read().split("graph")[1:]:
                G = nx.readwrite.parse_gml("graph" + gml_str)
                for n in G.nodes():
                    if "coltype" in G.node[n]:
                        if G.node[n]["coltype"] == "None":
                            G.node[n]["coltype"] = None
                    for k in ["scanids", "ioninjectiontimes", "mf", "template"]:
                        if k in G.node[n]:
                            G.node[n][k] = eval(G.node[n][k])
                graphs.append(remove_attr(sort_graph(G), "order"))
            return graphs
        else:
            raise ValueError("Incorrect graph format - json or gml")


def save_groups(groups: Sequence[nx.classes.ordered.OrderedDiGraph], filename: str, format: str = "json"):

    """

    :param groups:
    :param filename:
    :param format:
    :return:
    """

    save_trees(trees=groups, filename=filename, format=format)
    return


def load_groups(filename: str, format: str = "json"):

    """

    :param filename:
    :param format:
    :return:
    :rtype: list of NetworkX Graphs
    """

    return load_trees(filename=filename, format=format)

