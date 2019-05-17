
import networkx as nx
from networkx.readwrite import json_graph


def save_trees(trees, filename, format="json"):
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
                    if "scanids" in tc.node[n]:
                        tc.node[n]["scanids"] = str(tc.node[n]["scanids"])
                    if "mf" in tc.node[n]:
                        tc.node[n]["mf"] = str(tc.node[n]["mf"])
                for line in nx.readwrite.generate_gml(tc):
                    out.write((line + "\n").encode('ascii'))
            else:
                raise ValueError("Incorrect format - json or gml")


def load_trees(filename, format="json"):

    def sort_graph(G):
        G_sort = nx.OrderedDiGraph()
        G_sort.graph["id"] = G.graph["id"]
        G_sort.add_nodes_from(sorted(G.nodes(data=True), key=lambda x: x[1]['order']))
        G_sort.add_edges_from(sorted(G.edges(data=True), key=lambda x: x[2]['order']))
        return G_sort

    def remove_attr(G, atr):
        for n in G.nodes():
            del G.node[n][atr]
        for e in G.edges():
            del G[e[0]][e[1]][atr]
        return G

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
                    if "scanids" in G.node[n]:
                        G.node[n]["scanids"] = eval(G.node[n]["scanids"])
                    if "mf" in G.node[n]:
                        G.node[n]["mf"] = eval(G.node[n]["mf"])
                graphs.append(remove_attr(sort_graph(G), "order"))
            return graphs
        else:
            raise ValueError("Incorrect graph format - json or gml")


def save_groups(groups, filename, format="json"):
    save_trees(trees=groups, filename=filename, format=format)
    return


def load_groups(filename, format="json"):
    return load_trees(filename=filename, format=format)


