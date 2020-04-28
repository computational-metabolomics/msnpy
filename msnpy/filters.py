#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019-2020 Ralf Weber
#
# This file is part of MSnPy.
#
# MSnPy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MSnPy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MSnPy.  If not, see <https://www.gnu.org/licenses/>.
#


import networkx as nx
from networkx import isolates, get_node_attributes


def validate_injection_time_ms1(G: nx.OrderedDiGraph, max_injection_time: float, label: str = " ms "):

    """

    :param G:
    :param max_injection_time:
    :param label:
    :return:
    :rtype: nx.OrderedDiGraph
    """

    h = list(nx.topological_sort(G))[0]
    if len(G.nodes[h]["ioninjectiontimes"]) > 1:
        raise ValueError("Validation only valid for trees no replication for MS1")
    if label in h:
        # should only have a single scanid
        if G.nodes[h]["ioninjectiontimes"][0] > max_injection_time:
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
        G.nodes[n[0]]["flag"] = False
    # isolated edges / dependent scan events
    for e in G.edges(data=True):
        if len(G.nodes[e[0]]["scanids"]) == 0 and len(G.nodes[e[1]]["scanids"]) > 0:
            G.nodes[e[0]]["flag"] = False
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
            G.nodes[n[0]]["flag"] = False

    for e in G.edges():
        if len(G.nodes[e[0]]["scanids"]) < min_replicates[G.nodes[e[0]]["mslevel"]-1]:
            G.nodes[e[0]]["flag"] = False
            G.nodes[e[1]]["flag"] = False
    return G
