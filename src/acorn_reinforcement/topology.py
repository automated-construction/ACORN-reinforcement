from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import networkx as nx


__all__ = []


def connected_components_nodes(G):
	return list(nx.connected_components(G))


def connected_components_edges(G):
	edges = list(G.edges(keys=True))

	components = []
	for nodes in connected_components_nodes(G):
		component = []
		for edge in reversed(edges):
			u, v, c = edge
			if u in nodes:
				component.append(edge)
				edges.remove(edge)
		components.append(component)

	return components


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

	pass