from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from operator import itemgetter

from itertools import combinations, permutations

import networkx as nx

from compas.geometry import distance_point_point

from compas.utilities import pairwise

from acorn_reinforcement.misc import remove_pair_in_dict


__all__ = [
	'rpp',
]

def rpp(required_edges, optional_edges, key_to_xyz):
	key_to_xyz = {int(key): [float(xyz) for xyz in value] for key, value in key_to_xyz.items()}

	def edge_length(u, v, key_to_xyz):
		return distance_point_point(key_to_xyz[u], key_to_xyz[v])
	def path_length(path, key_to_xyz):
		return sum([edge_length(u, v, key_to_xyz) for u, v in pairwise(path)])
	
	woundable_graph = eulerize_path_general(required_edges, optional_edges, path_cost_func=path_length, path_cost_func_args=key_to_xyz)
	return list(nx.eulerian_path(woundable_graph))


def eulerize_circuit(G):
	return nx.eulerize(G)


def eulerize_path(G, Gtot, path_cost_func, path_cost_func_args):

	odd_degree_nodes = [n for n, d in G.degree() if d % 2 == 1]
	odd_deg_pairs_paths = [(m, {n: nx.shortest_path(Gtot, source=m, target=n)}) for m, n in combinations(odd_degree_nodes, 2)]

	min_path_solution = None
	min_cost = None
	for extremities in combinations(odd_degree_nodes, 2):
		# apply eulerize_circuit using a zero-cost edge between extremities
		Gp = nx.Graph()
		path_cost = {}
		for n, Ps in odd_deg_pairs_paths:
			if n not in extremities:
				for m, P in Ps.items():
					if m not in extremities:
						if n != m:
							cost = path_cost_func(P, path_cost_func_args)
							path_cost[(P[0], P[-1])] = cost	
							weight = 1 / cost
							# distortion due to 1/x instead of -x?
							Gp.add_edge(m, n, weight=weight, path=P)

		best_matching = nx.Graph(list(nx.max_weight_matching(Gp)))

		total_cost = 0
		for u, v in best_matching.edges():
			if (u, v) in path_cost:
				total_cost += path_cost[(u, v)]
			else:
				total_cost += path_cost[(v, u)]

		if min_cost is None or total_cost < min_cost:
			min_cost = total_cost
			min_path_solution = (Gp, best_matching)

	Gp, best_matching = min_path_solution

	print('edges added within components:', [Gp[m][n]["path"] for m, n in best_matching.edges()])

	for m, n in best_matching.edges():
		path = Gp[m][n]["path"]
		G.add_edges_from(pairwise(path))
	
	return G


def eulerize_path_general(required_edges, optional_edges, path_cost_func, path_cost_func_args):
	# greedy approach: each subpath is optimal but the general path may not be optimal

	G = nx.MultiGraph(required_edges)
	Gtot = nx.MultiGraph(required_edges + optional_edges)

	# get disconnected subgraphs
	connected_edges = connected_components_edges(G)
	connected_graphs = [nx.MultiGraph(edges) for edges in connected_edges]

	# eulerise each subgraph
	eulerian_subpaths = []
	for graph in connected_graphs:
		if nx.is_eulerian(graph):
			path = euler_paths_from_eulerian(graph)
		elif nx.is_semieulerian(graph):
			path = euler_paths_from_semieulerian(graph)
		else:
			eulerize_path(graph, Gtot, path_cost_func, path_cost_func_args)
			path = euler_paths_from_semieulerian(graph)
		eulerian_subpaths.append(path)
	G = nx.MultiGraph([edge for graph in connected_graphs for edge in graph.edges()])

	# eulerise all subgraphs
	# path format is list of edges (u, v, c)
	paths_extremities = [(path[0][0], path[-1][-2]) for path in eulerian_subpaths] 

	# change G for a simplified graph with one edge for each sub-graph
	eulerize_disconnected_paths(G, Gtot, paths_extremities, path_cost_func, path_cost_func_args)

	return G


def euler_paths_from_eulerian(G):
	circuit = list(nx.eulerian_circuit(G, keys=True))
	return circuit


def euler_paths_from_semieulerian(G):
	path = list(nx.eulerian_path(G, keys=True))
	return path	

def eulerize_disconnected_paths(G, Gtot, paths_extremities, path_cost_func, path_cost_func_args):

	def find_shortest_path_from_targets(source, targets, path_cost_func, path_cost_func_args):
		min_path = None
		min_cost = None
		for target in targets:
			path = nx.shortest_path(Gtot, source, target)
			cost = path_cost_func(path, path_cost_func_args)
			if min_cost is None or cost < min_cost:
				min_path = path
				min_cost = cost
		return min_path, min_cost


	added_paths = []
	# pair path extremities
	nodes = [node for extremities in paths_extremities for node in extremities]
	node_pairs = {}
	for u, v in paths_extremities:
		node_pairs[u] = v
		node_pairs[v] = u

	b1, b2 = paths_extremities[0]
	nodes.remove(b1)
	nodes.remove(b2)

	# iteratively select and remove shortest paths
	for _ in range(len(paths_extremities) - 1):
		path_a, cost_a = find_shortest_path_from_targets(b1, nodes, path_cost_func, path_cost_func_args)
		path_c, cost_c = find_shortest_path_from_targets(b2, nodes, path_cost_func, path_cost_func_args)

		if cost_a < cost_c:
			added_paths.append(path_a)
			G.add_edges_from(nx.utils.pairwise(path_a))
			a2 = path_a[-1]
			a1 = node_pairs[a2]

			# a1-(a2---b1)-b2
			node_pairs[a1] = b2
			node_pairs[b2] = a1
			if b1 != b2:
				del node_pairs[b1]
			if a1 != a2:
				del node_pairs[a2]		

			b1 = a1
			nodes.remove(a1)
			nodes.remove(a2)
		
		else:
			added_paths.append(path_c)
			G.add_edges_from(nx.utils.pairwise(path_c))

			c1 = path_c[-1]
			c2 = node_pairs[c1]

			# b1-(b2---c1)-c2
			node_pairs[b1] = c2
			node_pairs[c2] = b1
			if b1 != b2:
				del node_pairs[b2]
			if c1 != c2:
				del node_pairs[c1]

			b2 = c2
			nodes.remove(c1)
			nodes.remove(c2)
	
	print('edges added between components:', added_paths)

	return G

def eulerize_disconnected_paths_2(G, Gtot, paths_extremities, path_cost_func, path_cost_func_args):

	added_paths = []
	print('paths_extremities', len(paths_extremities))
	# pair path extremities
	nodes = [node for extremities in paths_extremities for node in extremities]
	node_pairs = {}
	for u, v in paths_extremities:
		node_pairs[u] = v
		node_pairs[v] = u
	print('nodes', len(nodes))
	print('nodes set', len(list(set(nodes))))
	# compute all shortest paths and their costs
	shortest_paths_pairs = {}
	for u, v in combinations(nodes, 2):
		if v != node_pairs[u]:
			path = nx.shortest_path(Gtot, source=u, target=v)
			cost = path_cost_func(path, path_cost_func_args)
			shortest_paths_pairs[(u, v)] = (path, cost)
		else:
			print('exception', v, node_pairs[u])
	print('shortest_paths_pairs', len(shortest_paths_pairs))
	# iteratively select and remove shortest paths
	for _ in range(len(paths_extremities) - 1):
		print(len(shortest_paths_pairs))
		# select the min-cost shortest path
		((a2, b1), (path, cost)) = min(shortest_paths_pairs.items(), key=lambda tup: itemgetter(1)(tup[1]))

		added_paths.append(path)
		G.add_edges_from(nx.utils.pairwise(path))

		# a1-a2---b1-b2
		a1 = node_pairs[a2]
		b2 = node_pairs[b1]

		to_delete = []
		for pair in shortest_paths_pairs:
			if a2 in pair or b1 in pair or (a1 in pair and b2 in pair):
				to_delete.append(pair)
		for pair in to_delete:
			remove_pair_in_dict(shortest_paths_pairs, pair)

		# a1-(a2---b1)-b2
		node_pairs[a1] = b2
		node_pairs[b2] = a1
		del node_pairs[a2]
		del node_pairs[b1]
	
	print('edges added between components:', added_paths)

	return G


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

	import itertools as it
	import networkx as nx

	key_to_xyz = {
		0: [0.0, 0.0, 0.0],
		1: [1.0, 0.0, 0.0],
		2: [2.0, 0.0, 0.0],
		3: [0.0, 1.0, 0.0],
		4: [1.0, 1.0, 0.0],
		5: [2.0, 1.0, 0.0],
		6: [0.0, 2.0, 0.0],
		7: [1.0, 2.0, 0.0],
		8: [2.0, 2.0, 0.0],
		9: [0.0, 3.0, 0.0],
		10: [1.0, 3.0, 0.0],
		11: [2.0, 3.0, 0.0],
	}

	required_edges_1 = [(0,3), (1,4), (2,5)]
	required_edges_2 = [(0,3), (1,4), (2,5), (6,3), (7,4), (8,5)]
	required_edges_3 = [(0,1), (1,4), (4,3), (3,0), (5,8)]
	required_edges_4 = [(0,1), (1,3), (3,0), (5,8), (8,7), (7,5)]
	required_edges_4 = [(0,1), (1,3), (3,0), (5,8), (8,7), (7,5), (6,9)]

	optional_edges = list(it.combinations(key_to_xyz.keys(), 2))
	
	path = rpp(required_edges_4, optional_edges, key_to_xyz)

	print('path:', path)
