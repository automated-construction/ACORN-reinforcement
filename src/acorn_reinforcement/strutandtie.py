from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import ceil

from compas.datastructures import Graph

from compas.geometry import scale_vector, normalize_vector, subtract_vectors, sum_vectors

from compas.utilities import geometric_key


__all__ = ['StrutAndTie']


class StrutAndTie(Graph):

	def __init__(self):
		super(StrutAndTie, self).__init__()
		self.update_default_edge_attributes({'force': {}, 'strength_tension': 1, 'strength_compression': 1, 'start_with': True, 'cross_section': None})
		self.update_default_node_attributes({'support': [1, 1, 1], 'load': {}})

	@classmethod
	def assemble(cls, nodes, edges, supports, loads):
		sat = cls()

		for i, xyz in enumerate(nodes):
			attr = {axis: value for axis, value in zip('xyz', xyz)}
			attr.update({'support': supports.get(i, [1, 1, 1])})
			attr.update({'load': {k: loads[k].get(i, [0.0, 0.0, 0.0]) for k in range(len(loads))}})
			sat.add_node(i, attr)
		
		for u, v, inc, st, sc in edges:
			attr = {'force': {k: 0.0 for k in range(len(loads))}, 'start_with': inc, 'strength_tension': st, 'strength_compression': sc}
			sat.add_edge(u, v, attr)

		return sat

	def disassemble(self):
		key_to_idx = {key: i for i, key in enumerate(self.nodes())}
		nodes = [self.node_coordinates(key) for key in self.nodes()]
		edges = [(key_to_idx[u], key_to_idx[v], self.edge_attribute((u, v), 'start_with'), self.edge_attribute((u, v), 'strength_tension'), self.edge_attribute((u, v), 'strength_compression')) for u, v in self.edges()]
		supports = {key: self.node_support(key) for key in self.nodes()}
		loads = [{key: self.node_load(key)[k] for key in self.nodes()} for k in range(len(self.node_load(list(self.nodes())[0])))]
		return nodes, edges, supports, loads

# ==============================================================================
# Statics system
# ==============================================================================

	def node_support(self, key, support=None):
		# support format is a list of 3 binaries, 0 being fixed and 1 being free
		if support:
			return self.node_attribute(key, 'support', support)
		else:
			return self.node_attribute(key, 'support')

	def node_load(self, key, load=None):
		# load format is a list of 3 floats
		if load:
			return self.node_attribute(key, 'load', load)
		else:
			return self.node_attribute(key, 'load')

# ==============================================================================
# Forces
# ==============================================================================

	def edge_force(self, edge, force=None):
		# forces in kN, compression < 0 and tension > 0
		if force:
			return self.edge_attribute(edge, 'force', force)
		else:
			return self.edge_attribute(edge, 'force')

	def is_in_equilibrium(self, tol=1e-6):
		for node in self.nodes():
			forces = [self.node_load(node)]
			for nbr in self.neighbors(node):
				u, v, eps = (node, nbr, -1) if nbr in self.edge[node] else (nbr, node, 1)
				force = self.edge_force((u, v))
				force = {k: scale_vector(normalize_vector(subtract_vectors(self.node_coordinates(u), self.node_coordinates(v))), eps * f) for k, f in force.items()}
				forces.append(force)
			for k in range(len(forces[0])):
				forces_k = [force[k] for force in forces]
				fr = sum_vectors(forces_k)
				supp = self.node_support(node)
				for i in range(3):
					if abs(supp[i] * fr[i]) > tol:
						return False
		return True

	def store_edge_forces(self, edge_to_force):
		# forces in kN, compression > 0 and tension < 0
		for edge, force in edge_to_force.items():
			self.edge_force(edge, force)

	def set_line_forces(self, line_to_force, precision=None):
		# forces in kN, compression > 0 and tension < 0
		vertex_map = {geometric_key(key, precision): key for key in self.nodes()}
		edge_to_force = {}
		for (u_xyz, v_xyz), force in line_to_force.items():
			u_key = geometric_key(u_xyz, precision)
			v_key = geometric_key(v_xyz, precision)
			u = vertex_map[u_key]
			v = vertex_map[v_key]
			edge_to_force[(u, v)] = force
		self.store_edge_forces(edge_to_force)
		return edge_to_force

	# def struts(self):
	# 	return [edge for edge in self.edges() if self.edge_force(edge) < 0]

	# def ties(self):
	# 	return [edge for edge in self.edges() if self.edge_force(edge) > 0]

	# def nulls(self):
	# 	return [edge for edge in self.edges() if self.edge_force(edge) == 0]

# ==============================================================================
# Load path
# ==============================================================================

	def load_path(self, compression=True, tension=True):
		lp = 0
		if compression:
			lp += sum([abs(self.edge_force(edge)) * self.edge_length(*edge) for edge in self.edges() if self.edge_force(edge) > 0])
		if tension:
			lp += sum([abs(self.edge_force(edge)) * self.edge_length(*edge) for edge in self.edges() if self.edge_force(edge) < 0])
		return lp

# ==============================================================================
# Cross-sections
# ==============================================================================

	def edge_cross_section(self, edge, cross_section=None):
		if cross_section is None:
			return self.edge_attribute(edge, 'cross_section')
		else:
			return self.edge_attribute(edge, 'cross_section', cross_section)

	def store_cross_sections(self, axial_strength, edges=None):
		# cross_section in m2 and axial_strength in MPa
		# if no edges are specified, apply to all edges
		if edges is None:
			edges = self.edges()
		for edge in edges:
			cross_section = abs(self.edge_force(edge)) / (axial_strength * 1000)
			self.edge_cross_section(edge, cross_section)

	def compute_number_elements(self, element_cross_section, edges=None):
		# cross_section in m2
		# if no edges are specified, apply to all edges
		if edges is None:
			edges = self.edges()
		return [int(ceil(self.edge_cross_section(edge) / element_cross_section)) for edge in edges]
