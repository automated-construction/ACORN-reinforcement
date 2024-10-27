from math import gcd, ceil, sqrt
import itertools as it

import json

from compas.utilities import linspace, pairwise
from compas.geometry import add_vectors, distance_point_point, distance_point_point_xy
from compas.datastructures import Mesh

from acorn_reinforcement.strutandtie import StrutAndTie
from acorn_reinforcement.eulerisation import rpp
from acorn_reinforcement.optimisation import lo, lo_plot
from acorn_reinforcement.winding import wind_polyline

### GEOMETRY AND TOPOLOGY INPUT ###

x0, y0, z0 = 0.0, 0.0, 0.0 # min corner
dx, dy, dz = 5, 5, 0.2 # length, width, thickness
nx, ny = 6, 6 # number of nodes along the X and Y directions
h = 0.0 # parabolic shape parameter
max_length = 3.0 # max length of lines in ground structure
init_max_length = 1.5 # max length of lines in initial ground structure

### NODES ###
N = nx * ny
prod = it.product(list(linspace(x0, x0 + dx, nx)), list(linspace(y0, y0 + dy, ny)))
bottom =  [[x, y, z0] for x, y in prod]
top = [add_vectors(point, [0.0, 0.0, dz]) for point in bottom]
nodes = bottom + top
d = sqrt(dx ** 2 + dy ** 2) / 2
for node in nodes:
	x = distance_point_point_xy(node[:2], [x0 + dx / 2, y0 + dy / 2])
	node[2] += (d + x) * (d - x) * h / d ** 2

### HULL INTRADOS & EXTRADOS ###
faces_bottom, faces_top = [], []
for j in range(ny - 1):
	for i in range(nx - 1):
		k = i + nx * j
		face = [k, k + 1, k + 1 + nx, k + nx]
		faces_bottom.append(list(reversed(face)))
		faces_top.append([k + N for k in face])
mesh = Mesh.from_vertices_and_faces(nodes, faces_bottom + faces_top)
mesh.quads_to_triangles()

### EDGES ###
edges = []
for i, j in it.combinations(range(2 * N), 2):
	l = distance_point_point(nodes[i], nodes[j])
	if l <= max_length:
		start_with = 1 if l <= init_max_length else 0
		tensile_strength = 1.0 if (i < N and j < N) or (i >= N and j >= N) or abs(i - j) == N else 0.0
		compressive_strength = 1.0 
		edges.append((i, j, start_with, tensile_strength, compressive_strength))

print('{} nodes and {} edges'.format(len(nodes), len(edges)))

### SUPPORT CONDITIONS ###
supports = {0: [1, 1, 0], nx - 1: [1, 1, 0], nx * (ny - 1): [1, 1, 0], N - 1: [1, 1, 0]}

### LOAD CONDITIONS ###
loads = [
	{i: [0.0, 0.0, -0.1] for i in range(int(len(nodes) / 2), len(nodes))},
	]

### ASSEMBLE STRUT-AND-TIE MODEL ###
sat = StrutAndTie.assemble(nodes, edges, supports, loads)

### OPTIMISE STRUT-AND-TIE MODEL ###
nodes, edges, supports, loads = sat.disassemble()
Nd, Cn, a, q = lo(nodes, edges, supports, loads)
for i in range(len(Cn)):
	sat.edge_force((Cn[i][0], Cn[i][1]), {k: q[k][i] for k in range(len(q))})

threshold_f = max([qij for qi in q for qij in qi]) * 1e-3
threshold_a = max(a) * 1e-3

# plot layout
lo_plot(Nd, Cn, a, q, threshold_a, str='Optimal layout', update=False)

lines_and_forces = [(sat.node_coordinates(u), sat.node_coordinates(v), list(sat.edge_force((u, v)).values())) for u, v in sat.edges()]

load_paths = {i: 0.0 for i in sat.edge_force(list(sat.edges())[0])}
tensile_load_paths = {i: 0.0 for i in sat.edge_force(list(sat.edges())[0])}
for u, v in sat.edges():
	l = sat.edge_length((u, v))
	for i, fi in sat.edge_force((u, v)).items():
		load_paths[i] += abs(fi) * l
		tensile_load_paths[i] += max(0, fi) * l
print('load path, tensile load path and ratio per load case:', {i: (lp, tensile_load_paths[i], tensile_load_paths[i] / lp) for i, lp in load_paths.items()})

### WINDING ###

filament_tensile_strength = 1
filament_cross_section = 1

edge_req_passes = {(u, v): ceil(max(sat.edge_force((u, v)).values()) / (filament_tensile_strength * filament_cross_section) - threshold_f) for u, v in sat.edges()}
sorted_edge_req_passes = sorted(edge_req_passes.values())
min_req_passes, max_req_passes = sorted_edge_req_passes[0], sorted_edge_req_passes[-1]

# wind along intrados and extrados
key_to_xyz = {node: sat.node_coordinates(node) for node in sat.nodes()}
req_bot, opt_bot, req_top, opt_top = [], [], [], []
for (u, v), passes in edge_req_passes.items():
	if u < N and v < N:
		req, opt = req_bot, opt_bot
	elif u >= N and v >= N:
		req, opt = req_top, opt_top
	else:
		continue
	if passes > 0:
		req += [(u, v)] * passes
	else:
		opt.append((u, v))

required_filament_length = sum([sat.edge_length((u, v)) for u, v in req_bot + req_top])

print('min and max required number of passes before semi-eulerisation', min_req_passes, max_req_passes)
print('required filament length before semi-eulerisation', required_filament_length)
print('required filament volume before semi-eulerisation', required_filament_length * filament_cross_section)

path_topo = []
for req, opt in [(req_bot, opt_bot), (req_top, opt_top)]:
	edges = rpp(req, opt, key_to_xyz) if len(req) != 0 else []
	if len(edges) != 0:
		path_topo.append([sat.node_coordinates(edges[0][0])] + [sat.node_coordinates(edge[1]) for edge in edges])
	else:
		path_topo.append([])

filament_length = sum([distance_point_point(u, v) for path in path_topo for u, v in pairwise(path)])

print('min and max required number of passes after semi-eulerisation')
print('filament length after semi-eulerisation', filament_length)
print('filament volume after semi-eulerisation', filament_length * filament_cross_section)

### WINDING PATH ###
offset = dx / nx / 4
path_geom = [wind_polyline(path, offset) for path in path_topo]
winding_length = sum([distance_point_point(u, v) for path in path_geom for u, v in path])
print('winding path length', winding_length)

### EXPORT DATA ###
# data = [lines_and_forces, path_topo, path_geom]
# repository = '/Users/robinoval/Desktop/playground/'
# filename = 'reinforcement'
# with open(repository + filename + '.json', 'w') as path:
#     json.dump(data, path)
