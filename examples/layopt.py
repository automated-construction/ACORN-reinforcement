import itertools as it
from acorn_reinforcement.strutandtie import StrutAndTie
from acorn_reinforcement.optimisation import lo, lo_plot

# input data: geometry, connectivity and statics system
st, sc = 1, 1
nodes = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
edges = [(i, j, st, sc, 1) for i, j in it.combinations(range(len(nodes)), 2)] # the tensile and compressive strength st and sc can be custom for each edge
supports = {0: [0, 0, 0], 1: [1, 0, 0]}
loads = [{2: [1.0, 0.0, 0.0]}, {1: [1.0, 0.0, 0.0]}]

# assemble SAT object
sat = StrutAndTie.assemble(nodes, edges, supports, loads)

# disassemble data
nodes, edges, supports, loads = sat.disassemble()

# optimise layout
Nd, Cn, a, q = lo(nodes, edges, supports, loads)

# plot layout
lo_plot(Nd, Cn, a, q, max(a) * 1e-3, str='Optimal layout', update=False)

# update force data
for i in range(len(Cn)):
	sat.edge_force((Cn[i][0], Cn[i][1]), {k: q[k][i] for k in range(len(q))})

# check equilibrium
print('is in equilibrium?', sat.is_in_equilibrium())

# output data
output = [(u, v, sat.edge_force((u, v))) for u, v in sat.edges()]