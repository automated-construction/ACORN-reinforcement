from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import asin

from compas.geometry import distance_point_point
from compas.geometry import midpoint_point_point
from compas.geometry import add_vectors
from compas.geometry import subtract_vectors
from compas.geometry import scale_vector
from compas.geometry import dot_vectors
from compas.geometry import cross_vectors
from compas.geometry import normalize_vector
from compas.geometry import rotate_points

from compas.utilities import pairwise


__all__ = []


def wind_polyline(polyline, offset):

	lines = {i: (start, end) for i, (start, end) in enumerate(pairwise(polyline))}
	lines_vectors = {i: subtract_vectors(end, start) for i, (start, end) in lines.items()}

	# get changes of orientation at nodes
	nodes_orientations = {}
	for i in range(len(polyline)):
		if i == 0 or i == len(polyline) - 1:
			nodes_orientations[i] = None
			continue
		u = normalize_vector(lines_vectors[i - 1])
		v = normalize_vector(lines_vectors[i])
		cross = cross_vectors(u, v)
		dot = dot_vectors(cross, [0.0, 0.0, 1.0])
		nodes_orientations[i] = dot

	# offset lines
	line_circle_points = {}
	for i, (a, b) in lines.items():
		ab = scale_vector(normalize_vector(lines_vectors[i]), offset)
		ab_r = cross_vectors(ab, [0.0, 0.0, 1.0])
		a1 = add_vectors(a, ab_r)
		a2 = add_vectors(a, scale_vector(ab_r, -1))
		b1 = add_vectors(b, ab_r)
		b2 = add_vectors(b, scale_vector(ab_r, -1))
		line_circle_points[i] = (a1, a2, b1, b2)

	path = []
	cross = []
	side = 1
	for i in lines:
		# start
		if i == 0:
			next_turn = nodes_orientations[i + 1]
			if next_turn >= 0:
				path.append((0, 2))
			elif next_turn < 0:
				side *= -1
				path.append((1, 3))
			cross.append(False)
		# end
		elif i == len(lines) - 1:
			x = path[-1][-1]
			path.append((x - 2, x))
			cross.append(False)
		# regular
		else:
			turns = nodes_orientations[i] * nodes_orientations[i + 1]
			x = path[-1][-1]
			# same side
			if turns == 0:
				if nodes_orientations[i] == 0 and nodes_orientations[i + 1] != 0:
					turns = side * nodes_orientations[i + 1]
				else:
					turns = 1
			if turns > 0:
				path.append((x - 2, x))
				cross.append(False)
			# change side
			elif turns < 0:
				side *= -1
				if x == 2:
					path.append((0, 3))
				elif x == 3:
					path.append((1, 2))
				cross.append(True)

	# alpha = arcsin(offset * 2 / line length)
	# depending on 0-3 or 1-2, negative or positive rotation around initial point
	wound_lines = []
	for i in lines:
		print(path)
		a = line_circle_points[i][path[i][0]]
		b = line_circle_points[i][path[i][1]]
		if path[i] == (0, 3) or path[i] == (1, 2):
			a0 = lines[i][0]
			b0 = lines[i][1]
			length = distance_point_point(a0, b0)
			angle = asin(offset * 2 / length)
			print(angle)
			if path[i] == (1, 2):
				angle *= -1
			a = rotate_points([a], angle, axis=[0.0, 0.0, 1.0], origin=a0)[0]
			b = rotate_points([b], angle, axis=[0.0, 0.0, 1.0], origin=b0)[0]
		wound_lines.append((a, b))

	return wound_lines
	
	# return [[line_circle_points[i][path[i][0]], line_circle_points[i][path[i][1]]] for i in lines]


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':

	from compas_plotters.plotter import Plotter

	points = {
		0: [0.0, 0.0, 0.0],
		1: [1.0, 0.0, 0.0],
		2: [2.0, 0.0, 0.0],
		3: [0.0, 1.0, 0.0],
		4: [1.0, 1.0, 0.0],
		5: [2.0, 1.0, 0.0],
		6: [0.0, 2.0, 0.0],
		7: [1.0, 2.0, 0.0],
		8: [2.0, 2.0, 0.0],
	}
	point_indices = [0, 4, 8, 5, 7, 6]
	polyline = [points[i] for i in point_indices]
	offset = .3
	wound_lines = wind_polyline(polyline, offset)

	plot_points = [{'pos': xyz, 'radius': 0.02, 'facecolor': '#ffffff'} for xyz in points.values()]
	plot_lines = [{'start': u, 'end': v, 'width': 1.0} for u, v in wound_lines]

	plotter = Plotter(figsize=(5, 5))
	plotter.draw_points(plot_points)
	plotter.draw_lines(plot_lines)
	plotter.show()
	