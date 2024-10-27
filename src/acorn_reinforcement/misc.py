from __future__ import print_function
from __future__ import absolute_import
from __future__ import division


__all__ = []


def remove_pair_in_dict(dictionary, pair):
	u, v = pair
	if (u, v) in dictionary:
		del dictionary[(u, v)]
	elif (v, u) in dictionary:
		del dictionary[(v, u)]
	else:
		print('missing', (u, v))


def json_lines(lines):
	return [xyz for line in lines for point in line for xyz in point]


def unjson_lines(lines):
	lines_2 = []
	line = []
	point = []
	for i, xyz in enumerate(lines):
		point.append(xyz)
		if i % 3 == 2:
			line.append(point)
			point = []
		if i % 6 == 5:
			lines_2.append(line)
			line = []
	return lines_2
