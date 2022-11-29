import math
import random
import numpy as np
from sklearn.neighbors import KDTree
from sklearn.neighbors import DistanceMetric


class BoxFinder:
	def __init__(self, X, size=1, n_iter=1000, start_step=0.5, min_step=0.1, split_coeff=0.9, verbose=0):
		'''
			size -- given cube size (float)
			start_step -- gradient ascent start step (float)
			min_step -- threshold to stop gradient ascent (float)
			split_coeff -- step decreasing coefficient
		'''
		self.s = size
		self.start_step = start_step
		self.min_step = min_step
		self.split_coeff = split_coeff
		self.verbose = verbose
		self.n_iter = n_iter

		if self.verbose:
			print("")

		self.X = np.array(X)
		self.pts_quant = len(self.X)
		self.dim = len(self.X[0])  # Space dimension

		# Turn s-cube to 1-cube
		self._normalize_space()

		self.KDTree = KDTree(self.X,
								metric=DistanceMetric.get_metric('chebyshev'))


	def _random_starts_quantity(self):
		'''
			How many times randomly generate
			start point near start_points
		'''
		return self.n_iter

	def fit(self):
		# Centers
		base_ctr = self.get_base_cube()
		max_ctr = self._find_max_enclosing_cube(base_ctr)
		# Coverage
		base_cov = self._count_points(base_ctr)
		max_cov = self._count_points(max_ctr)
		if base_cov > max_cov:
			max_ctr = base_ctr
			max_cov = base_cov

		return (max_ctr - 0.5)*self.s, max_cov, (base_ctr - 0.5)*self.s, base_cov

	#def _get_start_points(self):
	def get_base_cube(self):
		'''
			Optimal cube from integer lattice
		'''
		# Create grid and calculate how many points are in each cell
		grid = self._calc_grid()
		# Now need to find starting points for algorithm
		s = sorted(grid.items(), key=lambda x: x[1], reverse=True)
		return np.array(s[0][0])

	def _normalize_space(self):
		'''
			This method turns s-cube to 1-cube
		'''
		self.X /= self.s


	def _calc_grid(self):
		'''
			Creates integer-tick grid (cells are 1-cubes)
		'''
		grid = {}
		for p in self.X:
			# Grid unit cell (cube) center
			cube = tuple(np.floor(p) + 0.5)
			grid[cube] = grid.get(cube, 0) + 1
		return grid

	def _find_max_enclosing_cube(self, center):
		'''
			Do gradient ascenc from randomly generated
			points "near" center
		'''
		max_count = 0
		max_center = None

		if self.verbose:
			msg = "Doing gradient ascent near base cube (it has {} points)"
			print(msg.format(self._count_points(center)))
			print("#########################################")

		for i in range(self._random_starts_quantity()):
			rnd_start = self._random_point_near_center(center)

			if self.verbose:
				msg = "Generating next random point near current start point (it has {} points)"
				print(msg.format(self._count_points(rnd_start)))
				print("****************")

			# Make gradiend ascent from rnd_start
			center_, count = self._gradient_ascent(start=rnd_start)
			if count > max_count:
				max_count = count
				max_center = center_

			if self.verbose:
				print("Local maximum for this point -- {} points".format(count))
				print("****************")

		if self.verbose:
			print("#########################################")

		return np.array(max_center)

	def _random_point_near_center(self, center):
		'''
			Get random point near to center in 1-cube
		'''
		center = np.array(center)
		R = 0.5
		rnd_point = center + (2*np.random.rand(self.dim) - 1)*R
		if self._count_points(rnd_point) != 0:
			return rnd_point
		else:
			# Garantee, that it have > 0 points
			#print("Please decrease min_step")
			pts, dists = self.KDTree.query_radius([center], 0.5, return_distance=True, sort_results=True)
			R = dists[0][0]
			return center + (2*np.random.rand(self.dim) - 1)*R

	def _gradient_ascent(self, start):
		center = np.array(start)
		step = self.start_step

		while 1:
			shift, count = self._get_gradient_shift(center, step)
			if shift is None:
				step *= self.split_coeff
				if step < self.min_step:
					break

				if self.verbose:
					print("Gradient ascent: nowhere to go with current step, decreasing it")
					print("~~~~~~~")
			else:
				center += shift
				if self.verbose:
					print("Gradient ascent: moving to positon with {} points, step -- {}".format(count, step))
					print("~~~~~~~")

		return center, count

	def _get_gradient_shift(self, center, step):
		'''
			Checks 2*dim neighbor cubes and moves to it
			if there is more points, that in center
		'''
		center_count = self._count_points(center)
		shift = np.array([0.]*self.dim)
		max_neighbor = 0
		max_shift = None
		for i in range(self.dim):
			shift[i] = step
			for k in (-1, 1):  # Forward or backward
				shift[i] *= k
				count = self._count_points(center + shift)
				if count > max_neighbor:
					max_neighbor = count
					max_shift = np.array(shift)

			shift[i] = 0
		if max_neighbor > center_count:
			return max_shift, max_neighbor

		return None, center_count

	def _count_points(self, center):
		'''
			Counts how may points are in 1-cube with given center
		'''
		return self.KDTree.query_radius([center], 0.5, count_only=True)[0]
