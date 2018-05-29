# Find examples where the function: sum of the absolute values of the
# greatest two eigenvalues, is nonconvex for general (nonsymmetric)
# matrices.

# Such an example would show that we can't hope to minimize this sum
# in a semidefinite program because SDPs are limited to expressing
# convex minimization. Thus we would either have to find a convex 
# relaxation of the original problem, or find a more suitable model.

# I haven't found any yet.

import numpy as np
from scipy.optimize import brentq

# Use a variant of Webster to round off the values of a vector so that
# 	- each value has only a few significant digits (easier to transcribe)
#	- it all sums up to 1
def normalize_to_sign_digits(vector, significant_digits):
	round_factor = 10 ** significant_digits
	# We now need to find some factor x so sum(round(vector*x)) = 
	# round_factor
	f = lambda x: np.sum(np.round(vector*x))-round_factor
	vecsum = np.sum(vector)

	lower_bound = (round_factor - 0.5*len(vector))/vecsum
	upper_bound = (round_factor + 0.5*len(vector))/vecsum

	x = brentq(f, lower_bound, upper_bound)

	return np.round(vector*x)/round_factor

# Note, the result isn't uniformly distributed on the (n-1) simplex,
# but for our purposes, that doesn't matter.
def get_normalized_vector(length, significant_digits):
	round_factor = 10 ** significant_digits

	random_unif_vector = np.random.rand(1, length).flatten()
	return normalize_to_sign_digits(random_unif_vector,
		significant_digits)

def get_markov_matrix(rows, cols):
	vectors = [[get_normalized_vector(cols, 2) for x in xrange(rows)]]
	return np.array(vectors)

def eigenvalues(matrix):
	v, w = np.linalg.eig(matrix)
	# Get absolute value of eigenvalues, and sort them
	v = -np.abs(v.flatten())
	v.sort()
	return -v 

def sum_of_k_eigenvalues(matrix, k):
	return np.sum(eigenvalues(matrix)[0:k-1])

k = 2

record_violation_magnitude = 1e-2

for i in xrange(10**18):

	dim = np.random.randint(2, 10)

	A = get_markov_matrix(dim, dim)
	B = get_markov_matrix(dim, dim)
	lambda_var = np.round(np.random.random()*100)/100.0
	convex_AB = (lambda_var) * A + (1 - lambda_var) * B

	f_of_combination = sum_of_k_eigenvalues(convex_AB, k)
	combination_of_f = lambda_var * sum_of_k_eigenvalues(A, k) + \
						(1-lambda_var) * sum_of_k_eigenvalues(B, k)

	# If this is positive, we have proof the function is not convex.
	violation_magnitude = f_of_combination - combination_of_f

	if violation_magnitude > record_violation_magnitude:
		record_violation_magnitude = violation_magnitude
		print "Convexity violation found!"
		print record_violation_magnitude
		print A
		print B
		print convex_AB
		print f_of_combination, combination_of_f
		print lambda_var
		print "Eigenvalues:"
		print "A:", eigenvalues(A)
		print "B: ", eigenvalues(B) 
		print "combination: ", eigenvalues(convex_AB)
