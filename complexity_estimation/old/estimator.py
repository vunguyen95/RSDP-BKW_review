import numpy as np
from itertools import combinations_with_replacement as cwr
from math import log2, comb, floor, ceil, factorial
from scipy.integrate import quad
from numpy import inf, exp, pi, sqrt
import mpmath

mpmath.mp.dps = 100
def list_size(z, length):
	return z**length

def mem(z, length):
	return  list_size(z,length) * length * log2(z)
	
def list_cost(p, length, z, l):
	return list_size(z,length) * (log2(z)*length + l*log2(p))

def num_sol(p, z, error_length , nbr_parity):
	return 1 + z**(error_length) * p**(-nbr_parity)
def coll_cost(L1, L2, p, l, secret_length):
	return L1* L2 * p**(-l) * (secret_length+l) * log2(p)

def stern(p,z,secret_length, nbr_samples, MODE):
	#solutions = num_sol(p, z, secret_length + nbr_samples, nbr_samples)
	solutions = 1
	min = 10000
	cost_l = 0
	l_min = 0
	for l in range (nbr_samples+1):
		len1 = floor((secret_length+l)/2)
		len2 = ceil((secret_length+l)/2)
		assert len1 + len2 == secret_length + l
		#List 1
		L1 = list_size(z, len1)
		cost_L1 = list_cost(p, len1, z, l)* len1 * log2(z) 
		#List 2
		L2 = list_size(z, len2)
		cost_L2 = list_cost(p, len2, z, l)
		#Memory requirement
		MEM = mem(z, len1)
		
		C_bitop = (cost_L1 + cost_L2 + coll_cost(L1, L2, p, l, secret_length)) / solutions
		if MODE == 0:
			C = log2(C_bitop)
		if MODE == 1:
			C = log2( C_bitop * log2(MEM))
			 
		if C < min:
			min = C
			cost_l = cost_L1
			l_min = l
	return [min,log2(cost_l), l_min]
	

		
def count_arrangements(W, n, w):
    # Initialize a 2D array for dynamic programming
    dp = [[0] * (n + 1) for _ in range(W + 1)]

    # Base cases
    for j in range(n + 1):
        dp[0][j] = 1

    for i in range(1, W + 1):
        for j in range(1, n + 1):
            for k in range(0, min(i, w) + 1):
                dp[i][j] += dp[i - k][j - 1]

    return dp[W][n]

def probability_arrangements(W, n, w, probabilities):
    dp = [[[0.0] * (n + 1) for _ in range(W + 1)] for _ in range(len(probabilities) + 1)]

    for j in range(n + 1):
        dp[0][0][j] = 1.0

    for i in range(1, W + 1):
        for k in range(1, len(probabilities) + 1):
            for j in range(1, n + 1):
                for l in range(0, min(i, w) + 1):
                    if j - 1 >= 0 and i - l >= 0:
                        dp[k][i][j] += dp[k - 1][i - l][j - 1] * probabilities[k - 1][l]

    total_prob = 0.0
    for l in range(0, min(W, w) + 1):
        if n - 1 >= 0 and W - l >= 0:
            total_prob += dp[len(probabilities)][W - l][n - 1] * probabilities[-1][l]

    return total_prob
from itertools import combinations_with_replacement

def generate_combinations(n, W, w, probabilities):
    # Generate combinations of length 'n' without replacement
    valid_combinations = [combo for combo in cwr(range(0, w + 1), n) if sum(combo) == W]
    total_prob =0.0
    total_count = 0;
    for combo in valid_combinations:
    	denominator = 1
    	prob = 1.0
    	for element in set(combo):
    		count = combo.count(element)
    		denominator *= factorial(count)
    		#prob *= probabilities[element]**count
    		#print(prob)
    	
    	for element in combo:
    		prob *= probabilities[element]
    	permutation = factorial(n) // denominator
    	#print(permutation)
    	total_prob += prob * permutation
    	total_count += permutation
    return total_prob, total_count

   

def integrand(x):
	return(exp(-(x**2)/2))

def finite_field_convolution(P, Q, field_size=127):
    """
    Compute the convolution of two distributions P and Q over a finite field of given size.
    
    Parameters:
    P (list): Distribution P as a list of probabilities. Must have length equal to field_size.
    Q (list): Distribution Q as a list of probabilities. Must have length equal to field_size.
    field_size (int): The size of the finite field, default is 127.
    
    Returns:
    list: The convolution result as a list of probabilities.
    """
    # Initialize an empty list to store the convolution result
    R = [mpmath.mpf(0)] * field_size
    
    # Loop through each element in the finite field
    for x in range(field_size):
        # Calculate the convolution sum for each x
        for y in range(field_size):
            # Perform the modulo operation for (x - y)
            z = (x - y) % field_size
            
            # Update the value of R(x)
            R[x] += mpmath.mpf(P[y]) * mpmath.mpf(Q[z])
            
    total_prob = sum([mpmath.mpf(R[_]) for _ in range(field_size)])
    normalized_dist = [mpmath.mpf(R[_])/total_prob for _ in range(field_size)]
            
    return normalized_dist

def square_euclidean_imbalance(P, field_size):
	"""
    Compute the convolution of two distributions P and Q over a finite field of given size.
    
    Parameters:
    P (list): Distribution P as a list of probabilities. Must have length equal to field_size.
    field_size (int): The size of the finite field.
    
    Returns:
    The Square Euclidean Imbalance of P w.r.t the Uniform distribution."""
    
	U = [mpmath.mpf(1)/field_size] * field_size
	SEI = field_size*sum([(mpmath.mpf(P[i]) - mpmath.mpf(U[i]))**2 for i in range(p)])
	return SEI

def find_minimum_mem(SEI, t, delta, z, p):
	exponent_lower = 5.0
	delta_probability = p**(-delta)
	while True:
		M = 2**(exponent_lower)
		total_t_comb = (2*z + 1)**t * comb(int(M), t)
		samples = total_t_comb*(total_t_comb-1)* delta_probability/2
		if samples < 1.0:
			exponent_lower += 2.0
		elif samples >=  1.0 and log2(samples) < - log2(SEI):
			exponent_lower += 0.10
		else:
			break
	return exponent_lower

def t_comb(M, z, t):
	return (2*z + 1)**t * comb(int(M), t)

def delta_samples(list_size, delta_probability):
	return list_size*(list_size -1) * delta_probability/2
	
if __name__=='__main__':
	p = 127
	z = 7
	secret_length = 76
	M_stern = 51
	stern_complexity = stern(p,z,secret_length,M_stern,0)
	stern_bitop = stern_complexity[0]
	stern_mem = stern_complexity[1]
	print("Stern Complexity: ", stern_bitop, "/ Stern Memory: ", stern_mem)
	
	""" To find the best parameters for the new approach, first we loop through values of delta
	 such that cost_DFT < stern_bitop"""
	
	delta_lower = secret_length % 3

	while log2(p**( (secret_length - delta_lower)/3) * log2( p**( (secret_length - delta_lower)/3))) >= stern_bitop:
		delta_lower += 3
	
	print("Upper bound of delta: ", delta_lower)
	
	upper_W = 95
	"""E = [mpmath.mpf(0) for _ in range(p)]
	for i in range(7):
		E[2**i] = mpmath.mpf(1)/mpmath.mpf(7)
	P = [0 for _ in range(p)]
	for i in range(7):
		P[2**i] = 1./7
	for _ in range(upper_W - 1):
		P = finite_field_convolution(P, E, p)
	SEI = square_euclidean_imbalance(P, p)
	print("SEI: ", log2(SEI))"""
	
	print("Upper_W: ", upper_W)
		
	
	"""
	
	E = [mpmath.mpf(0) for _ in range(p)]
	for i in range(7):
		E[2**i] = mpmath.mpf(1)/mpmath.mpf(7)
	params = [0,0,0,0,0]
	complexities = [stern_bitop, stern_mem, 1]
	for delta in range(delta_lower, delta_lower + 9, 3): #could be from delta_lower to secret_length
		 #compute the upper bound for t so that -log2(SEI) < stern_cost 
		 complexity_DFT = p**( (secret_length - delta)/3) * log2( p**( (secret_length - delta)/3))
		 upper_t = floor((upper_W - ceil(3.93*(secret_length-delta)/3))/2)
		 
		 #delta probability
		 probability_delta = p**(-delta)
		 
		 #loop through different values of t
		 for t in range(2, upper_t + 1):
		 	#compute W:
		 	W = int(ceil(3.93*(secret_length - delta))/3 + 2*t)
		 	
		 	# compute SEI for each W
		 	P = [0 for _ in range(p)]
		 	for i in range(7):
		 		P[2**i] = 1./7
		 	for _ in range(W - 1):
		 		P = finite_field_convolution (P, E, p)
		 	SEI = square_euclidean_imbalance(P, p)
		 	
		 	#Get the minimum values for memory
		 	minimum_mem_exp = find_minimum_mem(SEI, t, delta, z, p)
		 	M = 2**(minimum_mem_exp)
		 	total_t_comb = t_comb(M,z, t)
		 	complexity_t_comb = total_t_comb*((secret_length)*log2(p) + M*log2(z)  + delta*log2(p))
		 	samples = delta_samples(total_t_comb, probability_delta)
		 	complexity_samples = samples*(secret_length) * log2(p) + M
		 	 
		 	#overall complexity
		 	complexity = log2(complexity_samples + complexity_t_comb + complexity_DFT)
		 	mem = log2(total_t_comb*((secret_length)*log2(p) + M*log2(z) ))
		 	ratio = samples*SEI
		 	success_probability = 1- 2 * quad(integrand, -inf, -sqrt(ratio)/2)[0]/sqrt(2*pi)
		 	complexity += log2(1/success_probability)
		 	if complexity < complexities[0]:
		 		complexities[0] = complexity
		 		complexities[1] = mem
		 		complexities[2] = success_probability
		 		params[0] = W
		 		params[1] = log2(SEI)
		 		params[2] = minimum_mem_exp
		 		params[3] = delta
		 		params[4] = t
	
	print("Complexity: ", complexities[0], "Memory: ", complexities[1], "Success probability: ", complexities[2])
	print("W: ", params[0], "SEI: ", params[1], "queries: ", params[2], "delta: ", params[3], "t: ", params[4])	 	
	"""
		 	
		 	
	M = 2**(12.6)
	delta = 19
	t = 9
	W = ceil(3.93*(secret_length - delta))/3 + 2*t
	
	print("W: " , W)
	
	print("delta probability:",log2(p**-delta))
	
	#print("count:" , count_arrangements(W-2*t, 15, 3))
	total_t_comb = t_comb(M,z, t)
	print("total t combinations: ", log2(total_t_comb))
	
	complexity_t_comb = log2(total_t_comb*((secret_length)*log2(p) + delta*log2(p)))
	print("complexity t combinations: ", complexity_t_comb)

	
	total_2t_comb = total_t_comb*(total_t_comb-1)/2  
	
	#print("2t created samples: ", log2(total_2t_comb))
	
	
	
	
	samples = delta_samples(total_t_comb, p**(-delta))
	print("samples: ", log2(samples))
	complexity_samples = log2(samples*(secret_length) * log2(p) * secret_length )
	print ("complexity samples: ", complexity_samples)
	Memory = log2(total_t_comb*((secret_length)*log2(p) + M ))
	print("MEMORY: " , Memory)
	print ("DFT: ", log2(p**( (secret_length - delta)/3) * log2( p**( (secret_length - delta)/3))))
	
	ratio = samples* 2**(-125.34)
	print("PE:  ",   quad(integrand, -inf, -sqrt(ratio)/2)[0]/sqrt(2*pi))
	print("Success probability:  ",  1- 2 * quad(integrand, -inf, -sqrt(ratio)/2)[0]/sqrt(2*pi))
	#print(log2((2*z + 1)**t))
	
	
	
	
	"""
	# Example usage
	n = 15 # Number of integers in each combination
	s = 22  # Target sum
	w = 3
	probabilities = [1.0/127.0,14.0/127.0,42.0/127.0,70.0/127.0]
	result = generate_combinations(n, s, w, probabilities)

	print(f"Combinations of {n} positive integers that sum up to {s}:")
	print(log2(result[0]), result[1]) 
	
	print(set([1,2,1]))"""
