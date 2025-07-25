import numpy as np
from itertools import combinations_with_replacement as cwr
from math import log2, comb, floor, ceil, factorial
from scipy.integrate import quad
from numpy import inf, exp, pi, sqrt
import mpmath

mpmath.mp.dps = 100
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

#Convolution of probability mass function.
def convolve_pmf(pmf1, pmf2):
    """Convolve two probability mass functions (PMFs)."""
    len1, len2 = len(pmf1), len(pmf2)
    result = np.zeros(len1 + len2 - 1)
    for i in range(len1):
        for j in range(len2):
            result[i + j] += pmf1[i] * pmf2[j]
    return result

	
"""---------------------------------------------"""

def dual_cost(k, m, p, z, t, delta1, delta2, delta3, u, sei, verb):
	"""
	DEPTH3
	z = 14
	E = {(-2)^i, i =0,...,13} 
	
	----------------params------------------
	k: secret length
	m: oracle calls
	p: field size
	z: restricted set size
	t: t-error sample
	delta1: matching on delta position of t-error sample
	delta2: matching on delta position of 2t-error sample
	delta3:
	u: average weight of covering codes part ( w  = 8t + u)
	sei: square euclidean imbalance
	"""
	global pmf_x
	if (k-delta1 -delta2 - delta3) % 3 != 0: #l =3
		return 6000
	#b = blocks
	b = int((k - delta1-delta2 -delta3)/3) #l'
	if verb:
		print("b = ", b)
	
	
	#t-error sample list. 8 Lists
	L_t = (z)**t * comb(int(m),t)
	#2t-error sample list. 4 lists. Devided by z due to dependencies.
	L_2t = L_t ** 2 * p**(-delta1) /z
	#4t- error sample list. 2 list
	L_4t = L_2t ** 2 * p **(-delta2)
	
	#Final list, 8t-error sample
	L = L_4t ** 2 * p **(-delta3)
	
	#Cost of building tree
	C_t = 8* L_t * (k * log2(p) + delta1 * log2(p))
	C_2t = 4*L_2t * (k*log2(p) + delta2 * log2(p)) 
	C_4t = 2*L_4t * (k*log2(p) + delta3 * log2(p))
	
	
		
	C_coll = L * k * log2(p)
	C_sample = C_t + C_2t + C_4t + C_coll
	
	if verb:
		print('delta1 prob = 2^', log2(p**(-delta1)))
		print('delta2 prob = 2^', log2(p**(-delta2)))
		print('delta3 prob = 2^', log2(p**(-delta3)))
		print('t list= 2^',log2(L_t))
		print('2t list = 2^',log2(L_2t))
		print('4t list = 2^',log2(L_4t))
		print('L = 2^', log2(L) )
		print('building t-lists = 2^', log2(C_t))
		print('building 2t-lists = 2^', log2(C_2t))
		print('building 4t-lists = 2^', log2(C_4t))
		print('Final list cost = 2^', log2(C_coll))
		print('create_sample() cost = 2^',log2(C_sample))
	
	#Covering, 3 as the block code length
	#block length = 3
	C_cover = 2 * 3 * log2(p) * p**4 + b * L
	if verb:
		print('Covering codes', log2(C_cover))
	#Amount of samples kept for DFT
	#probability of having weight exactly u and probability
	# Convolve the PMF u times
	
	pmf_y = pmf_x  # Start with the PMF of X
	for _ in range(b - 1):
    		pmf_y = convolve_pmf(pmf_y, pmf_x)
    		
	L_sht = L * pmf_y[u] #pmf_y[u] is the probability of error having weight exactly u.  
	
	#C_sht = p**b * log2(p**b) 
	C_sht = p**(b+1) * (b+1) * log2(p) + p**(b + 1)
	if verb:
		print('ratio of kept samples', pmf_y[u])
		print('DFT',log2(C_sht) )
	
	#Success probability
	
	d = float(L_sht * sei) 
	if d < 50: #d = 50 gives roughly 9x10^-6 success probability (typically 7^8 guess).
		print("not enough vector")
		return 6001
	#Approximation, only use when d is large, e.g., d > 75
	#P_success =  np.exp(-z**b * np.exp(-d/4)/np.sqrt(2*pi))
	P_success = 1.0-  quad(integrand, -inf, -sqrt(d/2.0))[0]/sqrt(2.0*pi) #prob that the correct guess score higher than a wrong guess
	P_success = P_success**(z**b -1) #prob that the correct guess having the best score
	if verb:
		print('Sample for guessing = 2^', log2(L_sht))
		print('ratio = ', d)
		print('sucess probability = 2^', log2(P_success))
	lists =[log2(8*m), log2(L_t), log2(L_2t), log2(L_4t), log2(L)]
	memory_cost = max(lists)
	if verb:
		print('memory_access_cost', memory_cost)
	complexity = (C_sample + C_cover + C_sht)*memory_cost/P_success
	return log2(complexity)	


def opt_dual(k,p, z, delta1_below, delta2_below,delta3_below,exp_below, exp_above):
	
	#bjmm_shifted complexity
	complexity = 155
	weight_per_block = 3.95
	params = {}
	mem = 0
	global E

	for delta1 in range (delta1_below,13):
		for delta2 in range(delta2_below,delta1+1):
			for delta3 in range(delta3_below, delta2+ 1):
				print("deltas", delta1, delta2, delta3)
				if (k- delta1 - delta2 - delta3) % 3 != 0:
					#print('no')
					continue
				num_blocks = (k-delta1 - delta2 - delta3) / 3
				dft = p**num_blocks * log2(p ** num_blocks)

				if log2(dft) >= complexity:
					#print('big')
					continue
			
				for t in range(3,8): #First half (7-11), later half (2,7)
					print("t = ",t)
					u = int(round(weight_per_block * (k-delta1 -delta2- delta3)/3))
					w = u + 8*t
				
				#Estimate SEI
				
					P = [0 for _ in range(p)]
					for i in range(int(z/2)):
						P[2**i] = 1./z
						P[p-2**i] = 1./z
					for _ in range(w-1):
						P = finite_field_convolution(P,E,p)
					sei = square_euclidean_imbalance(P,p)
					
					print("sei =", log2(sei))
					if -log2(sei) > complexity:
						print('sei too small')
						break
					#Test if t is too small
					m_max = 2**exp_above
					#L_tmax = z**t * comb(int(m_max),t)
					L_tmax = log2(z**t * comb(int(m_max),t))
					#L_2tmax = L_tmax **2 * p**(-delta1)
					L_2tmax = 2 * L_tmax - delta1*log2(p) / z
					#L_4tmax = L_2tmax **2 * p **(-delta2)
					L_4tmax = 2 * L_2tmax - delta2*log2(p)
					#L_max = L_4tmax ** 2 * p ** (-delta3)
					L_max = 2* L_4tmax - delta3*log2(p)
					#L_shtmax = L_max * comb(u, int(u/2)) * 2**(-u)
					L_shtmax = L_max + log2(comb(u, int(u/2))) - u 
					#d_max = float(L_shtmax * sei)
					d_max = L_shtmax + log2(sei)
					if d_max < log2(50):
						print("t too small, increase")
						continue
					
					exp = exp_below
					while True:
						m = 2**exp
						L_t = z**t * comb(int(m),t)
						L_2t = L_t ** 2 * p**(-delta1)
						L_4t = L_2t ** 2 * p**(-delta2)
						L = L_4t ** 2 * p**(-delta3)
						#skip the other prob for now, too cumbersome.
						L_sht = L * comb(u, int(u/2)) * 2**(-u)
						d = float(L_sht * sei)
						if d < 30:
							exp += 0.2
						else:
							break
					print("exp_start", exp)
					#Also test if t is too big
					m_test = 2**exp
					L_test = z**t * comb(int(m_test),t)
					if log2(L_test) > complexity:
						print("Building cost too much, skip")
						break
					while exp <= exp_above:
						m = 2**exp
						temp = dual_cost(k, m, p, z, t, delta1, delta2, delta3, u, sei, verb = False)
						#print("temp", temp)
						if temp < complexity:
							complexity = temp
							params = {"delta1": delta1, "delta2": delta2, "delta3": delta3, "t": t, "u" : u, "w": w, "sei": log2(sei), "m": log2(8*m)}
							mem = log2(8*m)
						else:
							exp += 0.01
	
	return complexity, mem, params
if __name__=='__main__':
	p = 127
	z = 14 #bigger E set
	k = 50 #shorter secret_length
	E = [mpmath.mpf(0) for _ in range(p)]
	for i in range(int(z/2)):
		E[2**i] = mpmath.mpf(1)/mpmath.mpf(z)
		E[p-2**i] = mpmath.mpf(1)/mpmath.mpf(z)
	P = [0 for _ in range(p)]
	for i in range(int(z/2)):
		P[2**i] = 1./z
		P[p-2**i] = 1./z
	NIST1 = 150
	NIST2 = 250		
	
	#Number of oracle called needed for algebraic attacks
	max_calls = 0
	for i in range(1,z+1):
		max_calls += comb(k,i)
	max_exp = int(floor(log2(max_calls)))
	print('Oracle calls for algebraic attacks', max_exp)
	x_values = [0,1,2, 3, 4, 5, 6]  # Possible values of error weight
	probabilities = [6.20001e-05, 0.00260401, 0.0364561, 0.218736, 0.502573, 0.226548, 0.01302]  # Corresponding probabilities

	# Create the initial PMF of X
	max_x = max(x_values)
	pmf_x = np.zeros(max_x + 1)
	for i, prob in zip(x_values, probabilities):
		pmf_x[i] = prob
	
	
	C = []
	M = []
	PARAM = []
	
	exp_below =  13
	delta1_below = 6
	delta2_below = 3
	delta3_below = 3
	

	for exp_above in range(15,max_exp - 5,1): 
		print("exp_below",exp_below)
		print("exp_above",exp_above)
		print("delta1_below", delta1_below)
		print("delta2_below", delta2_below)
		print("delta3_below", delta3_below)
		comp, mem, params = opt_dual(k, p, z,delta1_below, delta2_below, delta3_below, exp_below,exp_above)
		print("complexity:", comp)
		print("mem:" , mem)
		print("params:", params)
		
		if comp < NIST1:
			exp_below = mem - 3
			#delta1_below = params["delta1"]
			#delta2_below = params["delta2"]
			#delta3_below = params["delta3"]2
			if round(comp,2) not in C:
				C.append(round(comp,2))
				M.append(round(mem,2))
				PARAM.append([params["delta1"], params["delta2"], params["delta3"], params["t"]])
				with open("output1.txt", "a") as file:
    				# Writing the C list
					file.write("C: " + str(round(comp,2)) + "\n")
    
    				# Writing the M list
					file.write("M: " + str(round(mem,2)) + "\n")
    
    				# Writing the PARAMS list
					file.write("PARAMS:\n")
					extract = ["delta1", "delta2", "delta3", "t"]
					for key in extract:
						file.write(str(params[key]) + ", ")
					file.write("\n")
			
	print(C)
	print(M)
	print(PARAM)
	exit()
	
	#increment for memory is 0.01
	 #Test SEI, with E \cup -E
	upper_W = 56
	print('weight=', upper_W)
	for _ in range(upper_W - 1):
		P = finite_field_convolution(P, E, p)
	SEI = square_euclidean_imbalance(P, p)
	print("SEI: ", log2(SEI))
	print('-----------------')
	print(dual_cost(k=50, m = 2**(19.9), p = 127, z = 14, t = 3, delta1 = 9, delta2 = 9, delta3 = 8, u = 32, sei = SEI, verb = True))

