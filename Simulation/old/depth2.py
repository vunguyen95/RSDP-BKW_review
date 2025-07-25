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


def dual_cost(k, m, p, z, t, delta, u, sei, verb):
	"""
	Cost of the dual approach:
	k: secret length
	m: oracle calls
	p: field size
	z: restricted set size
	t: t-error sample
	delta: matching on delta position of t-error sample
	u: average weight of covering codes part ( w  = 2t + u)
	sei: square euclidean imbalance
	"""
	if (k-delta) % 3 != 0:
		return 600
	b = (k - delta)/3
	#b = 19
	
	#t-error sample list
	L_t = z**t * comb(int(m),t)
	#Cost of building t-error sample list, sorting by delta position
	C_t = 2* L_t * (k * log2(p) + delta * log2(p))
	
	#2t-error sample, throw away L_t useless, all zero samples, also assuming two batches
	#of M samples
	L = L_t ** 2 * p**(-delta) 
	
	C_coll = L * k * log2(p)
	C_sample = C_t + C_coll
	
	if verb:
		print('Building lists',log2(C_sample))
		print('building t-lists', log2(C_t))
		print('collision', log2(C_coll))
		print('delta prob = 2^', log2(p**(-delta)))
		print('t_list= 2^',log2(L_t))
		print('2t_list',log2(L))
	if L < 0:
		return 601
	#Covering, 3 as the block code length
	
	C_cover = 2 * 3 * log2(p) * p**4 + b * L
	if verb:
		print('Covering codes', log2(C_cover))
	#Amount of samples kept for DFT
	L_sht = L * comb(u, int(u/2)) * 2**(-u) 
	
	C_sht = p**b * log2(p**b)
	if verb:
		print('ratio of kept samples', comb(u, int(u/2)) * 2**(-u))
		print('DFT',log2(C_sht))
	
	#Success probability
	d = float(L_sht * sei)
	P_success =  1.0- 2.0 * quad(integrand, -inf, -sqrt(d)/2.0)[0]/sqrt(2.0*pi)
	if verb:
		print('ratio = ', d)
		print('sucess probability =', P_success)
	complexity = (C_sample + C_cover + C_sht)/P_success
	return log2(complexity)

def dual_cost_v1(k, m, p, z, t, delta, u, sei, verb):
	"""
	This version, we in the create_sample step, we multiply each sample with E\cup -E. Only keep the desired pattern
	Cost of the dual approach:
	k: secret length
	m: oracle calls
	p: field size
	z: restricted set size
	t: t-error sample
	delta: matching on delta position of t-error sample
	u: average weight of covering codes part ( w  = 2t + u)
	sei: square euclidean imbalance
	"""
	if (k-delta) % 3 != 0:
		return 600
	b = (k - delta)/3
	#b = 19
	
	#t-error sample list. Here we go from z**t to (2*z)**t
	L_t = (2*z)**t * comb(int(m),t)
	
	#Cost of building t-error sample list, sorting by delta position
	C_t = 2* L_t * (k * log2(p) + delta * log2(p))
	
	#2t-error sample, also assuming two batches
	#of M samples
	# Now, we only pick the desired pattern, for example, only keep equal contribution from E and -E
	#Each 2t sample contains 2t coefficients from E and -E, so on average, we keep comb(2t,t)*2**(-2t)
	L = L_t ** 2 * p**(-delta) * comb(2*t, t) * 2**(-2*t)
	
	
	C_coll = L * k * log2(p)
	C_sample = C_t + C_coll
	
	if verb:
		print('Building lists',log2(C_sample))
		print('building t-lists', log2(C_t))
		print('collision', log2(C_coll))
		print('delta prob = 2^', log2(p**(-delta)))
		print('t_list= 2^',log2(L_t))
		print('2t_list',log2(L))
		print('L ratio', comb(2*t, t) * 2**(-2*t))
	if L < 0:
		return 601
	#Covering, 3 as the block code length
	
	C_cover = 2 * 3 * log2(p) * p**4 + b * L
	if verb:
		print('Covering codes', log2(C_cover))
	#Amount of samples kept for DFT
	L_sht = L * comb(u, int(u/2)) * 2**(-u) 
	
	C_sht = p**b * log2(p**b)
	if verb:
		print('ratio of kept samples', comb(u, int(u/2)) * 2**(-u))
		print('DFT',log2(C_sht))
	
	#Success probability
	d = float(L_sht * sei)
	P_success =  1.0- 2.0 * quad(integrand, -inf, -sqrt(d)/2.0)[0]/sqrt(2.0*pi)
	if verb:
		print('ratio = ', d)
		print('sucess probability =', P_success)
	complexity = (C_sample + C_cover + C_sht)/P_success
	return log2(complexity)
def opt_dual(k,p, z, delta_lower, exp_below, exp_above):
	complexity = 300
	#print("exp_above: ", exp_above)
	weight_per_block = 3.95
	params = {}
	mem = 0
	for delta in range (delta_lower, delta_lower + 4,3):
	#The general trend is, the larger M, the smaller t (other wise C_sample is high, no need to iterate through large t)
		for t in range (2,6): # t = 2,3,4,5 for extremely large m
			#print("delta, t: ", delta, t)
			u = int(round(weight_per_block*(k-delta)/3))
			w = u + 2*t
			#print("u:", u)
			#print("w:", w)
			
			# Estimate the SEI
			global E, EN
			plus = int(w/2)
			minus = w - plus
			P = [0 for _ in range(p)]
			for i in range(7):
				P[2**i] = 1./7
			for _ in range(plus-1):
				P = finite_field_convolution(P,E,p)
			for _ in range(minus):
				P = finite_field_convolution(P,EN,p)
			sei = square_euclidean_imbalance(P,p)
			
			#Test if t is too small
			m_max = 2**exp_above
			L_tmax = z**t * comb(int(m_max),t)
			L_max = L_tmax**2 * p **(-delta)
			L_shtmax = L_max * comb(u, int(u/2)) * 2**(-u)
			d_max = float(L_shtmax * sei)
			if d_max < 0.1:
				#print("sei:",log2(sei))
				#print("m_max:", log2(m_max))
				#print("d_max:", d_max)
				continue
			# if t_max is still too small, it return complexity, no params, could be error.
				
			# Find the reasonable lower bound for m
			exp = exp_below
			while True:
				m = 2**exp
				L_t = z**t * comb(int(m),t)
				L = L_t ** 2 * p**(-delta)
				L_sht = L * comb(u, int(u/2)) * 2**(-u)
				d = float(L_sht * sei)
				if d < 0.001:
					#print(d)
					exp += 0.2
				else:
					break
			#Start iterating from exp to exp_above
			while exp <= exp_above:
				m = 2**exp
				temp = dual_cost(k,m, p, z, t, delta, u, sei, verb = False)
				if temp < complexity:
					complexity = temp
					params = {"delta": delta, "t": t, "u": u, "w": w, "sei": log2(sei), "m": log2(m)}
					mem = exp + 1
				else:
					exp += 0.01
	return complexity, mem, params
	
def opt_dual_v1(k,p, z, delta_lower, exp_below, exp_above):
	complexity = 300
	#print("exp_above: ", exp_above)
	weight_per_block = 3.95
	params = {}
	mem = 0
	for delta in range (delta_lower, delta_lower + 4,3):
	#The general trend is, the larger M, the smaller t (other wise C_sample is high, no need to iterate through large t)
		for t in range (2,11): # t = 2,3,4,5 for extremely large m
			#print("delta, t: ", delta, t)
			u = int(round(weight_per_block*(k-delta)/3))
			w = u + 2*t
			#print("u:", u)
			#print("w:", w)
			
			# Estimate the SEI
			global E, EN
			plus = int(w/2)
			minus = w - plus
			P = [0 for _ in range(p)]
			for i in range(7):
				P[2**i] = 1./7
			for _ in range(plus-1):
				P = finite_field_convolution(P,E,p)
			for _ in range(minus):
				P = finite_field_convolution(P,EN,p)
			sei = square_euclidean_imbalance(P,p)
			
			#Test if t is too small
			m_max = 2**exp_above
			L_tmax = (2*z)**t * comb(int(m_max),t)
			L_max = L_tmax**2 * p **(-delta) * comb(2*t,t) * 2**(-t)
			L_shtmax = L_max * comb(u, int(u/2)) * 2**(-u)
			d_max = float(L_shtmax * sei)
			if d_max < 0.1:
				#print("sei:",log2(sei))
				#print("m_max:", log2(m_max))
				#print("d_max:", d_max)
				continue
			# if t_max is still too small, it return complexity, no params, could be error.
				
			# Find the reasonable lower bound for m
			exp = exp_below
			while True:
				m = 2**exp
				L_t = (2*z)**t * comb(int(m),t)
				L = L_t ** 2 * p**(-delta) * comb(2*t,t) * 2**(-t)
				L_sht = L * comb(u, int(u/2)) * 2**(-u)
				d = float(L_sht * sei)
				if d < 0.001:
					#print(d)
					exp += 0.2
				else:
					break
			#Start iterating from exp to exp_above
			while exp <= exp_above:
				m = 2**exp
				temp = dual_cost_v1(k,m, p, z, t, delta, u, sei, verb = False)
				if temp < complexity:
					complexity = temp
					params = {"delta": delta, "t": t, "u": u, "w": w, "sei": log2(sei), "m": log2(m)}
					mem = exp + 1
				else:
					exp += 0.01
	return complexity, mem, params
	
"""---------------------------------------------"""
def dual_cost_v2(k, m, p, z, t, delta1, delta2, u, sei, verb):
	"""
	This version, we do many create_samples() step to keep the error profile more balance.
	Assume two steps, delta1, delta2,
	 
	
	----------------params------------------
	k: secret length
	m: oracle calls
	p: field size
	z: restricted set size
	t: t-error sample
	delta1: matching on delta position of t-error sample
	delta2: matching on delta position of 2t-error sample
	u: average weight of covering codes part ( w  = 4t + u)
	sei: square euclidean imbalance
	"""
	if (k-delta1 -delta2) % 3 != 0:
		return 600
	#b = blocks
	b = int((k - delta1-delta2)/3)
	#b = 19
	
	#t-error sample list. 4 list here L1,2,3,4
	L_t = (z)**t * comb(int(m),t)
	#2t-error sample list. 2 list here
	L_2t = L_t ** 2 * p**(-delta1)
	
	#Final list, 4t-error sample
	L = L_2t ** 2 * p **(-delta2)
	
	#Cost of t-error sample list, sorting by delta1 position
	C_t = 4* L_t * (k * log2(p) + delta1 * log2(p))
	#Cost of 2t-error sample list, sorting by delta2 position
	C_2t = 2*L_2t * (k*log2(p) + delta2 * log2(p)) 
	
	
		
	C_coll = L * k * log2(p)
	C_sample = C_t + C_2t + C_coll
	
	if verb:
		print('delta1 prob = 2^', log2(p**(-delta1)))
		print('delta2 prob = 2^', log2(p**(-delta2)))
		print('t list= 2^',log2(L_t))
		print('2t list = 2^',log2(L_2t))
		print('4t list = 2^',log2(L))
		print('building t-lists = 2^', log2(C_t))
		print('building 2t-lists = 2^', log2(C_2t))
		print('Final list cost = 2^', log2(C_coll))
		print('create_sample() cost = 2^',log2(C_sample))
	
	#Covering, 3 as the block code length
	#block length = 3
	C_cover = 2 * 3 * log2(p) * p**4 + b * L
	if verb:
		print('Covering codes', log2(C_cover))
	#Amount of samples kept for DFT
	pmf_y = pmf_x  # Start with the PMF of X
	for _ in range(b - 1):
    		pmf_y = convolve_pmf(pmf_y, pmf_x)
	L_sht = L * comb(u, int(u/2)) * 2**(-u) * pmf_y[u] #pmf_y[u] being the probability of the error weigh u
	

	C_sht = p**(b+1) * (b+1) * log2(p) + p**(b + 1) 
	if verb:
		print('ratio of kept samples', comb(u, int(u/2)) * 2**(-u),  pmf_y[u])
		print('DFT',log2(C_sht) )
	
	#Success probability
	d = float(L_sht * sei) 
	if d < 90: #Depth 2 gives more testing candidates (typically 7**15, hence require lots of samples)
		print("not enough sample")
		break
	P_success = 1.0-  quad(integrand, -inf, -sqrt(d/2.0))[0]/sqrt(2.0*pi) #prob that the correct guess score higher than a wrong guess
	P_success = P_success**(z**b -1) #prob that the correct guess having the best score
	if verb:
		print('ratio = ', d)
		print('sucess probability = 2^', log2(P_success))
	complexity = (C_sample + C_cover + C_sht)/P_success
	return log2(complexity)	



def opt_dual_v2(k,p, z, delta1_below, delta2_below,exp_below, exp_above):
	#bjmm_shifted complexity
	complexity = 150
	#print("exp_above: ", exp_above)
	weight_per_block = 3.95
	params = {}
	mem = 0
	global E, EN
	#Too slow, test in different regimes delta1 = 16,20 and delta2 = 12-16?
	for delta1 in range (delta1_below,20):
		for delta2 in range(delta2_below,delta1+1):
			print("deltas", delta1, delta2)
			if (k- delta1 - delta2) % 3 != 0:
				#print('no')
				continue
			num_blocks = (k-delta1 - delta2) / 3
			dft = p**num_blocks * log2(p ** num_blocks)

			if log2(dft) >= complexity:
				#print('big')
				continue
			
			for t in range(4,12):
				print("t = ",t)
				u = int(round(weight_per_block * (k-delta1 -delta2)/3))
				w = u + 4*t
				
				#Estimate SEI
				
				plus = int(u/2) + 4*t
				minus = w - plus
				P = [0 for _ in range(p)]
				for i in range(7):
					P[2**i] = 1./7
				for _ in range(plus-1):
					P = finite_field_convolution(P,E,p)
				for _ in range(minus):
					P = finite_field_convolution(P,EN,p)
					sei = square_euclidean_imbalance(P,p)
				print("sei =", log2(sei))
				#Test if t is too small
				m_max = 2**exp_above
				L_tmax = z**t * comb(int(m_max),t)
				L_2tmax = L_tmax **2 * p**(-delta1)
				L_max = L_2tmax **2 * p **(-delta2)
				L_shtmax = L_max * comb(u, int(u/2)) * 2**(-u)
				d_max = float(L_shtmax * sei)
				if d_max < 80:
					print("t too small, increase")
					continue
					
				exp = exp_below
				while True:
					m = 2**exp
					L_t = z**t * comb(int(m),t)
					L_2t = L_t ** 2 * p**(-delta1)
					L = L_2t ** 2 * p**(-delta2)
					L_sht = L * comb(u, int(u/2)) * 2**(-u)
					d = float(L_sht * sei)
					if d < 80:
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
					temp = dual_cost_v2(k, m, p, z, t, delta1, delta2, u, sei, verb = False)
					#print("temp", temp)
					if temp < complexity:
						complexity = temp
						params = {"delta1": delta1, "delta2": delta2, "t": t, "u" : u, "w": w, "sei": log2(sei), "m": log2(4*m)}
						mem = log2(4*m)
					else:
						exp += 0.02
	
	return complexity, mem, params


if __name__=='__main__':
	p = 127
	z = 7
	k = 76
	E = [mpmath.mpf(0) for _ in range(p)]
	for i in range(7):
		E[2**i] = mpmath.mpf(1)/mpmath.mpf(7)
	# Minus E
	EN = [mpmath.mpf(0) for _ in range(p)]
	for i in range(7):
		EN[p-2**i] = mpmath.mpf(1)/mpmath.mpf(7)
	P = [0 for _ in range(p)]
	for i in range(7):
		P[2**i] = 1./7
	NIST1 = 160
	NIST2 = 250		

	#Number of oracle called needed for algebraic attacks
	max_calls = 0
	for i in range(1,8):
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
	
	exp_below = 10
	delta1_below = 15 #15
	delta2_below = 10 #10
	
	for exp_above in range(13,max_exp - 2,1): #range (exp_below+3, max_exp-2)
		print("exp_below",exp_below)
		print("exp_above",exp_above)
		print("delta1_below", delta1_below)
		print("delta2_below", delta2_below)
		comp, mem, params = opt_dual_v2(k, p, z,delta1_below, delta2_below, exp_below,exp_above)
		print("complexity:", comp)
		print("mem:" , mem)
		print("params:", params)
		
		if comp < NIST1:
			exp_below = mem - 2
			#delta1_below = params["delta1"]
			#delta2_below = params["delta2"]
			#delta3_below = params["delta3"]
			if round(comp,2) not in C:
				C.append(round(comp,2))
				M.append(round(mem,2))
				PARAM.append([params["delta1"], params["delta2"], params["t"]])
			
	print(C)
	print(M)
	print(PARAM)
	exit()
	
	
	
	"""TEST SPACE"""
	
	 #Test SEI, remember to run with E \cup -E later
	
	upper_W = 8 + 4*2
	print('weight=', upper_W)
	plus = int(8/2) + 2*11
	minus = upper_W - plus
	P = [0 for _ in range(p)]
	for i in range(7):
		P[2**i] = 1./7
	for _ in range(plus - 1):
		P = finite_field_convolution(P, E, p)
	for _ in range(minus):
		P = finite_field_convolution(P,EN,p)
	SEI = square_euclidean_imbalance(P, p)
	print("SEI: ", log2(SEI))
	
	#delta = 27: W = 125, -169.78, W = 127, -172.55, W = 119, -161.45
	#delta = 26: W = 115, -155.9
	#print("Upper_W: ", upper_W)
	
	#print(dual_cost(k = 111, m = 2**(39.9), p = 127, z = 7, t = 5, delta = 30, u = 107, sei = 2**(-158), verb = True))
	#print(dual_cost(k = 76, m = 2**(17.45), p = 127, z = 7, t = 7, delta = 19, u = 75, sei = 2**(-119.8), verb = True))
	print('-----------------')
	#print(dual_cost_new(k = 76, m = 2**(16.56), p = 127, z = 7, t = 7, delta = 19, u = 75, sei = SEI, verb = True))
	
	#print(dual_cost_v2(k=76, m = 2**(10.5), p = 127, z = 7, t = 11, delta1 = 18, delta2 = 16, u = 55, sei = SEI, verb = True))
	#print(dual_cost_v3(k=76, m = 2**(19.1), p = 127, z = 7, t = 5, delta1 = 15, delta2 = 14, delta3 = 14, u = 43, sei = SEI, verb = True))

