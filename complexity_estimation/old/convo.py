import numpy as np
import matplotlib.pyplot as plt
import mpmath
from math import log2
#Average error weight: 3.95288
#Error Distribution: 
#Weight: 0: 6.20001e-05; Weight: 1: 0.00260401; Weight: 2: 0.0364561; Weight: 3: 0.218736; Weight: 4: 0.502573; Weight: 5: 0.226548; Weight: 6: 0.01302; Weight: 7: 0; Weight: 8: 0; Weight: 9: 0; Expected weight: 3.95288
#14 42 70
def convolve_pmf(pmf1, pmf2):
    """Convolve two probability mass functions (PMFs)."""
    len1, len2 = len(pmf1), len(pmf2)
    result = np.zeros(len1 + len2 - 1)
    for i in range(len1):
        for j in range(len2):
            result[i + j] += pmf1[i] * pmf2[j]
    return result

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
	SEI = field_size*sum([(mpmath.mpf(P[i]) - mpmath.mpf(U[i]))**2 for i in range(field_size)])
	return SEI
def square_euclidean_imbalance1(P,Q, field_size):
	"""
    Compute the convolution of two distributions P and Q over a finite field of given size.
    
    Parameters:
    P (list): Distribution P as a list of probabilities. Must have length equal to field_size.
    field_size (int): The size of the finite field.
    
    Returns:
    The Square Euclidean Imbalance of P w.r.t Q."""
	SEI = field_size*sum([(mpmath.mpf(P[i]) - mpmath.mpf(Q[i]))**2 for i in range(field_size)])
	return SEI
"""DISTRIBUTION OF ERROR WEIGHT"""
# Define the PMF of X
x_values = [0,1,2, 3, 4, 5, 6]  # Possible values of error weight
probabilities = [6.20001e-05, 0.00260401, 0.0364561, 0.218736, 0.502573, 0.226548, 0.01302]  # Corresponding probabilities

# Create the initial PMF of X
max_x = max(x_values)
pmf_x = np.zeros(max_x + 1)
for i, prob in zip(x_values, probabilities):
    pmf_x[i] = prob

# Convolve the PMF n times
n = 2  # Number of variables to sum
pmf_y = pmf_x  # Start with the PMF of X
for _ in range(n - 1):
    pmf_y = convolve_pmf(pmf_y, pmf_x)

print(pmf_y[8])
# Compute the range of possible values for Y
y_values = np.arange(len(pmf_y))

# Plot the resulting PMF
plt.bar(y_values, pmf_y, alpha=0.7, color="skyblue", edgecolor="black")
plt.xlabel("Sum (Y)")
plt.ylabel("Probability")
plt.title(f"Distribution of Y = sum(X_j), n={n}")
plt.savefig('covno.png', format='png', dpi=300)
#plt.show()
exit()
"""DISTRIBUTION OF NOISE WITH CORRECT GUESS"""
p = 127
E = [mpmath.mpf(0) for _ in range(p)]
for i in range(7):
	E[2**i] = mpmath.mpf(1)/mpmath.mpf(7)
# Minus E
EN = [mpmath.mpf(0) for _ in range(p)]
for i in range(7):
	EN[p-2**i] = mpmath.mpf(1)/mpmath.mpf(7)
upper_W = 8 + 2*2
print('weight=', upper_W)
plus = 6 
minus = upper_W - plus
P = [0 for _ in range(p)]
for i in range(7):
	P[2**i] = 1./7
for _ in range(plus - 1):
	P = finite_field_convolution(P, E, p)
for _ in range(minus):
	P = finite_field_convolution(P,EN,p)
for prob in P:
	print(f"{prob}")
SEI = square_euclidean_imbalance(P, p)
print("SEI: ", log2(SEI))

P_correct = np.loadtxt('guessingd1.txt')
P_random = np.loadtxt('rguessingd1.txt')
#print(P_correct)
print("-------------------")
#print(P_random)
SEI_correct_to_Uniform = square_euclidean_imbalance(P_correct,p)
SEI_random_to_Uniform = square_euclidean_imbalance(P_random,p)
SEI_random_to_Theory = square_euclidean_imbalance1(P_random, P, p)

SEI_correct_to_Theory = square_euclidean_imbalance1(P_correct,P,p)
print("SEI_correct_to_PTheory:" , log2(SEI_correct_to_Theory))
print("SEI_correct_to_PUniform:" , log2(SEI_correct_to_Uniform))
print("SEI_random_to_Uniform:" , log2(SEI_random_to_Uniform))
print("SEI_random_to_Theory:" , log2(SEI_random_to_Theory))
