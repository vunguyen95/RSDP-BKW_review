import numpy as np
from concurrent.futures import ThreadPoolExecutor
import mpmath
from math import log2
# Define the sets
MODULO = 127
E = [mpmath.mpf(0) for _ in range(MODULO)]
for i in range(7):
	E[2**i] = mpmath.mpf(1)/mpmath.mpf(7)
# Minus E
E_minus = [mpmath.mpf(0) for _ in range(MODULO)]
for i in range(7):
	E_minus[MODULO-2**i] = mpmath.mpf(1)/mpmath.mpf(7)

# Parameters
NUM_NUMBERS = 2**15  # Number of numbers to generate

NUM_BINS = 127
NUM_THREADS = 12  # Number of threads
NUM_RUNS = 12  # Number of test runs (one per thread)

# Function to perform one run of the test
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
def perform_test():
    # Initialize the histogram
    counts = np.zeros(NUM_BINS, dtype=np.float64)
    
    # Generate random numbers
    for _ in range(NUM_NUMBERS):
        # Select 6 random elements from E and E_minus
        sum_E = np.sum(np.random.choice(E, size=6, replace=True))
        sum_E_minus = np.sum(np.random.choice(E_minus, size=6, replace=True))
        
        # Compute the modulo
        random_number = (sum_E + sum_E_minus) % MODULO
        
        # Update the histogram
        counts[random_number] += 1

    # Normalize to get a probability distribution
    return counts / NUM_NUMBERS

# Multithreaded execution
def compute_average_distribution():
    with ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        # Run the test multiple times in parallel
        results = list(executor.map(perform_test, range(NUM_RUNS)))
    
    # Average the distributions
    average_distribution = np.mean(results, axis=0)
    return average_distribution

# Main execution
if __name__ == "__main__":
    #final_distribution = compute_average_distribution()
    Psim = perform_test()
    
    # Print the results
    for i, prob in enumerate(Psim):
        print(f"Value {i}: {prob:.6f}")
    
    upper_W = 8 + 2*2
    print('weight=', upper_W)
    plus = 6 
    minus = upper_W - plus
    P = [0 for _ in range(MODULO)]
    for i in range(7):
        P[2**i] = 1./7
    for _ in range(plus - 1):
        P = finite_field_convolution(P, E, MODULO)
    for _ in range(minus):
        P = finite_field_convolution(P,E_minus,MODULO)
    SEI = square_euclidean_imbalance(P, MODULO)
    print("SEI of P_theory and PU:", log2(SEI))
    SEI1 = square_euclidean_imbalance1(Psim, P, MODULO)
    print("SEI of P_sim and P_theory:", log2(SEI1))

"""
    # Optionally, visualize the distribution
    import matplotlib.pyplot as plt
    plt.bar(range(NUM_BINS), final_distribution)
    plt.title("Averaged Probability Distribution")
    plt.xlabel("Modulo Value")
    plt.ylabel("Probability")
    plt.show()
"""
