# RSDP-BKW_review
Repository for the reviewing process of the paper "A BKW solver for the Restricted Syndrome Decoding Problem"

./Simulation (C++)
Simulations for small parameters.

g++ -o main main_d1 utils.cc misc.cc 

Some tests:
- bias.cc (testing the change Square Euclidean Distance when varying minus and plus, -pthread needed)
- code_statistic.cc (average error weight for covering codes, probability of keeping desirable codewords)
- divergence.cc (testing the divergence with the theoretical distribution, when generating random samples, -pthread needed)


./complexity_estimation (Python)
- depth2.py/depth3.py/biggerE.py (parameters optimization and complexity for depth_2, depth_3, and larger restricted set)
- convo.py (convolution of probability mass functions)
- plot.py (generating plot)

python3 [...].py
