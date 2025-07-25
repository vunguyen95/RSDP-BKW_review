#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <thread>
#include <mutex>
#include <fstream>
#include <numeric> // For std::iota
#include <cmath>

const int MODULO = 127;
const int NUM_NUMBERS = pow(2,20); // 2^20 numbers
const int NUM_THREADS = 1;      // Number of threads
//const int NUM_RUNS = 12;         // Number of runs (one per thread)

// Sets E and E_minus
std::vector<int> E = {1, 2, 4, 8, 16, 32, 64};
std::vector<int> E_minus = {63, 95, 111, 119, 123, 125, 126};

// Mutex for thread-safe updates
std::mutex result_mutex;

// Function to perform a single test run
void perform_test(std::array<double, MODULO>& result) {
    // Thread-local random engine
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distrib_E(0, E.size() - 1);
    std::uniform_int_distribution<int> distrib_E_minus(0, E_minus.size() - 1);

    // Local histogram
    std::array<int, MODULO> local_counts = {0};

    // Generate random numbers
    for (int i = 0; i < NUM_NUMBERS; ++i) {
        int sum_E = 0, sum_E_minus = 0;

        // Sum 6 random elements from E
        for (int j = 0; j < 6; ++j) {
            sum_E += E[distrib_E(gen)];
            sum_E_minus += E_minus[distrib_E_minus(gen)];
        }

        // Compute the modulo
        int random_number = (sum_E + sum_E_minus) % MODULO;
        local_counts[random_number]++;
    }

    // Normalize to probability distribution
    std::array<double, MODULO> local_probs;
    for (int i = 0; i < MODULO; ++i) {
        local_probs[i] = static_cast<double>(local_counts[i]) / NUM_NUMBERS;
    }

    // Thread-safe update of the shared result
    std::lock_guard<std::mutex> lock(result_mutex);
    for (int i = 0; i < MODULO; ++i) {
        result[i] += local_probs[i];
    }
}

int main() {
    // Final probability distribution
    std::array<double, MODULO> final_distribution = {0.0};

    // Threads
    std::vector<std::thread> threads;

    // Launch threads
    for (int t = 0; t < NUM_THREADS; ++t) {
        threads.emplace_back(perform_test, std::ref(final_distribution));
    }

    // Wait for all threads to finish
    for (auto& th : threads) {
        th.join();
    }

    // Average the results
    for (double& prob : final_distribution) {
        prob /= NUM_THREADS;
    }

    // Print the final distribution
    for (int i = 0; i < MODULO; ++i) {
        std::cout << "Value " << i << ": " << final_distribution[i] << std::endl;
    }

    std::fstream input("ptheory.txt");
    std::vector<float> Ptheory;
    if(!input){
        std::cerr << "Error, can not read file" << std::endl;
    }
    float prob;
    while(input >> prob){
        Ptheory.emplace_back(prob);
    }

    
    double sei{0};
    for(auto i = 0; i < final_distribution.size();++i){
        sei+= pow(final_distribution[i] - Ptheory[i],2);
    }
    sei *= MODULO;
    std::cout << "SEI to Theory= 2^" << log2(sei) <<  std::endl;
    return 0;
}