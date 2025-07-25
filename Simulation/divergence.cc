/* Verification of Bias, using multi_thread to spead up*/
// Test the divergence of distribution for different plus and minus sign.
#include<cmath>
#include<iostream>
#include<cmath>
#include<set>
#include<vector>
#include<random>
#include<cassert>
#include<iomanip>
#include<thread>
#include<future>
#include<fstream>

using namespace std;

//inline num_gen();

float compute_sei(const int& plus, const int& minus, const float& nbr, const vector<float>& Ptheory){
    int p = 127;
    // vector<float> U(p);
    // for(auto i = 0; i < U.size(); i++){
    //     U[i] = 1.0/127.0;
    // }
    set<int> E = {1, 2, 4, 8, 16, 32, 64};
    set<int> E_prime = {126, 125, 123, 119, 111, 95, 63};
    vector<float> P(p);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(0, E.size() - 1);
    for(auto i = 0; i < nbr; i++){
        int sum{0};
        //sum up
        for(auto j = 0; j < plus; j++){
            auto it = std::next(E.begin(), distrib(gen));
            sum = (sum + *it) % p;
        }
        // minus
        for(auto k = 0; k < minus; k++){
            auto it = std::next(E_prime.begin(), distrib(gen));
            sum = (sum + *it) % p;
        }
        P[sum] += 1.0;
    }
    for(auto& prob: P){
        prob = prob/nbr;
    }
    float sei{0};
    assert(Ptheory.size() == p);
    for(auto i = 0; i < P.size(); i++){
        sei += pow(P[i]-Ptheory[i],2);
    }
    sei *= static_cast<float>(p);
    return sei;
    //return P;

}
int main(){
    float runs = 10;
    int threads = 12;
    float nbr = pow(2,20.156);
    cout << "Nbr of samples" << nbr << endl;
    int plus = 6;
    int minus = 6;
    fstream input("ptheory.txt");
    vector<float> Ptheory;
    if(!input){
        cerr << "Error, can not read file" << endl;
    }
    float prob;
    while(input >> prob){
        Ptheory.emplace_back(prob);
    }
    // for(auto prob: Ptheory){
    //     cout << prob << endl;
    // }
    float sei{0};
    cout << plus << minus << endl;
    for(auto r = 0; r < runs; r++){
        vector<float> results(threads);
        vector<future<float>> futures(threads);
        for(decltype(futures)::size_type i = 0; i< threads; i++){
            futures[i] = std::async(compute_sei, plus, minus, nbr, Ptheory);
        }
        for(decltype(futures)::size_type i = 0; i < threads; i++){
            results[i] = futures[i].get();
        }
        // for(auto sei : results){
        //     cout << log2(sei) << endl;
        // }
        
        for(auto i = 0; i< threads; i++){
            sei += results[i];
        }
    }
    sei /= (threads*runs);
    cout << "SEI: " << setprecision(10) <<  log2(sei) << endl;

    cout << "SEI theory to Uniform random" << endl;
    float seicheck{0};
    for(auto prob: Ptheory){
        seicheck += pow(prob - 1.0/127.0,2);
    }
    seicheck *= 127;
    cout << log2(seicheck);
    // vector<float> distribution(127);
    // for(auto i = 0; i < distribution.size(); i++){
    //     for(auto j = 0; j < threads; j++){
    //         distribution[i] += results[j][i];
    //     }
    //     distribution[i]/= run;
    // }
    // for(auto k = 0; k < distribution.size(); k++){
    //     cout << distribution[k] << endl;
    // }
}