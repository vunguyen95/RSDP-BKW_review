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

vector<float> compute_sei(const int& plus, const int& minus, const float& nbr){
    int p = 127;
    //int p = 31;
    vector<float> U(p);
    for(auto i = 0; i < U.size(); i++){
        U[i] = 1.0/127.0;
        //U[i] = 1.0/31.0;
    }
    //set<int> E = {1, 2, 4, 8, 16};
    //set<int> E_prime = {30, 29, 27, 23, 15};
    set<int> E = {1, 2, 4, 8, 16, 32, 64};
    set<int> E_prime = {126, 125, 123, 119, 111, 95, 63 };
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
    /*float sei{0};
    for(auto i = 0; i < P.size(); i++){
        sei += pow(P[i]-U[i],2);
    }
    sei *= static_cast<float>(p);
    return sei;*/
    return P;

}
int main(){
    float run = 12.0;
    int threads = 12;
    float nbr = pow(2,20.156)/12;
    int plus = 6;
    int minus = 6;



    cout << plus << minus << endl;
    vector<vector<float>> results(threads);
    vector<future<vector<float>>> futures(threads);
    for(decltype(futures)::size_type i = 0; i< threads; i++){
        futures[i] = std::async(compute_sei, plus, minus, nbr);
    }
    for(decltype(futures)::size_type i = 0; i < threads; i++){
        results[i] = futures[i].get();
    }
    /*float sei{0};
    for(auto i = 0; i< threads; i++){
        sei += results[i];
    }
    sei /= run;
    cout << "SEI: " << setprecision(10) <<  log2(sei) << endl;*/

    vector<float> distribution(127);
    for(auto i = 0; i < distribution.size(); i++){
        for(auto j = 0; j < threads; j++){
            distribution[i] += results[j][i];
        }
        distribution[i]/= threads;
    }
    // for(auto k = 0; k < distribution.size(); k++){
    //     cout << distribution[k] << endl;
    // }
    fstream input("ptheory.txt");
    vector<float> Ptheory;
    if(!input){
        cerr << "Error, can not read file" << endl;
    }
    float prob;
    while(input >> prob){
        Ptheory.emplace_back(prob);
    }
    float check_sei{0};
    for(auto k = 0; k < distribution.size(); k++){
        check_sei += pow(distribution[k] - Ptheory[k],2);
    }
    check_sei *= 127;
    cout << log2(check_sei);

}

/*int main(){
    int p = 127;
    float run = 10.0;
    set<int> E = {1, 2, 4, 8, 16, 32, 64};
    set<int> E_prime = {126, 125, 123, 119, 111, 95, 63 };
    vector<float> U(127);
    for(auto i = 0; i < U.size(); i++){
        U[i] = 1.0/127.0;
    }

    float nbr = pow(10,7);
    int comb_plus = 2;
    int comb_minus = 6;
    cout << comb_plus << endl;
    cout << comb_minus << endl;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, E.size() - 1);
    float sei_avg{0};
    for(auto count = 0; count < run; count++){
    vector<float> P(127);
    

    for(auto i = 0; i < nbr; i++){
        int sum{0};
        for(auto j = 0; j < comb_plus; j++){
            auto it = std::next(E.begin(), distribution(gen));
            sum = (sum + *it) % p;
            //cout << *it << endl;
            //cout << sum << endl;
        }
        for(auto k = 0; k < comb_minus; k++){
            auto it = std::next(E_prime.begin(), distribution(gen));
            sum = (sum + *it) % p;
        }
        assert(sum < p);
        P[sum] += 1.0;
    }
    for(auto &prob: P){
        prob = prob/nbr;
        //cout << prob << endl;
    }
    cout << "--------------SEI--------------" << endl;
    float sei{0};
    for(auto i = 0; i< P.size(); i++){
        sei += pow(P[i]-U[i],2);
    }
    sei*= p;
    sei_avg += sei;
    
    }
    cout << setprecision(10) << log2(sei_avg/run) << endl;
}*/