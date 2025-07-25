#include<iostream>
#include<vector>
#include<random>
#include<set>
#include"headers/misc.h"
#include<math.h>
#include <iomanip>
using namespace std;

using Set = set<int>;
using Vec = vector<int>;

//GLOBAL restricted set E = {2^i}
const int prime = 127;
const int NUM_BIT = 7;
const Vec E_POS = {1, 2, 4, 8, 16, 32, 64};
const Vec E_NEG = {-1, -2, -4, -8, -16, -32 , -64};

const Vec E_prime = {-1, -2, -4, -8, -16, -32 , -64, 1, 2, 4, 8, 16, 32, 64};


Vec bin_Rep(const int& num, const int& NUM_BIT){
    Vec res(NUM_BIT);
    int num_copy = num;
    int i = 0;
    while(num_copy > 0){
        res[i] = num_copy % 2;
        num_copy /= 2;
        i++;
    }
    return res;
}

int weight_ELEM(const int& num, const int& NUM_BIT){
    if(num == 0){
        return 0;
    }
    int num_POS = num;
    int num_NEG = (-num) + prime;
    Vec bin_POS = bin_Rep(num_POS, NUM_BIT);
    Vec bin_NEG = bin_Rep(num_NEG, NUM_BIT);
    int weight_POS{0}, weight_NEG{0};
    for(auto i = 0; i < NUM_BIT; i++){
        weight_POS += (bin_POS[i] == 1)? 1 : 0;
        weight_NEG += (bin_NEG[i] == 1)? 1 : 0;
    }

    return (weight_POS >= weight_NEG)? weight_NEG : weight_POS;
}
int weight_VEC(const Vec& vec, const int& NUM_BIT){
    int weight{0};
    for(auto i = 0; i < vec.size(); i++){
        weight += weight_ELEM(vec[i], NUM_BIT);
    }
    return weight;
}

Vec vec_gen(const int& length, const Vec& E){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distrib(0,E.size()-1);
    Vec res(length);
    for(auto i = 0; i < length; i++){
        auto val = distrib(gen);
        res[i] = E[val];
    }
    return res;
}
void generateVectors(int n, int w, std::vector<int>& currentVector, unsigned long long int& count, const vector<double>& probabilities, double& total_prob) {
    if (n == 0) {
        if (w == 0) {
            // count the probability for each arrangement
            double prob = 1.0;
            for (int element : currentVector) {
                prob *= probabilities[element];
                //std::cout << element << ' ';
            }
            count++;
            total_prob += prob;
            //cout << prob << endl;
        }
        return;
    }

    for (int i = 0; i <= 3; ++i) {
        if (w - i >= 0) {
            currentVector.push_back(i);
            generateVectors(n - 1, w - i, currentVector, count, probabilities, total_prob);
            currentVector.pop_back();
        }
    }
    return;
}

int main(){
    cout << bin_Rep(8, NUM_BIT); 
    cout << weight_ELEM(125, NUM_BIT) << endl;
    cout << weight_VEC({2,127}, NUM_BIT) << endl;

    auto test = vec_gen(10, E_POS);
    cout << test;
    cout << weight_VEC(test, NUM_BIT);
    int weight_max = 1;
    int num = 0;
    for(int i = 1; i <= 126; i++){
        if(weight_ELEM(i, NUM_BIT) > weight_max){
            weight_max = weight_ELEM(i, NUM_BIT);
            num = i;
        }
    }
    cout << "Max weight:" << weight_max << endl;
    cout << "Number: " << num << endl;
    vector<int> counts = {0, 0, 0};
    for(int i = 1; i<= 126; i++){
        auto pos = weight_ELEM(i, NUM_BIT) - 1;
        counts[pos]++;
    }
    cout << counts;
    int n = 10;
    int w = 15;
    unsigned long long int count = 0;
    double total_prob = 0.0;
    //vector<float> probabilities = {1.00000/127.00000 ,14.00000/127.00000, 42.00000/127.00000, 70.00000/127.000000};
    vector<double> probabilities = {0.1, 0.2,0.3, 0.4};
    //Vec currentVector;
    //generateVectors(n,w, currentVector, count, probabilities, total_prob);

    //cout << count << endl;
    //cout << setprecision(15) <<log2(total_prob) << endl;
    
    cout << "------------------" << endl;
    /*long long int max_error_weight = 0;
    Vec max_word;
    for(auto first = 0; first < 127; first ++){
        for(auto second = 0; second < 127; second++){
            for(auto third = 0; third < 127; third++){
                Vec test1 = {first, second, third};
                int min_weight = 10;
                Vec closet_codeword(test1.size());
                Vec correct_error(test1.size());
                for(int i = 0; i < 127; i++){
                    Vec codeword = {i, i, i};
                    Vec error(codeword.size());
                    for(int j = 0; j< error.size(); j++){
                            auto temp = test1[j] - codeword[j];
                            error[j] = (temp >= 0) ? temp : (temp +127)%127;   
                    }
                    //cout << error;
                    //cout << weight_VEC(error,NUM_BIT) << endl;
                    if(weight_VEC(error,NUM_BIT) < min_weight){
                        correct_error = error;
                        min_weight = weight_VEC(error, NUM_BIT);
                    } 
                }
                max_error_weight += min_weight;
            }
        }
    }
    //cout << max_word;
    cout << max_error_weight;*/
    Vec test1 = {126, 105, 84};
    int min_weight = 10;
    Vec closet_codeword(test1.size());
    Vec correct_error(test1.size());
    for(int i = 0; i < 127; i++){
        Vec codeword = {i, i, i};
        Vec error(codeword.size());
        for(int j = 0; j< error.size(); j++){
                auto temp = test1[j] - codeword[j];
                error[j] = (temp >= 0) ? temp : (temp +127)%127;   
        }
        cout << error;
        cout << weight_VEC(error,NUM_BIT) << endl;
        if(weight_VEC(error,NUM_BIT) < min_weight){
            correct_error = error;
            min_weight = weight_VEC(error, NUM_BIT);
            closet_codeword = codeword;
        } 
    }
    cout << "--------------: " << endl;
    cout << closet_codeword;
    cout << correct_error;
    cout << min_weight;
}