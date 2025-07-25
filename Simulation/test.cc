#include<iostream>
#include<vector>
#include<random>
#include<set>
#include"headers/misc.h"
#include<math.h>
#include"headers/utils.h"
#include <iomanip>
using namespace std;


int main(){
    /*Params*/
    Set E = {1, 2, 4, 8, 16,32,64}; 
    
    Vec currentCoefficient;
    vector<Vec> allCoefficients;
    int t = 2;
    int nbr_samples = 15;
    int field_size = 127;
    int length = 3;
    int steps = 2;
    //int delta1=2;
    //int delta2=1;
    //int delta = delta1 + delta2;
    int delta = 1;
    
    // Initialization
    generate_coefficients(E, t, currentCoefficient, allCoefficients);
    //cout << allCoefficients;
    cout << allCoefficients.size() << endl;
    Vec secret = {0,0,0};
    cout << "Secret: " << secret;
    cout << "-------------------------------------" << endl;
    //Test again if the noisy products are combined correctly
    vector<Pair> samples = {{{1,1,1}, 1}, {{2,1,1}, 2}};
    Map test;
    Map test2;
    combine_t_samples(t, delta, samples, allCoefficients, test, field_size);
    combine_t_samples(t,delta, samples, allCoefficients, test2, field_size);
    auto test3 = merge_list(test2,test,field_size);
}


//Test again if the noisy products in merge list are combine correctly.

/*--------------------TEST SPACE----------------------------------------------------------------*/

 // Test scalar product
    /*for(auto scalar : E){
        cout << scalar_product(secret, scalar, 127);
    }*/

/* 
    //Test generation (OK.)
    auto samples = sample_gen(secret, E, 1000, 127);
    for(auto i = 0; i < samples.size(); i++){
        auto parity = samples[i].first;
        auto noisy_product = samples[i].second;
        auto error = (noisy_product - product(parity, secret, 127)) % 127;
        if (error < 0){
            error = error + 127;
        }
        assert(E.find(error) != E.end());
    }

    */
/*
 //Test linear combinations
    Vec currentCombination;
    vector<Vec> allCombinations;
    Vec indices = {0, 1,2,3,4,5,6,7,8,9};
    generate_combinations(indices, t, currentCombination, allCombinations, 0);
    for(auto combinations: allCombinations){
        cout << combinations;
    }
    cout << allCombinations.size() << endl;
    */

/*
    // Test sort samples
    Map table;
    auto samples = sample_gen(secret, E, 5, 127);
    for(auto sample: samples){
        cout << sample.first;
        cout << sample.second << endl;
        sort_sample(table ,3, sample);
    }
    for(auto it = table.begin(); it!= table.end(); it++){
        cout << "Tag: " << (*it).first;
        for(auto j : (*it).second){
            cout << "sample: " << j.first;
            cout << "noisy product: " << j.second << endl;
            auto error = (j.second - product(j.first, secret, 127)) % 127;
            if (error < 0){
                error = error + 127;
            }
            cout << "Check: " << error << endl;
        }
    }
    */


/* //TEST weight function
    for(auto i = 0; i < 127; i++){
        //cout << codeword_gen(i, 3, 127);
    }
    cout << hamming_weight(127) << endl;
    int weight_max = 1;
    for(auto i = 0; i <= 126; i++){
        if(weight_ELEM(i, 127) > weight_max){
            weight_max = weight_ELEM(i, 127);
        }
    }
    cout << weight_max << endl;

    vector<int> counts = {0, 0, 0};
    for(int i = 1; i<= 126; i++){
        auto pos = weight_ELEM(i, 127) - 1;
        counts[pos]++;
    }
    cout << counts;
   cout << weight_VEC({17, 80, 32}, 127);
    */


// Test covering codes
    /*Vec test = {126, 105, 84};  
    cout << closest_codeword(test, 127);   

   long long int max_error_weight = 0;
    for(auto first = 0; first < 127; first ++){
        for(auto second = 0; second < 127; second++){
            for(auto third = 0; third < 127; third++){
                Vec test1 = {first, second, third};
                Vec codeword = closest_codeword(test1, 127);
                Vec error(test1.size());
                for(auto i = 0; i < error.size(); i++){
                    auto temp = test1[i] - codeword[i];
                    error[i] = (temp >= 0) ? temp : (temp + 127) % 127;
                }
                max_error_weight += weight_VEC(error, 127);
            }
        }
    }
    //cout << max_word;
    cout << max_error_weight;
    auto avg = static_cast<float>(max_error_weight)/static_cast<float>(pow(127,3));
    cout << avg;
    */


  /*
    // Test combine t samples
    cout << "Test combine:" << endl;
    auto samples = sample_gen(secret,E, 4, 127);
    for(auto sample: samples){
        cout << sample.first;
        cout << sample.second << endl;
    }
    Map table;
    combine_t_samples(t, delta, samples, allCoefficients, table, 127);
    cout << "---------------table----------------------: " << endl;
    for(auto it = table.begin(); it!= table.end(); it++){
        cout << "Tag: " << (*it).first;
        for(auto j : (*it).second){
            cout << "sample: " << j.first;
            cout << "noisy product: " << j.second << endl;
            auto error = (j.second - product(j.first, secret, 127)) % 127;
            if (error < 0){
                error = error + 127;
            }
            cout << "Check: " << error << endl;
        }
    }
*/

/*    
  //Test 2t_samples BRUTE_FORCE
    cout << "Brute-force 2t samples to compare: " << endl;
    set<Pair> dft_sample_test;
    for(auto i = 0; i < t_sample_test.size() -1; i++){
        auto minus = scalar_product(t_sample_test[i].first, -1, field_size);
        auto minus_checksum =(-1) * t_sample_test[i].second;
        for(auto j = i + 1; j < t_sample_test.size(); j++){
            // a- b
            auto dft_vector = addition(minus, t_sample_test[j].first, field_size);
            Vec delta_tag(dft_vector.end()-delta, dft_vector.end());
            Vec zeros(delta);

            if(delta_tag == zeros){
                auto dft_checksum = (minus_checksum + t_sample_test[j].second) % field_size;
                if(dft_checksum < 0){
                    dft_checksum += field_size;
                }
                auto dft_sample = make_pair(dft_vector, dft_checksum);
                dft_sample_test.insert(dft_sample);
            }
            // a+ b
            auto dft_vector_sum = addition(t_sample_test[i].first, t_sample_test[j].first, field_size);
            Vec delta_tag_sum(dft_vector_sum.end()-delta, dft_vector_sum.end());
            if(delta_tag_sum == zeros){
                auto dft_checksum_sum = (t_sample_test[i].second + t_sample_test[j].second)% field_size;
                if(dft_checksum_sum < 0){
                    dft_checksum_sum += field_size;
                }
                auto dft_sample_sum = make_pair(dft_vector_sum, dft_checksum_sum);
                dft_sample_test.insert(dft_sample_sum);
            }
            

        }
    }
    cout << dft_sample_test.size() << endl;  

*/