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
    Set E_prime = {-1, -2, -4, -8, -16, -32, -64, 1, 2, 4, 8, 16, 32, 64};
    //Set E_prime = {1,2,4,8,16};
    //Set E = {1, 2, 4, 8, 16};
    //Set E_prime = {-1, -2, -4, -8, -16, 1, 2, 4, 8, 16};
    //Set E_prime = {0, 1, 2, 4, 8, 16, 30, 29, 27, 23, 15};
    
    Vec currentCoefficient;
    vector<Vec> allCoefficients;
    int t = 2;
    int nbr_samples = 35;
    cout << "NBR SAMPLES:" << nbr_samples << endl;
    int field_size = 127;
    int length = 9;
    int steps = 2;
    int delta1=2;
    int delta2=1;
    int delta = delta1 + delta2;
    assert((length - delta) % 3 == 0);
    // Initialization
    generate_coefficients(E, t, currentCoefficient, allCoefficients);
    //cout << allCoefficients;
    cout << allCoefficients.size() << endl;
    Vec secret = {1,4,16,8,32, 2,64, 1, 16};
    cout << "Secret: " << secret;
    cout << "-------------------------------------" << endl;
//----------------------------------------------------------------------
    //Inital lists, 4 batches of M samples. From each batch, create a list at depth 0 in the tree

    //Each list in the tree is sorted according to deltai (except for the last list). At depth 0, return 4 hash table.
    vector<Map> listsDepth0(4);
    cout << listsDepth0.size() << endl;
    for(int i = 0; i< 2; i++){
        //Create two batches of M samples
        auto list1 = sample_gen(secret, E, nbr_samples, field_size);
        auto list2 = sample_gen(secret, E, nbr_samples, field_size);
        //From each batches, make the list in Depth 0.
        combine_t_samples(t, delta1, list1, allCoefficients, listsDepth0[2*i], field_size);
        combine_t_samples(t, delta1, list2, allCoefficients, listsDepth0[2*i+1], field_size);
        
    }
    cout << "Checking depth 0: " << endl;
    for(auto list: listsDepth0){
        int count{0};
        for(auto entries: list){
            count+= entries.second.size();
        }
        cout << count << endl;
    }
//----------------------------------------DEPTH 1--------------------------

// In depth 1, merge again two obtain the final list in depth 2
    vector<Map> listsDepth1(2);
    for(int i = 0; i < 2; i++){
        //Store unique combinations from depth 0 lists
        auto delta1_samples = merge_list(listsDepth0[2*i], listsDepth0[2*i + 1], field_size, E);
        cout << "Checking depth 1: " << endl;
        cout << "Depth1 size: " << delta1_samples.size() << endl;
        //auto it = delta1_samples.begin();
        // for(auto j = 0; j < 10; j++){
        //     cout << (*it).first;
        //     it++;
        // }
        cout << "Merge completed" << endl;
        for(auto sample : delta1_samples){
            //Populate the lists in depth 1
            sort_sample(listsDepth1[i], delta, sample);
        }
        cout << "Depth 1 populated" << endl;
    }
    
    
    // for(auto list: listsDepth1){
    //     int count{0};
    //     for(auto entries: list){
    //         count+= entries.second.size();
    //         for(auto sample: entries.second){
    //             auto parity = sample.first;
    //             Vec delta1_tag(delta1);
    //             Vec parity_tag(parity.end() - delta1, parity.end());
    //             assert(delta1_tag == parity_tag);
    //         }
    //     }
    //     cout << count << endl;
    // }
//-----------------------------------------Final List----------------------------

    auto final_samples = merge_list(listsDepth1[0], listsDepth1[1], field_size, E);
    cout << "Checking final list size: " << endl;
    cout << final_samples.size() << endl;
    auto it1 = final_samples.begin();
    for(auto j = 0; j < 10; j++){
        cout << (*it1).first;
        it1++;
    }


//----------------------------------------Covering------------------

    set<Pair> Lsht;
    //final_sample = {samples = {parity, bit}}
    //Lsht = {{codewords, bit}} 
    auto word_length = secret.size() - delta1 - delta2;
    
    for (auto sample: final_samples){
        auto parity = sample.first;
        auto sum = sample.second;
        Vec word(parity.begin(), parity.begin() + word_length);
        auto error = decode(word, field_size);
        //only keep the error with desired pattern and weight.
        if(weight_VEC(error, field_size) == 8){
            auto dec = decomp(error, field_size);
            if(dec.first == dec.second){
                auto codeword = subtraction(word, error, field_size);
                //Find the info word, [3,1] repetition codes
                auto nbr_codes = codeword.size()/3;
                Vec info;
                for(auto j = 0; j< nbr_codes;j++){
                    info.emplace_back(codeword[3*j]);
                }

                Lsht.insert(make_pair(info, sum));
            }
        }

        
    }
    cout << "Sample For guessing: " << Lsht.size() << endl;
    // auto it2 = Lsht.begin();
    // for(auto j = 0; j < 100; j++){
    //     cout << (*it2).first;
    //     // auto codeword = (*it2).first;
    //     // auto nbr_codes = codeword.size()/3;
    //     //         Vec info;
    //     //         for(auto j = 0; j< nbr_codes;j++){
    //     //             info.emplace_back(codeword[3*j]);
    //     //         }
    //     // cout << info;
    //     it2++;
    // }
//-----------------------------------------GUESSING---------------
    Vec correct_guess = {secret[0]+secret[1] + secret[2], secret[3]+secret[4]+secret[5]};
    Vec random_guess = {0,0};
    cout << "Correct guess " << correct_guess;
    vector<float> probabilities(field_size);
    for(auto sample: Lsht){
        auto info = sample.first;
        //assert(info.size() == correct_guess.size());
        //auto product = (info[0]*correct_guess[0] + info[1]*correct_guess[1]) % field_size;
        auto product = (info[0]*random_guess[0] + info[1]*random_guess[1]) % field_size;
        
        auto noise = modSub(sample.second, product, field_size);
        probabilities[noise] = probabilities[noise] + 1.0;
    }
    for(auto prob: probabilities){
        prob /= static_cast<float>(Lsht.size());
        cout << prob << endl;
    }
    

}



