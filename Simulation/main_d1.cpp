#include<iostream>
#include<vector>
#include<random>
#include<set>
#include"headers/misc.h"
#include<math.h>
#include"headers/utils.h"
#include <iomanip>
#include<cmath>
#include<fstream>
#include<string>
using namespace std;


int main(){
    /*Params*/
    Set E = {1, 2, 4, 8, 16,32,64};
    Set E_prime = {-1, -2, -4, -8, -16, -32, -64, 1, 2, 4, 8, 16, 32, 64};
    //Set E_prime = {1,2,4,8,16};
    //Set E = {1, 2, 4, 8, 16};
    
    Vec currentCoefficient;
    vector<Vec> allCoefficients;
    int t = 2;
    int nbr_samples = 55;
    cout << "NBR SAMPLES:" << nbr_samples << endl;
    int field_size = 127;
    int length = 7;
    int steps = 1;
    int delta = 1;
    assert((length - delta) % 3 == 0);
    // Initialization
    generate_coefficients(E, t, currentCoefficient, allCoefficients);
    //cout << allCoefficients;
    cout << allCoefficients.size() << endl;
    Vec secret = {1,4,16,8,32, 2,64};
    cout << "Secret: " << secret;
    cout << "-------------------------------------" << endl;
//----------------------------------------------------------------------
    //Inital lists, 2 batches of M samples. From each batch, create a list at depth 0 in the tree

    //Each list in the tree is sorted according to deltai (except for the last list). At depth 0, return 4 hash table.
    vector<Map> listsDepth0(2);
    cout << listsDepth0.size() << endl;
    for(int i = 0; i< 1; i++){
        //Create two batches of M samples
        auto list1 = sample_gen(secret, E, nbr_samples, field_size);
        auto list2 = sample_gen(secret, E, nbr_samples, field_size);
        //From each batches, make the list in Depth 0.
        combine_t_samples(t, delta, list1, allCoefficients, listsDepth0[2*i], field_size);
        combine_t_samples(t, delta, list2, allCoefficients, listsDepth0[2*i+1], field_size);
        
    }
    cout << "Checking depth 0: " << endl;
    for(auto list: listsDepth0){
        int count{0};
        for(auto entries: list){
            count+= entries.second.size();
        }
        cout << count << endl;
    }

//-----------------------------------------Final List----------------------------

    auto final_samples = merge_list(listsDepth0[0], listsDepth0[1], field_size, E);
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
    auto word_length = secret.size() - delta;
    
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
    fstream input("ptheory.txt");
    vector<float> ptheory;
    if(!input){
        cerr << "Error, can not read file" << endl;
    }
    float prob;
    while(input >> prob){
        ptheory.emplace_back(prob);
    }
    assert(ptheory.size() == field_size);
    // for(auto prob : ptheory){
    //     cout << prob << endl;
    // }


    vector<float> U(field_size);
    for(auto i = 0; i < U.size(); i++){
        U[i] = 1.0/static_cast<float>(field_size);
    }
    cout << "CORRECT GUESS STATISTIC: " << endl;
    Vec correct_guess = {secret[0]+secret[1] + secret[2], secret[3]+secret[4]+secret[5]};
    cout << "Correct guess " << correct_guess;
    vector<float> probabilities(field_size);
    for(auto sample: Lsht){
        auto info = sample.first;
        //assert(info.size() == correct_guess.size());
        auto product = (info[0]*correct_guess[0] + info[1]*correct_guess[1]) % field_size;
        //auto product = (info[0]*random_guess[0] + info[1]*random_guess[1]) % field_size;
        
        auto noise = modSub(sample.second, product, field_size);
        probabilities[noise] = probabilities[noise] + 1.0;
    }
    for(auto& prob: probabilities){
        prob /= static_cast<float>(Lsht.size());
        // cout << prob << endl;
    }
    float sei{0};
    for(auto i=0; i< field_size; i++){
        sei += pow(probabilities[i] - U[i],2);
    }
    sei*= field_size;
    

    float sei2theory{0};
    for(auto j = 0; j < field_size;j++){
        sei2theory += pow(probabilities[j] - ptheory[j],2);
    }
    sei2theory *= field_size;

    

//-----------------------------------------------------------------
    cout << "Random Guess statistic: " << endl;
    float sei_avgwrong_uniform{0};
    float sei_avgwrong_theory{0};
    for(auto first = 0; first < field_size; first++){
        for(auto second = 0; second < field_size; second++){
            Vec random_guess = {first, second}; //cout << "Random guess: " << random_guess ;
            if(random_guess == correct_guess){
                continue;
            }

            vector<float> P_wrong(field_size);

            for(auto sample: Lsht){
                auto info = sample.first;
                auto product = (info[0]*random_guess[0] + info[1]*random_guess[1]) % field_size;
                auto noise = modSub(sample.second,product, field_size);
                P_wrong[noise] = P_wrong[noise] + 1.0;
            }
            
            for(auto& prob: P_wrong){
                prob /= static_cast<float>(Lsht.size());
                // cout << prob << endl;
            }
            float sei_wrong{0};
            float sei_wrong_theory{0};
            for(auto i = 0; i< field_size; i++){
                sei_wrong += pow(P_wrong[i]- U[i],2);
                sei_wrong_theory += pow(P_wrong[i] - ptheory[i],2);
            }
            sei_wrong *= field_size;
            sei_wrong_theory *= field_size;
            sei_avgwrong_uniform += sei_wrong;
            sei_avgwrong_theory += sei_wrong_theory;


            //cout << "SEI WRONG guess to Puniform = 2^" << setprecision(10) << log2(sei_wrong) << endl;
            if(sei_wrong >= sei){
                cout << "----------------" << endl;
                cout << "Wrong guess is farther from Uniform than Correct guess" << endl;
                cout << random_guess;
                cout << log2(sei_wrong) << endl;
                cout << "--------------" << endl;
                
            }
            // if(sei_wrong_theory <= sei2theory){
            //     cout << "Wrong guess is closer to Theory than Correct guess" << endl;
            //     cout << random_guess;
            //     cout << log2(sei_wrong_theory) << endl;
            //     cout << "---------------" << endl;
            // }
        }

    }
    cout << "SEI CORRECT GUESS to PUniform: = 2^" <<  setprecision(10) << log2(sei) << endl;
    cout << "SEI CORRECT GUESS to PTheory = 2^" << setprecision(10) << log2(sei2theory) << endl;
    cout << "Average SEI of WRONG GUESS to PUniform = 2^" << setprecision(10) << log2(sei_avgwrong_uniform/(pow(127,2))) << endl;
    cout << "Average SEI of WRONG GUESS to PTheory = 2^" << setprecision(10) << log2(sei_avgwrong_theory/(pow(127,2))) << endl;
    cout << "NBR SAMPLES:" << nbr_samples << endl;
}



