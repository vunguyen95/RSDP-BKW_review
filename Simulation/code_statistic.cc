#include<iostream>
#include<vector>
#include<random>
#include<set>
#include"headers/misc.h"
#include<math.h>
#include<cassert>
//#include"headers/utils.h"
#include <iomanip>
#include<map>
#include<fstream>
#include <chrono>
using namespace std;
using Vec = vector<int>;
using Map = map<int, vector<Vec>>;
//For a random words (a,b,c) in F_127^3, find the nearest codeword from 1-3 rep code and the error weight.
// Compute how many words gives weight-[0 - 9] errors.
// If possible how many "nearest codewords" for any words. 
inline int hamming_weight(int n){
    int count = 0;
    while(n!= 0){
        n &= (n -1);
        count++;
    }
    return count;
}
inline int weight_ELEM(const int& n, const int& field_size){
    assert(n >= 0 && n < field_size);
    auto elem_POS = n;
    auto elem_NEG = (-n) + field_size;
    auto weight_POS = hamming_weight(elem_POS);
    auto weight_NEG = hamming_weight(elem_NEG);
    return (weight_POS >= weight_NEG) ? weight_NEG : weight_POS;
    
}

inline int weight_VEC(const Vec& a, const int& field_size){
    int weight = 0;
    for(auto i = 0; i < a.size(); i++){
        weight += weight_ELEM(a[i], field_size);
    }
    return weight;
}

inline int modAdd(const int& a, const int& b, const int field_size){
    int res = (a + b) % field_size;
    return (res < 0)? res + field_size : res;
}
inline Vec addition(const Vec&a, const Vec& b, const int& field_size){
    assert(a.size() == b.size());
    Vec res(a.size());
    for( auto i = 0; i < res.size(); i++){
        //res[i] = (a[i] + b[i]) % field_size;
        res[i] = modAdd(a[i], b[i], field_size);
    }
    return res;
}

Vec parity_gen(const int& k, const int& field_size){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, field_size - 1);
    Vec parity(k);
    for(auto i = 0; i < parity.size(); i++){
        parity[i] = distribution(gen);
    }
    assert(parity.size() == k);
    return parity;
}
inline Vec codeword_gen(const int& i, const int& length, const int& field_size){
    //assert(i < field_size && i >= 0);
    Vec res(length, i);
    return res;
}
int closest_codeword(const Vec& word, const int& field_size, Map& table){
    //assert(word.size() == 3);
    int min_weight = 10;
    Vec correct_codeword(word.size());
    //Vec correct_error;
    Vec error(word.size());
    for(auto i = 0; i < field_size; i++){
        //auto random_codeword = codeword_gen(i, word.size(), field_size);
        
        for(int j = 0; j < error.size(); j++){
            auto temp = word[j] - i;
            error[j] = (temp >= 0) ? temp : (temp + field_size) % field_size;
        }
        auto error_weight = weight_VEC(error, field_size);

        auto where = table.find(error_weight);
        Vec codeword(3,i);
        if(where == table.end()){
            vector<Vec> bucket;
            bucket.emplace_back(codeword);
            table.insert(make_pair(error_weight, bucket));
        }
        else{
            (*where).second.emplace_back(codeword);
        }


        if(error_weight < min_weight){
            min_weight = error_weight;
            //correct_codeword = codeword_gen(i, word.size(), field_size);
            //correct_error = error;
        }     
    }
    //cout << correct_error;
    //cout << weight_VEC(correct_error, field_size) << endl;
    return min_weight;
}

//new, take a random error instead of the first one.
vector<Vec> possible_errors(const Vec& word, const int& field_size){
    vector<Vec> res;
    vector<vector<Vec>> temp(10);
    //assert(word.size() == 3);
    int min_weight = 10;
    Vec correct_codeword(word.size());
    //Vec correct_error;
    Vec error(word.size());
    for(auto i = 0; i < field_size; i++){
        
        for(int j = 0; j < error.size(); j++){
            auto temp = word[j] - i;
            error[j] = (temp >= 0) ? temp : (temp + field_size) % field_size;
        }
        auto error_weight = weight_VEC(error, field_size);
        temp[error_weight].emplace_back(error);
    }
    for(auto i = 0; i < temp.size(); i++){
            if(temp[i].size()!= 0){
                for(auto error: temp[i]){
                    res.emplace_back(error);
                }
                break;
            }
        }
    return res;
}

pair<int,int> decomp(const Vec& word, const int field_size){
    int min_weight = 10;
    Map table;
    pair<int,int> res = {0,0};
    //Vec correct_codeword;
    //Vec correct_error;
    Vec error(word.size());
    for(auto i = 0; i < field_size; i++){
        for(int j = 0; j < error.size(); j++){
            auto temp = word[j] - i;
            error[j] = (temp>=0) ? temp : (temp+field_size) % field_size;
        }
        auto error_weight = weight_VEC(error,field_size);
        if(error_weight <= min_weight){
            min_weight = error_weight;
            //correct_codeword = codeword_gen(i, word.size(), field_size);
            //correct_error = error;

            auto where = table.find(error_weight);
            if(where == table.end()){
                vector<Vec> bucket;
                bucket.emplace_back(error);
                table.insert(make_pair(error_weight,bucket));
            }
            else{
                (*where).second.emplace_back(error);
            }
        }
    }
    //cout << correct_codeword;
    //cout << correct_error;
    auto where = table.find(min_weight);
    for(auto e : (*where).second){
        for(auto i = 0; i < e.size(); i++){
        auto weight = weight_ELEM(e[i], field_size);
        if(weight == hamming_weight(e[i])){
            res.first+= weight;
        }
        else{
            res.second += weight;
        }
    }
    }
    
    return res;
}

pair<int,int> decomp1(const Vec& word, const int field_size){
    pair<int, int> res;
    for(auto i = 0; i < word.size(); i++){
        auto weight = weight_ELEM(word[i], field_size);
        if(weight == hamming_weight(word[i])){
            res.first+= weight;
        }
        else{
            res.second += weight;
        }
    }
    return res;
}

Vec decode(const Vec& words, const int field_size){
    Vec res;
    assert(words.size() % 3 == 0);
    int nbr_codes = words.size()/3;
    Vec info;
    for(auto i = 0; i< nbr_codes; i++){
        Vec word(words.begin() + 3*i, words.begin() + 3*(i+1));
        auto errors = possible_errors(word,field_size);
        //Randomized the selected error
        auto now = std::chrono::system_clock::now();
        auto duration = now.time_since_epoch();
        size_t millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        size_t index = millis % errors.size();
        auto chosen_vector = errors[index];
        for(auto entry: chosen_vector){
            res.emplace_back(entry);
        }        
    }
    return res;
}

int main(){
    int p = 127;
    pair<int,int> probs = {0,0};
    /*for(auto first = 0; first < p ; first++){
        for(auto second = 0; second < p; second++){
            for(auto third = 0; third < p; third++){
                Vec random = {first,second, third};
                auto dec = decomp(random, p);
                probs.first += dec.first;
                probs.second += dec.second;
            }
        }
    }
    cout << probs.first << endl;
    cout << probs.second << endl;
   */ 
   Vec test = {64,25,63};
   auto errors = possible_errors(test, p);
   for(auto error: errors){
        cout << error;
        auto dec = decomp1(error, p);
        cout << dec.first << " " << dec.second << endl;
   }
   Vec test0 = {20, 102 ,96};
    errors = possible_errors(test0, p);
   for(auto error: errors){
        cout << error;
        auto dec = decomp1(error, p);
        cout << dec.first << " " << dec.second << endl;
   }
//    Vec test1 = {64, 25, 63, 20, 102 ,96};
//    cout << "-----Chosen-----------" << endl;
//    auto random = decode(test1, p);
//    cout << random;
//    auto dec = decomp1(random, p);
//    cout << dec.first << " " << dec.second;

// Test the probability of kept codeword when the length increase.
   int trials = pow(2,15);
   int count{0};
   for(auto i = 0; i < trials; i++){
        auto randomshit = parity_gen(12, p);
        auto error = decode(randomshit, p);
        if(weight_VEC(error,p) == 16){
            auto dec = decomp1(error, p);
            if(dec.first == dec.second){
                count++;
                // cout << error;
                // cout << dec.first << dec.second << endl;
            }
        }
        
        
   }
   cout << count << endl;

} 



/*int main(){
    //int p = 31;
    int p = 127;
    //string fileName = "code_distrib.txt";
    //ofstream outputFile(fileName);
    vector<float> error_distrib(10);
    long long int avg_error_weight = 0;
    for(auto first = 0; first < p; first ++){
        for(auto second = 0; second < p; second++){
            for(auto third = 0; third < p; third++){
                Map table;
                Vec random = {first, second,third};
                auto min_weight = closest_codeword(random, p, table);
                error_distrib[min_weight] += 1.0;
                avg_error_weight += min_weight;
                //outputFile << "Words:" << random;
                //outputFile << "Smallest distance: " << min_weight;
                //auto lookup = *(table.find(min_weight));
                //outputFile << "Number of nearest codewords: " << lookup.second.size() << endl;
                //for (auto codeword: lookup.second){
                //    outputFile << codeword;
                //}
                //outputFile << "-----------------------------------" << endl;

            }
        }
    }
    for(auto& prob: error_distrib){
        prob /= static_cast<float>(pow(p,3));
    }
    //cout << max_word;
    auto avg = static_cast<float>(avg_error_weight)/static_cast<float>(pow(p,3));
    float exp_weight{0};
    cout << "Average error weight: " << avg << endl;
    cout << "Error Distribution: " << endl;
    for(auto i = 0; i < error_distrib.size(); i++){
        cout << "Weight: " << i << ": " << error_distrib[i] << "; ";
        exp_weight += i*error_distrib[i];
    }
    
    
    cout << "Expected weight: " << exp_weight << endl;
    
   vector<int> counts(3);
   for(auto i = 1; i < p; i++){
    auto weight = weight_ELEM(i,p);
    counts[weight-1]++;
   }
   cout << counts;
}*/
/*

    Vec test = {2, 2, 4};
    Map table_test;
    cout << closest_codeword(test, 127, table_test);
    auto min_weight = closest_codeword(test, 127, table_test);
    outputFile << "Words: " << test;
    outputFile << "Smallest distance" << min_weight << endl;
    outputFile << "Number of nearest codewords: " << (*table_test.find(min_weight)).second.size() << endl;
    for(auto codeword: (*table_test.find(min_weight)).second){
        outputFile << codeword;
    }
    for(auto entry: table_test){
        cout << "Distance: " << entry.first << endl;
        cout << "Codewords : " << endl;
        for(auto& codeword: entry.second){
            cout << codeword;
        }
    }
    
    
}*/