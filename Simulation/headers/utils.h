#ifndef UTILS_H
#define UTILS_H
#include<iostream>
#include<vector>
#include<set>
#include<random>
#include<unordered_map>
#include<cassert>
#include"misc.h"
#include<utility>
#include<algorithm>
#include<chrono>
using namespace std;

using Set = set<int>;
using Vec = vector<int>;
using Pair = pair<Vec, int>;
/* Hashing utilities*/
struct seq_hash{
    std::size_t operator()(const vector<int>& v) const; 
};

using Map = unordered_map< vector<int>, vector<Pair>, seq_hash>;

/*Bernoulli distribution*/
inline bool filter(float p){
    random_device rd;
    mt19937 gen(rd());
    bernoulli_distribution distrib(p);
    return distrib(gen);
}

/*Arimethics, inline*/
inline int modAdd(const int& a, const int& b, const int field_size){
    int res = (a + b) % field_size;
    return (res < 0)? res + field_size : res;
}

inline int modSub(const int& a, const int& b, const int& field_size){
    int res = (a-b) % field_size;
    return (res < 0)? res + field_size : res;
}
int product(const Vec& a, const Vec& b, const int& field_size);

inline Vec addition(const Vec&a, const Vec& b, const int& field_size){
    assert(a.size() == b.size());
    Vec res(a.size());
    for( auto i = 0; i < res.size(); i++){
        //res[i] = (a[i] + b[i]) % field_size;
        res[i] = modAdd(a[i], b[i], field_size);
    }
    return res;
}
inline Vec subtraction(const Vec&a, const Vec&b, const int& field_size){
    Vec res(a.size());
    for(auto i = 0; i < res.size(); i++){
        res[i] = modSub(a[i], b[i], field_size);
    }
    return res;
}

inline Vec scalar_product(const Vec& a, const int& scalar, const int& field_size){
    Vec res(a.size());
    for(auto i = 0; i < res.size(); i++){
        res[i] = (scalar* a[i]) % field_size;
        if(res[i] < 0){
            res[i] = res[i] + field_size;
        }
    }
    return res;
}

/* Functions in the simulations*/

/* Secret generator
@k: length of the secret
@E: the domain of each element in the secret, i.e., the restricted group
*/
vector<int> secret_gen(const int& k, const Set& E);
/*------------------------------------------------*/

/* Parity generator
@k : length of the parity
@field_size: size of the field

Return a parity vector of length k in F_{field_size}^k
*/
Vec parity_gen(const int& k, const int& field_size);
/*------------------------------------------------*/

/* error generator
@E: error domain

Return an arbitrary element in E
*/
int error_gen(const Set& E);
/*------------------------------------------------*/


/*Samples generator
@secret: the secret vector
@k:= secret.size() length of the samples
@E: domain of the errors
@ nbr_samples: number of samples.
Returns a pair of (sample, noisy product)
*/
vector<Pair> sample_gen(const Vec& secret, const Set& E, const int& nbr_samples, const int& field_size);
/*------------------------------------------------*/

/* Sort pair of (sample, noisy product) by delta positions to a hash table
@delta: the number of positions (the index ranging from k-delta+1 to k)
@table: the table
@sample: vector of Pair
inline
Return the table
*/
inline Map sort_sample(Map& table, const int& delta, Pair& sample){
    Vec delta_tag(sample.first.end() - delta, sample.first.end());
    //cout << delta_tag;
    auto where_tag = table.find(delta_tag);
    if(where_tag == table.end()){
        vector<Pair> bucket;
        bucket.emplace_back(sample);
        table.insert(make_pair(delta_tag, bucket));
    }
    else{
        //cout << (*where_tag).first;
        // This is taking a lot of time (not sure why)
        if(find((*where_tag).second.begin(), (*where_tag).second.end(), sample) == (*where_tag).second.end()){
            (*where_tag).second.emplace_back(sample);
        }
        //(*where_tag).second.emplace_back(sample);
    }
    return table;
}
/*------------------------------------------------*/


/*Samples t combine
@:samples_set : set of pairs (sample, noisy product)
@t: the amount of pairs
@E_prime: domain of the coefficients for each samples

*/

void combine_t_samples(const int& t, const int& delta, const vector<Pair>& samples, const vector<Vec>& allCoefficients, Map& table, const int& field_size);
/*------------------------------------------------*/

/* Combine sample in the table to get 2t_samples
@: table

Return the samples in a set

/*------------------------------------------------*/

set<Pair> merge_list(const Map& list1, const Map& list2, const int& field_size, const Set& E);

// Recursive function to generate combinations of coefficients from a set (i.e., |E|^t)
void generate_coefficients(const Set& coefficients, int t, Vec& currentCoefficient,
                           vector<Vec>& allCoefficients);
/*------------------------------------------------*/

// Recursive function to generate all selections of t indices out of n, i.e., n choose t
void generate_combinations(const Vec& indices, int t, Vec& currentCombination, vector<Vec>& allCombinations, int index);

/*Covering codes, inline*/


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
inline Vec codeword_gen(const int& i, const int& length, const int& field_size){
    assert(i < field_size && i >= 0);
    Vec res(length, i);
    return res;
}

/* Find the errors corresponding to the closest codewords from a [1,3]- F_q repetition code
@word: length-3 word (h1, h2, h3)
@max_weight: maximum error weight = 3x3.
@field_size
Return all possible errors .
*/

inline vector<Vec> possible_errors(const Vec& word, const int& field_size){
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

// Decode a random words into blocks of repetition code. Each error block is selected pseudo randomly.
inline Vec decode(const Vec& words, const int field_size){
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
        for(auto value: chosen_vector){
            res.emplace_back(value);
        }        
    }
    return res;
}

inline pair<int,int> decomp(const Vec& word, const int field_size){
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

/* Find the closest codeword from nbr_codes*[1,3] - F_q concatenated repetition codes./*/
// inline Vec get_info(const Vec& words, const int& field_size){
//     assert(words.size() % 3 == 0);
//     int nbr_codes = words.size()/3;
//     Vec info;
//     for(auto i = 0; i< nbr_codes; i++){
//         Vec word(words.begin() + 3*i, words.begin() + 3*(i+1));
//         auto codeword = closest_codeword(word, field_size);
//         info.emplace_back(codeword[0]);
        
//     }
//     return info;
// }
#endif
