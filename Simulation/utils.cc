#include"headers/utils.h"
#include"headers/misc.h"
/*------------------------------------------------*/
std::size_t seq_hash::operator()(const vector<int>& v) const {
        std::size_t seed = v.size();
        for (const int& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

/*------------------------------------------------*/


vector<int> secret_gen(const int& k, const Set& E){
    Vec secret(k);
    if(E.empty() || E.size() == 0){
        return {};
    }
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, E.size() - 1);
    for(auto i = 0; i < secret.size(); i++){
        auto it = std::next(E.begin(), distribution(gen));
        secret[i] = (*it);
    }
    return secret;
}

/*------------------------------------------------*/
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
/*------------------------------------------------*/

int error_gen(const Set& E){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, E.size() - 1);
    auto it = std::next(E.begin(), distribution(gen));
    return (*it);

}
/*------------------------------------------------*/



vector<Pair> sample_gen(const Vec& secret, const Set& E, const int& nbr_samples, const int& field_size){
    vector<Pair> samples(nbr_samples);
    for(auto i  = 0; i < samples.size(); i++){
        auto parity = parity_gen(secret.size(), field_size);
        auto error = error_gen(E);
        auto noisy_product = (product(parity, secret, field_size) + error) % field_size;
        samples[i] = pair<Vec, int>(parity, noisy_product);
    }

    return samples;
}
/*------------------------------------------------*/



/*------------------------------------------------*/

//combine t samples without filter, no dependencies.

void combine_t_samples(const int& t, const int& delta, const vector<Pair>& samples, const vector<Vec>& allCoefficients, Map& table, const int& field_size){
    Vec currentCombination;
    vector<Vec> allCombinations;
    Vec indices(samples.size());
    iota(indices.begin(), indices.end(), 0);
    generate_combinations(indices, t, currentCombination, allCombinations, 0);
    int count{0};
    // Getting all linear combinations with all possible coefficients
    for(auto combination : allCombinations){
        //cout << "Combinations: " << combination;
        for(auto coefficient : allCoefficients){
            //cout << "Coefficients: " << coefficient;
            Vec sumParity(samples[0].first.size());
            int sumCheck = 0;
            assert(combination.size() == coefficient.size());
            for(auto i = 0; i < coefficient.size(); i++){
                //Get the samples according to the index in combination
                auto sample = samples[combination[i]];

                //multiply the parity in samples with the coefficient
                auto temp = scalar_product(sample.first, coefficient[i], field_size);
                sumParity = addition(sumParity, temp, field_size);

                //multiply the check in samples with the coefficient
                sumCheck = (sumCheck + coefficient[i]*sample.second) % field_size;
            }  
            //Sort the t_sample into a table, ready make 2t_sample
            if(sumCheck < 0){
                sumCheck += field_size;
            }
            
            
            auto t_sample = make_pair(sumParity, sumCheck);
            sort_sample(table, delta, t_sample);
            count++;
            
            //cout << "Sum: " << sumParity;
            //cout << "Check sum: " << sumCheck << endl;
        }
    }
    
 

}
/*------------------------------------------------*/

set<Pair> merge_list(const Map& list1, const Map& list2, const int& field_size, const Set& E){
    set<Pair> res;
    int total{0};
    for(auto entry1 : list1){
        auto tag = entry1.first;
        auto it = list2.find(tag);
        if(it != list2.end()){
            auto entry2 = *it;
            for(auto sample1: entry1.second){
                for(auto sample2: entry2.second){
                    auto merge_vec = subtraction(sample1.first, sample2.first, field_size);
                    auto merge_checksum = modSub(sample1.second, sample2.second, field_size);
                    auto merge_sample = make_pair(merge_vec, merge_checksum);
                    total++;
                    bool exist = false;
                    //checking dependencies.
                    for(auto e: E){
                        auto e_vec = scalar_product(merge_vec, e, field_size);
                        auto e_checksum = (merge_checksum * e) % field_size;
                        auto e_sample = make_pair(e_vec,e_checksum);
                        if(res.find(e_sample) != res.end()){
                            exist = true;
                            break;
                        }
                    }
                    if(!exist){
                        res.insert(merge_sample);
                    }
                    
                    // cout << sample1.first;
                    // cout << sample2.first;
                    // cout << merge_vec;
                    // cout << "---Checksum-----" << endl;
                    // cout << sample1.second << endl;
                    // cout << sample2.second << endl;
                    // cout << merge_checksum << endl;;
                    // cout << "---------------------" << endl;
                }
            }
        }
    }
    cout << "total:" << total << "/ Effective: " << res.size() << " / Ratio: " << total/res.size() << endl; 
    return res;
}

/*------------------------------------------------*/
/*Arimethics*/
int product(const Vec& a, const Vec& b, const int& field_size){
    int res{0};
    assert(a.size() == b.size());
    for(auto i = 0; i< a.size(); i++){
        res = (res + a[i]*b[i] % field_size) % field_size;
    }
    return res;
}





void generate_coefficients(const Set& coefficients, int t, Vec& currentCoefficient,
                           vector<Vec>& allCoefficients) {
    if (t == 0) {
        allCoefficients.push_back(currentCoefficient);
        return;
    }

    for (int coefficient : coefficients) {
        currentCoefficient.push_back(coefficient);
        generate_coefficients(coefficients, t - 1, currentCoefficient, allCoefficients);
        currentCoefficient.pop_back();
    }
}
void generate_combinations(const Vec& indices, int t, Vec& currentCombination, vector<Vec>& allCombinations, int index){
    if(currentCombination.size() == t){
        allCombinations.push_back(currentCombination);
        return;
    }
    for(int i = index; i < indices.size(); i++){
        currentCombination.push_back(indices[i]);
        generate_combinations(indices, t, currentCombination, allCombinations, i + 1);
        currentCombination.pop_back();
    }
}

