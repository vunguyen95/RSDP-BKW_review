#ifndef VECTOR_H
#define VECTOR_H
#include<bits/stdc++.h>
#include<vector>
using namespace std;

set<vector<int>> samplesGen(set<vector<int>>& res, int& samples, int& length , int& weight);
vector<vector<int>> getComb(vector<int>& v, int n, int& r);
void combination(vector<int>& v, vector<int>& data, int start, int end, int index, int r, vector<vector<int>>& res);

/* vector in binary form*/
vector<char> represent(vector<int>& v, int& n);

/* vector in representation form*/
vector<int> bucketRepresent(vector<char>& v);
vector<char> add(vector<char>& u, vector<char>& v);
int weight(vector<char>& v);
vector<vector<char>> parityGen(int& l, int& n);
/*Scalar product with vectors in representation form*/
int rscalar(const vector<int>& v1, const vector<int>& v2);

/*XOR in representation form*/
vector<int> radd(vector<int>& v1, vector<int>& v2);

/*Multiplication with vectors in representation form*/
vector<int> mul(const vector<vector<int>>& a, const vector<int>& b);

/*Multiplication in binary form form*/
vector<char> mul1( vector<vector<char>>& a, vector<char>& b);
vector<vector<int>> get(vector<vector<int>>& parity, int& i);
#endif
