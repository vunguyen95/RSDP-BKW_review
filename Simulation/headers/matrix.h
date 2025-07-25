// matrix operations
#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
#include<vector>
#include<cstdlib>
#include<ctime>
using namespace std;
vector<vector<char>> mmul( vector<vector<char>>& u, vector<vector<char>>& v);
vector<char> mul( vector<char>& b,  vector<vector<char>>& a );
char scalar(vector<char>& u, vector<char>&v);
vector<vector<char>> keygen(int& n, int& m, int& b);
vector<char> add(vector<char>& u, vector<char>& v, int& length);
int weight(vector<char>& v, int& length);
#endif
