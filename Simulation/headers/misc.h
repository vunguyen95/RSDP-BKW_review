#ifndef MISC_H
#define MISC_H
#include<iostream>
#include<vector>
using namespace std;
std::ostream& operator<< (ostream& os, const vector<vector<char>>& v);
std::ostream& operator<< (ostream& os, const vector<vector<int>>& v);
std::ostream& operator<< (ostream& os, const vector<char>& v);
std::ostream& operator<< (ostream& os, const vector<int>& v);
#endif
