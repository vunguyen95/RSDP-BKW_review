#include"headers/misc.h"
std::ostream& operator<< (ostream& os, const vector<vector<char>>& v)
{
        for (int i = 0; i < v.size(); i++)
        {
                for(int j = 0; j < v[i].size(); j++)
                {
                        os << static_cast<int>(v[i][j]) << " ";
                }
                os << endl;
        }
        return os;
}
std::ostream& operator<< (ostream& os, const vector<char>& v)
{
	for(int i = 0; i < v.size(); i++)
	{
		os << static_cast<int>(v[i]) << " ";
	}
	os << endl;
	return os;
}
std::ostream& operator<< (ostream& os, const vector<int>& v)
{
	for(int i = 0; i < v.size(); i++)
	{
		os << v[i] << " ";
	}
	os << endl;
	return os;
}
std::ostream& operator<< (ostream& os, const vector<vector<int>>& v)
{
        for (int i = 0; i < v.size(); i++)
        {
                for(int j = 0; j < v[i].size(); j++)
                {
                        os << (v[i][j]) << " ";
                }
                os << endl;
        }
        return os;
}
