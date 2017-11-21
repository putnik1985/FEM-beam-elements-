#ifndef TRANS_H
#define TRANS_H

#include "stl.h"

struct FEM_Node{
string filename;
int number;
double x;
}; //responsible for the operations with the nodes


class StringList{
vector<string> v;
char delim;

public:
StringList(const string& str, int n, char c = ' '):delim(c) { for(int i=0;i<n;++i) v.push_back(str); }

StringList(const vector<string>& in): delim(' ')
                                     {
                                     for(int i=0;i<in.size();++i)
                                     v.push_back(in[i]);
                                     }
StringList(string& str, char c = ' '): // this function changes the input string
           delim(c)
           {
           for(int i=0;i<str.size();i++)
           if (str[i]==c) str[i] = ' ';
           
           //cout << " const: " << str <<'\n';
           istringstream iss(str,ios_base::in);
           string word;
           while (iss >> word) v.push_back(word);

           }
void set_string(string& str, char c = ' '){
delim =c;
           v.clear();
           for(int i=0;i<str.size();i++)
           if (str[i]==c) str[i] = ' ';
           
           //cout << " const: " << str <<'\n';
           istringstream iss(str,ios_base::in);
           string word;
           while (iss >> word) v.push_back(word);

} 

string& operator[](int i){ return v[i];}
const string& operator[](int i) const {return v[i];}

int size() const {return v.size();}

void print() { for(int i=0;i<v.size();++i) cout << v[i] << delim ; }
void insert(const int& pos,const StringList& list){
                                                   /////cout << list.size() << '\n';
                                                   for(int i=pos, j=0;
                                                       i<v.size() && j<list.size();
                                                       ++i, ++j
                                                       ) v[i] = list[j];

                                                  } 
vector<double>  numbers() const { vector<double> vec;
                                  for(int i=0;i<v.size();++i){
                                  // cout << " number: " << atof(v[i].c_str()) << '\n';  
                                  vec.push_back( atof(v[i].c_str()) );
                                  }
                                  return vec;
                                }

};



#endif
