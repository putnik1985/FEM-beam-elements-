
#ifndef MODES_H
#define MODES_H

#include "Transform.h"
#include "Matrix.h"

struct Mode{
double frequency;
vector<double> mode;
void print() {
             cout << " frequency: " << frequency <<'\n';
             for(int i=0;i<mode.size();++i)
                 cout << mode[i] <<'\n';  
             } 
void clear() { mode.clear(); }




}; // responsible for the operation with the mode

bool operator<(const Mode& m1, const Mode& m2){
return m1.frequency < m2.frequency;
}

istream& operator>>(istream& ifs, vector<Mode>& modes){
string words;
Mode mode;

     while (getline(ifs,words) && !isdigit(words[0]) && words[0]!='-');
     ///cout << words << '\n';
     double freq;
     StringList list(words,';');
     //list.print();
     //exit(99);         
     mode.frequency = atof(list[0].c_str());    
     modes.push_back(mode);

         while (getline(ifs,words) && 
                                  ( isdigit(words[0]) || words[0] == '-' )
               )
        {
        list.set_string(words,';');
        //cout << list[0] << '\n';
        mode.frequency = atof(list[0].c_str());
        modes.push_back(mode);
        }          

       ////cout << modes.size() << '\n';

       while (getline(ifs,words) && !isdigit(words[0]));
       //cout << words << '\n';
       
       list.set_string(words,';');
       // list.print();
       //#####################################3
       //cout << list.size() <<'\n';
       //cout << modes.size() << '\n';
       //#####################################3

       assert(list.size()==modes.size());
       for(int i=0;i<list.size();++i)
       modes[i].mode.push_back( atof(list[i].c_str()) );

         while (getline(ifs,words) && 
                                  ( isdigit(words[0]) || words[0] == '-' )
               )
        {

       list.set_string(words,';');
       for(int i=0;i<list.size();++i)
       modes[i].mode.push_back( atof(list[i].c_str()) );

        } 

return ifs;
}


#endif
