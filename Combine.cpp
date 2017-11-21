
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Matrix.h"
#include "Transform.h"

#include <fstream>


int main(int argc, char *argv[]){
string filename;
vector<string> files;
ifstream is; 
string title;
int dimension = 0;
int dim;
string str;
string line;
int Kd, Td;


      if(argc < 2) {
      cerr << " Matrix [-Stiffness] [-Mass] [-Transform] is not specified " <<'\n';
      return 1;
      }


      for(int i=1;i<argc;++i){
      if (strcmp(argv[i],"-Mass")==0) title = "Mass";
      if (strcmp(argv[i],"-Stiffness")==0) title = "Stiffness";
      if (strcmp(argv[i],"-Transformation")==0) title = "Transformation";
      }
      
      // it is crucial consider if it wrong input names

      if (title.length()==0){
      cerr << " Matrix [-Stiffness] [-Mass] [-Transform] is not specified " <<'\n';
      return 2;
      }

       while (cin >> filename)
       files.push_back(filename);


       //cout << " number of files read: " << files.size() << '\n';
       for(
          vector<string>::iterator p = files.begin();
          p != files.end(); 
          ++p
          ){
              is.open((*p).c_str(),ios_base::in);
              if(is){
              is >> Kd >> Td;
              cout << " file: " << *p << " dimension: " << Kd << '\n';
              if(strcmp(title.c_str(),"Transformation")==0) dimension+=Td;
              else  dimension+=Kd;  
              is.close();
               } else {
                 cout << " can not open file: " << *p << '\n';
               }
       }

       cout << '\n';
       cout << title << ": " << dimension << '\n';

int index = 0; // indicates where to insert submatrix
       for(
          vector<string>::iterator p = files.begin();
          p != files.end(); 
          ++p
          ){
              is.open((*p).c_str(),ios_base::in);
              if(is){
              //cout << " process file: " << *p << '\n';
              while( getline(is,line) ){
                                       //cout << line << '\n';
                                       StringList list(line,':');
                                       if (list[0]==title) { 
                                                            dim = atoi(list[1].c_str());
                                                            //cout << " dimension: " << dim << '\n';
                                                            //cout << " lets process the lines " << '\n';
                                                            if(strcmp(title.c_str(),"Transformation")==0) {
                                                             
                                                            //cout << " process transformation file: " << *p << '\n';
                                                            for(int i=1;i<=dim;++i){
                                                            getline(is,line);
                                                            cout << *p <<';'<< index + i <<';'<<line<<'\n';
                                                            }
                                                            index+=dim;
                                                            continue;
                                                            }
                                                            // scaning Stiffness and Mass
                                                            //cout << line << '\n';
                                                                for(int i=0;i<dim;++i){
                                                                getline(is,line); // read line from matrix;
                                                                StringList strlist(line,';');
                                                                StringList zeros("0",dimension,';'); // create zero string only 0;
                                                                zeros.insert(index,strlist); // insert line read into zero list
                                                                zeros.print();
                                                                cout << '\n'; 
                                                                //cout << line << '\n';   
                                                                }
                                                           index+=dim; // next matrix should be inserted  
                                                           }   
                                       }   
              is.close();
              } 
       }

/*
StringList list("0.0",6,';');
StringList list2("1.0",9,';');

list.print(); 
cout << '\n';
list2.print();

list.insert(2,list2);
cout << '\n';

list.print();
   */
   

return 0;
}
