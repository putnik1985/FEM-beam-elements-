
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Matrix.h"


int main(int argc, char *argv[]){
string title="";

      if(argc > 1){
      if (strcmp(argv[1],"-mult")==0) { title = "mult"; }
      if (strcmp(argv[1],"-sum") == 0) {  title = "sum";        } 
      }

      string filename;
      double d;
      
      cin >> d >> filename;
      ifstream ifs(filename.c_str(), ios_base::in);
          if(!ifs) { cerr << " can not open file: " << filename << '\n'; return 1; }
          cout << filename << " is processed " << fixed << d << '\n';
          Matrix<double> M;
          ifs >> M;
          M = d * M;
          ifs.close();

      //M.print();
      //exit(3);
      

      while(cin >> d >> filename){
      // cout << d << " -- " << filename << '\n';
          
          ifs.open(filename.c_str(), ios_base::in);
          if(!ifs) { cerr << " can not open file: " << filename << '\n'; return 1; }
          Matrix<double> A;
          cout << filename << " is processed " << fixed << d << '\n';
          ifs >> A;
         // cout << " Matrix:\n ";
         // A.print();
          ifs.close();
             
            (title != "mult")   ? 
            M = M + d * A       :
            M = d * M * A       ;

            //cout << " Matrix:\n "; 
            //M.print();
      }
       
      cout << " dimension " <<  M.size() << '\n';
      M.print();
      //M.transpose();
      //M.print();



return 0;
}
