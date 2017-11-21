
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "Matrix.h"
#include "Modes.h"


int main(int argc, char *argv[]){

   if(argc < 2) { cerr << " Transform file is required \n"; return 1; }

   string filename;

   vector<string> Trans;
 
          ifstream in(argv[1]);
          if(!in) { cerr << " can not open file: " << argv[1] << '\n'; return 1; }

          while(getline(in,filename)){
          string word;
          istringstream iss(filename);
          iss >> word;
          if (word=="Transformation:") break;
          }

          while(getline(in,filename)){
          StringList list(filename,';');
          Trans.push_back(list[0]+'.'+list[3]);
          }

          for(int i=0;i<Trans.size();++i)
          cout << Trans[i] <<'\n';
          
          in.close();
      
return 0;
}
