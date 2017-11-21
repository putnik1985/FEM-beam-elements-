
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "Matrix.h"
#include "Modes.h"
#include "Read_Translat_Mesh.h"

int main(int argc, char *argv[]){

    if (argc < 2) { cerr << " Transform indexes should be specified'\n"; return 1; }

    vector<string> Trans;
           
    ifstream in(argv[1]);
    if (!in) { cerr << " can not open file: " << argv[1] << '\n'; }

    string word;
    while(in >> word) Trans.push_back(word);

      //for(int i=0;i<Trans.size();++i)
      //cout << Trans[i] << '\n';
     in.close();   
     string name, dofs;
     int node;
                        while (cin >> name >> dofs){
                        ////cout << name << " " << dofs << '\n'; 
                        vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),name);
                        (p==Trans.end()) ? node = 0 : node = p - Trans.begin() + 1;
                                if ( node > 0 ) {
                                              for(int k=0;k<dofs.size();++k){
                                                  int dof = dofs[k] - '0';
                                                  cout << name << " " << dof << " " << 6 * node -6 + dof << '\n';
                                              } 
                                } else cout << name << " " << " not found\n ";
                        }
 


return 0;
}
