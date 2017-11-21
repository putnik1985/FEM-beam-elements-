
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Matrix.h"


int main(int argc, char *argv[]){
string filename;
string title;
int dimension = 0;
int dim;
string str;
double tolerance = 0.000001;


      for(int i=1;i<argc;++i){
      //if (strcmp(argv[i],"-file")==0) { title = "file"; filename = argv[i+1]; }
      if (strcmp(argv[i],"-tolerance") == 0) { title = "Tolerance"; 
                                               tolerance = atof(argv[i+1]); // must be sure i+1 < argc nad number
                                             } 
      }

       Matrix<double> A;
       cin >> A;
       if (!A.symmetric()) { cerr << " Matrix is not symmetric \n";
                             return 1;
                           } 

       Matrix<double> Vectors = A.apply_Jacobi_rotations(tolerance);




       cout << "dimension " << A.size() << '\n';
       cout << " Eigenvalues:\n";
       A.print_diag();
       cout << " Eigenvectors:\n";
       cout << Vectors << '\n';

return 0;
}
