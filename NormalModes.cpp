
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "Matrix.h"
#include "Modes.h"


int main(int argc, char *argv[]){

  if (argc<3) { cerr << " Mass eigenvectors  & Stiffness eigenvectors [-print] must be specified\n"; return 1; }

  //read Mass matrix eigenvectors
  // Matrix mass must be positive defiened
 
     ifstream ifs(argv[1]);
     if(!ifs) { cerr << " can not open file: " << argv[1] << '\n'; return 2; }
     Mode mode;
   
     cout << " read Mass vectors \n";

     vector<Mode> modes;
     ifs >> modes;
     ifs.close(); // complete Mass modes reading

   
     vector<Mode> K_modes;

       ifs.open(argv[2]);
       if(!ifs) { cerr << " can not open file: " << argv[2] << '\n'; return 3; }

       cout << " read Stiffness vectors \n";
       ifs >> K_modes;
       ifs.close(); // compltete K modes reading

       Matrix<double> F(modes.size());
           for(int i=0;i<modes.size();++i)
               F.fill_column(modes[i].mode,i+1); 


       Matrix<double> KF(K_modes.size());
           for(int i=0;i<K_modes.size();++i)
               KF.fill_column(K_modes[i].mode,i+1); 
 
       // complete reading Mass eigen values
       //F.print(); 
      // return 8;


       Matrix<double> FT = F;
       FT.transpose();

       Matrix<double> KFT = KF;
       KFT.transpose();
       //FT.print();

       Matrix<double> D( modes.size() ); 
            for(int i=0;i<modes.size();++i){
            
            D(i+1, i+1) = modes[i].frequency;
            } 

       Matrix<double> invKD( K_modes.size() ); 
            for(int i=0;i<K_modes.size();++i){
             K_modes[i].frequency = abs(K_modes[i].frequency);
             assert(K_modes[i].frequency > 0);
             invKD(i+1, i+1) = 1/ sqrt(K_modes[i].frequency);
            }

       //KD.print();
 
       
       assert(modes.size() == F.size());

       Matrix<double> A(modes.size());
      /*
       assert(invD.symmetric());
       F.balance();
       FT.balance();
       KF.balance();
       KFT.balance();
      */
       // K = F * D * FT
       // M = F * D * FT

       A = invKD * KFT * F * D * FT * KF * invKD;  // Usually this matrix is not exactly symmetric
       
       if (argc >= 4 && strcmp(argv[3],"-print")==0 ) { A.print(); return 0; }

       /////A = invD;


       //assert(A.symmetric());
       //A.print();
       //return 8;

       modes.clear();
       F.clear();

       mode.clear();
       K_modes.clear();

       KFT.clear();
       D.clear();


       F = A.apply_Jacobi_rotations();
      // normalisation
       for(int i=1;i<=F.size();++i)
         for(int k=1;k<=F.size();++k)
         F(k,i)/=sqrt( abs( A(i,i) ) );

   
      /// F.print();
      /// return 8;
       KF = KF * invKD * F; // calculate actual modes

        for(int i=1;i<=A.size();++i){
         mode.frequency = (1/(2*M_PI)) * sqrt(1 / abs(A(i,i)) );
         modes.push_back(mode);
       }

       for( int i=0; i < modes.size(); ++i)
       modes[i].mode = KF.column(i+1);

       sort(modes.begin(),modes.end());

       cout << "dimension " << modes.size() << '\n';
       cout << " Eigenvalues:\n";
       for(int i=0; i<modes.size();++i)
       cout << fixed << setprecision(WRITE_PRECISION) << modes[i].frequency << ";\n"; 
       cout << " Eigenvectors:\n";
       for( int i=0; i < modes.size(); ++i)
       KF.fill_column(modes[i].mode, i+1);
       KF.print();
      
return 0;
}
