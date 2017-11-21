
#ifndef READ_TRANSLAT_H
#define READ_TRANSLAT_H

#include <assert.h>
#include "Matrix.h"

#define DENS_TOLERANCE 0.00001

struct Station{
int station;
double  x;
double Ro;
double Ri;
double E;
double Rho;
double Nu;
double Mass;
double Ip;
double Id;

     double JD() const { return (M_PI/2)*(Ro*Ro*Ro*Ro-Ri*Ri*Ri*Ri); }
     double JY() const { return (M_PI/4)*(Ro*Ro*Ro*Ro-Ri*Ri*Ri*Ri); }
     double JZ() const { return (M_PI/4)*(Ro*Ro*Ro*Ro-Ri*Ri*Ri*Ri); }
     double F() const { return (M_PI)*(Ro*Ro-Ri*Ri); }
     double G() const { return E / (2 * (1 + Nu)); }

Station& operator=(const Station& st2){
x = st2.x;
Ro = st2.Ro;
Ri = st2.Ri;
E = st2.E;
Rho = st2.Rho;
Nu = st2.Nu;
Mass = st2.Mass;
Ip = st2.Ip;
Id = st2.Id;
return *this;
}

}; //responsible for the calculation of the station parameters 

istream& operator>>(istream& is, Station& st){
     is >> st.station 
        >> st.x
        >> st.Ro
        >> st.Ri
        >> st.E
        >> st.Rho
        >> st.Nu 
        >> st.Mass
        >> st.Ip
        >> st.Id; 

  if( st.Rho < DENS_TOLERANCE ) st.Rho = DENS_TOLERANCE; //can not allow zero density 

return is;
}

ostream& operator<<(ostream& os, const Station& st){

  os << st.station <<";"<<st.x<<";";

return os;
}

Matrix<double> FORM_LOCAL_KE_BEAM(
                          double F,
                          double JD,
                          double JY,
                          double JZ,
                          double L,
                          double E,
                          double G
                         );

Matrix<double> FORM_LOCAL_ME_BEAM(double F, double L, 
                                  double RO, double JD
                                  );

int Add_Global(const int& i,const int& j, const Matrix<double>& KE ,Matrix<double>& K);

int Create_FEM(istream& is){
char c;
string s;
vector<Station> Trans;
Matrix<double> K;
Matrix<double> M;
int nodes = 0;
Station st1, st2;
double Jd1, Jd2, Jy1, Jy2, Jz1, Jz2;


     //first line is titles
     getline(is,s);  
     //cout << " title: " << s << '\n';

     //read first line of the mesh
     is >> st1;
     Trans.push_back(st1);
     ++nodes;

     while(is >> st2){
     ++nodes; 
     //cout << "x= " << st2.x << '\n';

     Trans.push_back(st2);
     assert( st2.x - st1.x > 0 );

     Matrix<double> ke = FORM_LOCAL_KE_BEAM(
                          st1.F(),
                          st1.JD(),
                          st1.JY(),
                          st1.JZ(),
                          st2.x - st1.x,
                          st1.E,
                          st1.G()
                         );
/*
     int direction = 3;
     double s =0.0;
     for(int k = 1; k <= 12;++k){
     s = ke(k,direction) + ke(k,direction +6);
     //cout << " a = " << ke(k,direction) << '\n'; 
     //cout << " b = " << ke(k,direction +6) << '\n'; 
     cout << " sum = " << s << '\n'; 
     } 
*/
     //$%#$%#$%#$%#$%#$%#$%#$%$#
     //  cout << " dimension: " << ke.size() 
     //       << " node: " << Trans.size()-1
     //       << " node: " << Trans.size() << '\n';
     //       ke.print();
     //34%#$%#$%#$%#%#%#%#$%#$%



       Matrix<double> me = FORM_LOCAL_ME_BEAM(
                                       st1.F(), st2.x - st1.x, 
                                       st1.Rho, st1.JD()
                                       );

     //cout << "size: " << Trans.size() <<'\n';
     //cout << " nodes: " << nodes << '\n';

     //$%#$%#$%#$%#$%#$%#$%#$%$#
       //cout << " dimension: " << K.size() << '\n';
       //K.print();
     //34%#$%#$%#$%#%#%#%#$%#$%


     K.set_dimension( 6 * Trans.size());
     M.set_dimension( 6 * Trans.size()); 

     //$%#$%#$%#$%#$%#$%#$%#$%$#
      //cout << " dimension: " << K.size() << '\n';
     /// K.print();
     //34%#$%#$%#$%#%#%#%#$%#$%

     //$%#$%#$%#$%#$%#$%#$%#$%$#
     ////cout << " dimension: " << K.size() << '\n';
     ////K.print();
     //34%#$%#$%#$%#%#%#%#$%#$%     
  
     Add_Global(Trans.size()-1,Trans.size(),ke ,K);
     Add_Global(Trans.size()-1,Trans.size(),me ,M);

     //$%#$%#$%#$%#$%#$%#$%#$%$#
     ////cout << " dimension: " << K.size() << '\n';
     ///K.print();
     //34%#$%#$%#$%#%#%#%#$%#$%

     st1 = st2;
     //cout << " eof: " << is.eof() << '\n';

     }
     cout << 6 * Trans.size() << " " 
          << Trans.size() << '\n';

     cout << " Stiffness: " << 6 * Trans.size() << '\n';
     K.print();
     cout << " Mass: " << 6 * Trans.size() << '\n';
     M.print();
     cout << " Transformation: " << Trans.size() << '\n';
     for(int i = 0; i<Trans.size();++i)
     cout << i + 1 <<";"<< Trans[i] << '\n';

return 0;
}


Matrix<double> FORM_LOCAL_KE_BEAM(
                          double F,
                          double JD,
                          double JY,
                          double JZ,
                          double L,
                          double E,
                          double G
                         )
{
       Matrix<double> KE(12);
       int I, J;

       KE(1,1) = E * F / L ; KE(1,7) = - E * F / L;

       KE(2,2) = 12 * E * JZ / (L*L*L);
       KE(2,6) = 6 * E * JZ / (L*L) ;
       KE(2,8) = -12 * E * JZ / (L*L*L); 
       KE(2,12) = 6 * E * JZ / (L*L) ;
	
       KE(3,3) = 12 * E * JY / (L*L*L);
       KE(3,5) = -6 * E * JY / (L*L);
       KE(3,9) = -12 * E * JY / (L*L*L);
       KE(3,11) = -6 * E * JY / (L*L);
	
       KE(4,4) = G * JD / L;
       KE(4,10) = - G * JD / L;

       KE(5,5) = 4 * E * JY / L;
       KE(5,9) = 6 * E * JY / (L*L);
       KE(5,11) = 2 * E * JY / L;

       KE(6,6) = 4 * E * JZ / L;
       KE(6,8) = -6 * E * JZ / (L*L);
       KE(6,12) = 2 * E * JZ / L;

       KE(7,7) = E * F / L;

       KE(8,8) = 12 * E * JZ / (L*L*L);
       KE(8,12) = -6 * E * JZ / (L*L);

       KE(9,9) = 12 * E * JY / (L*L*L);
       KE(9,11) = 6 * E * JY / (L*L);

       KE(10,10) = G * JD / L;
       
       KE(11,11) = 4 * E * JY / L;
       KE(12,12) = 4 * E * JZ / L;

 
       for (I=1; I<=12; ++I)
       for ( J=I; J<=12; ++J)
       KE(J,I)=KE(I,J);


return KE;							 
}


Matrix<double> FORM_LOCAL_ME_BEAM(double F, double L, 
                                  double RO, double JD
                                  )
{
       Matrix<double> ME(12);
       int I,J;

       ME(1,1) = F * L * RO / 3;
       ME(1,7) = F * L * RO / 6;

	ME(2,2) = 13 * F * L * RO / 35;  
        ME(2,6) = 11 * F * (L*L) * RO / 210 ;
        ME(2,8) = 9 * F * L * RO / 70 ;
	ME(2,12) = -13 * RO * F * (L*L) / 420 ;
	
	ME(3,3) = 13 * F * L * RO / 35 ;
        ME(3,5) = -11 * F * (L*L) * RO / 210; 
	ME(3,9) = 9 * F * L * RO / 70 ;
        ME(3,11) = 13 * F * (L*L) * RO / 420; 

        ME(4,4) = JD*RO*L*0.5; 
        
	ME(5,5) = F * (L*L*L) * RO / 105; 
        ME(5,9) = -13 * F * (L*L) * RO /420; 
	ME(5,11) = -F * (L*L*L) * RO / 140; 

	ME(6,6) = F * (L*L*L) * RO / 105 ;
        ME(6,8) = 13 * F * (L*L) * RO / 420;
	ME(6,12) = -F * (L*L*L) * RO / 140 ;

	ME(7,7) = F * L * RO / 3 ;

	ME(8,8) = 13 * F * L * RO / 35 ;
        ME(8,12) = -11 * F * (L*L) * RO / 210; 

	ME(9,9) = 13 * F * L * RO / 35 ;
        ME(9,11) = 11 * F * (L*L) * RO / 210 ;

        ME(10,10) = JD*RO*L*0.5;
 
	ME(11,11) = F * (L*L*L) * RO / 105 ;
	ME(12,12) = F * (L*L*L) * RO / 105 ;


	for (I=1; I<=12;++I)
        for (J=I;J<=12;++J)
        ME(J,I)=ME(I,J);
/*
        double m = RO * F * L;
        ME(1,1) = m/2;
        ME(2,2) = m/2;        
        ME(3,3) = m/2;

        ME(7,7) = m/2;
        ME(8,8) = m/2;        
        ME(9,9) = m/2;
*/
return ME;
}

int Add_Global(const int& i,const int& j, const Matrix<double>& KE ,Matrix<double>& K)
{
      int gi, gj; // global DOFS
      vector<int> trans(12); // for beam elements only

      for(int dof = 1;dof<=6;++dof){
          trans[dof-1] = 6 * i - 6 + dof;
          trans[dof-1 + 6] = 6 * j - 6 + dof;
      }

      //cout << " input: \n" ;
     // KE.print();
     // cout << " sum: \n";

      for(int r=1; r<=12; ++r){
          gi = trans[r - 1];
          for(int c = 1; c<=12; ++c){ 
          gj = trans[c - 1];
          //cout << fixed << K(gi, gj) 
          //     << " + " << fixed << KE(r,c) << " = ";
          K(gi, gj) = K(gi, gj) + KE(r,c);
          //cout << fixed << K(gi,gj) << '\n';            
          }                       
     }
     //cout << " after sum:\n";
     //K.print();
return 0;
}


#endif
