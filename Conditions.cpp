
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

#include "Matrix.h"
#include "Transform.h"

#include <fstream>
#include <algorithm>

#define STIFFNESS 5 //number of words requred in the _stiffness instruction
#define DAMPING 5 //number of words requred in the _stiffness instruction
#define MASS 5 //number of words requred in the _stiffness instruction
#define GYRO 3 //number of words requred in the _stiffness instruction
#define FORCE 4 //number of words requred in the _stiffness instruction
#define DENS_TOLERANCE 0.00001

void Add_Dofs(Matrix<double>& K, const string& Dofs, const double& c, const int& i, const int& j);

double str_to_double(const string& str);
void fill_stiff_damp(istream& is, Matrix<double>& M, const string& keyword, const vector<string>& Trans);
void fill_mass(istream& is, Matrix<double>& M, const string& keyword, const vector<string>& Trans);
void Add_Dofs(vector<double>& P, const string& Dofs, const double& c, const int& i);
void fill_force(istream& is, vector<double>& P, const string& keyword, const vector<string>& Trans);
void fill_gyro(istream& is, Matrix<double>& M, const string& keyword, const vector<string>& Trans);
void fill_unb(istream& is, vector<double>& uy, vector<double>& uz, const string& keyword, const vector<string>&Trans);

int main(int argc, char *argv[]){

string title;
string line;

      if(argc < 3) {
      cerr << " Transform file [-Stiffness | Mass | Force | Gyro | Damping ] is not specified " <<'\n';
      return 1;
      }


ifstream ifs(argv[1],ios_base::in);

      if(!ifs) {
      cerr << " can not open file: " << argv[1] << '\n';
      return 1;
      }

      while(ifs>>line && line!="Transformation:");
      if(!ifs){
               cerr << " could not find Transformation vector in file: "
                    <<  argv[1] << '\n';
               return 2;
      }
int dimension = 0;
       ifs >> dimension; //read dimension of the vector
       //cout << " dimension: " << dimension << '\n';

       vector<string> Trans;

       while( getline(ifs,line) ){
          StringList station(line,';');
          if(station.size()>0) Trans.push_back(station[0]+"."+station[3]);
      }

      /*
      cout << " Function:\n";
      for(vector<string>::iterator i=Trans.begin();i!=Trans.end();++i)
      cout << i - Trans.begin() + 1 << " --> " << *i << '\n';
      */

      Matrix<double> Array(6*Trans.size());

      if (strcmp(argv[2],"-Mass")==0){
      title = "Mass";
      // fill up matrix diagonal with low density
      //for(int k=1;k<=Array.size();++k)
      //    Array(k,k) = DENS_TOLERANCE; 

      fill_mass(cin, Array, "_mass", Trans);
      } //if Mass
      if (strcmp(argv[2],"-Stiffness")==0) {
                  title = "Stiffness";
                  fill_stiff_damp(cin, Array, "_stiffness", Trans);
                    
      }// if stiffness

      if (strcmp(argv[2],"-Force")==0) { title = "Force";
      //cout << " force: \n";
      vector<double> force( 6*Trans.size() );
      fill_force(cin, force, "_force", Trans);
                 for(int i=0;i<force.size();++i)
                 cout << fixed << force[i] << ";\n";
      }

      if (strcmp(argv[2],"-Gyro")==0){
      title = "Gyro";
      fill_gyro(cin, Array,"_gyro", Trans);
      }

      if (strcmp(argv[2],"-Damping")==0){
      title = "Damping";
      fill_stiff_damp(cin, Array, "_damping", Trans);
        
      } //if damping



      if (strcmp(argv[2],"-Unbalance")==0) { title = "Unbalance";
      
      vector<double> unb_y( 6*Trans.size() );
      vector<double> unb_z( 6*Trans.size() );

      fill_unb(cin, unb_y, unb_z, "_unbalance", Trans);
               for(int i=0;i<unb_y.size();++i)
               cout << fixed << unb_y[i] << ";" << unb_z[i] << ";\n";
      }


      
      if(title!="Force" && title!="Unbalance") Array.print();

return 0;
}

void Add_Dofs(Matrix<double>& K, const string& Dofs, const double& c, const int& i = 0, const int& j = 0){

                       /*
                       cout << " DOFS: " << Dofs
                            << " Stiffness: " << c
                            << " node1= " << i
                            << " node2= " << j <<'\n';
                       // return;
                       */

                       if ( i==0 && j == 0 ) return;

                       if( i==0 && j != 0 ){

                               for(int k=0;k<Dofs.size();++k){
                               int dof = Dofs[k] - '0';
                               K(6*j-6+dof,6*j-6+dof) = K(6*j-6+dof,6*j-6+dof) + c;
                               }

                       }

                       if( i != 0 && j == 0 ){

                               for(int k=0;k<Dofs.size();++k){
                               int dof = Dofs[k] - '0';
                               K(6*i-6+dof,6*i-6+dof) = K(6*i-6+dof,6*i-6+dof) + c;
                               }

                       }


                       if( i != 0 && j != 0 ){

                               for(int k=0;k<Dofs.size();++k){
                               int dof = Dofs[k] - '0';
                               //cout << dof << '\n'; 
                               
                               K(6*i-6+dof,6*i-6+dof) = K(6*i-6+dof,6*i-6+dof) + c;
                               K(6*j-6+dof,6*j-6+dof) = K(6*j-6+dof,6*j-6+dof) + c;
                               K(6*i-6+dof,6*j-6+dof) = K(6*i-6+dof,6*j-6+dof) - c;
                               K(6*j-6+dof,6*i-6+dof) = K(6*j-6+dof,6*i-6+dof) - c;
                               
                               }

                       }
}

double str_to_double(const string& str){
istringstream iss(str,ios_base::in);
double d;
       iss >> d;
     //  cout << DBL_MAX << '\n';
       cout << iss.str()<< " is " << d << '\n';
return d;
}


void fill_stiff_damp( 
                     istream& is, 
                     Matrix<double>& M, 
                     const string& keyword, 
                     const vector<string>& Trans
                    )
{
string line;
int node1 = 0, node2 = 0;

                  while(getline(is,line)){
                  StringList list(line);
                      if( list.size() == STIFFNESS
                                              && list[0]==keyword ){
                     
                      cout << "node 1: " << list[1] << " "
                           << "node 2: " << list[2] << " "
                           << "Dofs : "  << list[3] << " "
                           << keyword << " : " << list[4] << '\n';
                    
 
                   // find nodes
                      vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),list[1]);
                     (p==Trans.end()) ? node1 = 0 : node1 = p - Trans.begin() + 1;

                      p = find(Trans.begin(),Trans.end(),list[2]);
                     (p==Trans.end()) ? node2 = 0 : node2 = p - Trans.begin() + 1;

                    //cout << "n1= " << node1 << " n2= " << node2 << '\n';
                      double c = atof(list[4].c_str()); // does not convert correctly 500000.12234 --> 500000
                      //cout << "c= " << c << '\n';

                    Add_Dofs(M, list[3], c, node1, node2);

                      } // string with stiffness proccesing
                  }// while stiffness is found

}


void fill_mass( 
                     istream& is, 
                     Matrix<double>& M, 
                     const string& keyword, 
                     const vector<string>& Trans
                    )
{
string line;
int node = 0;

                  while(getline(is,line)){
                  StringList list(line);
                      if( list.size() == STIFFNESS
                                              && list[0]==keyword ){
                     
                      cout << "node  : " << list[1] << " "
                           << "Mass : " << list[2] << " "
                           << "Jp : "  << list[3] << " "
                           << "Jd : "  << list[4] << '\n';
                    
 
                   // find nodes
                      vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),list[1]);
                     (p==Trans.end()) ? node = 0 : node = p - Trans.begin() + 1;


                    //cout << "n1= " << node1 << " n2= " << node2 << '\n';
                    double c = atof(list[2].c_str()); // does not convert correctly 500000.12234 --> 500000
                    // add mass
                    Add_Dofs(M, "123", c, node);

                    c = atof(list[3].c_str()); // does not convert correctly 500000.12234 --> 500000
                    // add Jp
                    Add_Dofs(M, "4", c, node);


                    c = atof(list[4].c_str()); // does not convert correctly 500000.12234 --> 500000
                    // add Jd
                    Add_Dofs(M, "56", c, node);


                      } // string with stiffness proccesing
                  }// while stiffness is found

}


void Add_Dofs(vector<double>& P, const string& Dofs, const double& c, const int& i){

                       if( i != 0 )
                               for(int k=0;k<Dofs.size();++k){
                               int dof = Dofs[k] - '0';
                               P[6*i-6+dof - 1] = P[6*i-6+dof - 1] + c; // minus 1 because it is a vector
                               }

}


void fill_force(istream& is, vector<double>& P, const string& keyword, const vector<string>& Trans){

string line;
int node1 = 0;

                  while(getline(is,line)){
                  StringList list(line);
                      if( list.size() == FORCE
                                              && list[0]==keyword ){
                     
                      cout << "node  : " << list[1] << " "
                           << "Dofs : "  << list[2] << " "
                           << keyword << " : " << list[3] << '\n';
                    
 
                   // find nodes
                        vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),list[1]);
                       (p==Trans.end()) ? node1 = 0 : node1 = p - Trans.begin() + 1;

                    //cout << "n1= " << node1 << " n2= " << node2 << '\n';
                      double c = atof(list[3].c_str()); // does not convert correctly 500000.12234 --> 500000
                    //cout << "c= " << c << '\n';
                      if(node1 > 0)
                      Add_Dofs(P, list[2], c, node1);
                      else cout << "node: " << list[1] << " not found\n";
                      } // string with stiffness proccesing
                  }// while stiffness is found

}


void Add_Gyro(Matrix<double>& K, const double& c, const int& i){
                     K(6*i-6 + 5,6*i-6 + 6) = K(6*i-6 + 5,6*i-6 + 6) + c;
                     K(6*i-6 + 6,6*i-6 + 5) = K(6*i-6 + 6,6*i-6 + 5) - c;
}

void fill_gyro(istream& is, Matrix<double>& M, const string& keyword, const vector<string>& Trans){

string line;
int node1 = 0;

                  while(getline(is,line)){
                  StringList list(line);
                      if( list.size() == GYRO
                                              && list[0]==keyword ){
                     
                      cout << "node  : " << list[1] << " "
                           << keyword << " : " << list[2] << '\n';
                    
 
                   // find nodes
                        vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),list[1]);
                       (p==Trans.end()) ? node1 = 0 : node1 = p - Trans.begin() + 1;

                    //cout << "n1= " << node1 << " n2= " << node2 << '\n';
                      double c = atof(list[2].c_str()); // does not convert correctly 500000.12234 --> 500000
                    //cout << "c= " << c << '\n';
                      if(node1 > 0)
                      Add_Gyro(M, c, node1);
                      else cout << "node: " << list[1] << " not found\n";
                      } // string with stiffness proccesing
                  }// while stiffness is found

}


void fill_unb(istream& is, vector<double>& uy, vector<double>& uz, const string& keyword, const vector<string>&Trans){

string line;
int node1 = 0;

                  while(getline(is,line)){
                  StringList list(line);
                      if( list.size() == FORCE
                                              && list[0]==keyword ){
                     
                      cout << "node  : " << list[1] << " "
                           << "Phase : "  << list[3] << " "
                           << keyword << " : " << list[2] << '\n';
                    
 
                   // find nodes
                    vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),list[1]);
                    (p==Trans.end()) ? node1 = 0 : node1 = p - Trans.begin() + 1;

                    //cout << "n1= " << node1 << " n2= " << node2 << '\n';
                     double unb = atof(list[2].c_str()); // does not convert correctly 500000.12234 --> 500000
                     double phase = atof(list[3].c_str()); // does not convert correctly 500000.12234 --> 500000
                     phase *= M_PI / 180; // degrees in radians
                     // cout << "c= " << unb << " phase " << phase << '\n';
                     if(node1 > 0) // add unbalance to the Dofs 2 and 3
                     {
                      uy[ 6 * node1 - 6 + 2 - 1]+=unb * cos(phase);
                      /////cout <<  unb * sin(phase) << '\n';
                      uy[ 6 * node1 - 6 + 3 - 1]+=unb * sin(phase);
                      uz[ 6 * node1 - 6 + 2 - 1]+=-unb * sin(phase);
                      uz[ 6 * node1 - 6 + 3 - 1]+=unb * cos(phase);

                     }               
                      else cout << "node: " << list[1] << " not found\n";
                      } // string with stiffness proccesing
                  }// while stiffness is found



}
