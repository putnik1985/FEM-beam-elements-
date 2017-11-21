
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "Matrix.h"
#include "Modes.h"
#include "Read_Translat_Mesh.h"

void set_precision(Matrix<double>& m, const int precision){

     
     int n = m.size();
  
         for(int i=1;i<=n;++i)
          for(int j=1;j<=n;++j){
             {stringstream ss;
              ss << fixed << setprecision(precision) << m(i,j);
              cout << ss.str() << '\n';
              ss >> m(i,j);
              ss.clear(); }
          }
}


int main(int argc, char *argv[]){

    if (argc < 4) { cerr << " Transform indexes [ -static | -modal | -transient ] Solution file should be specified [-number number] [-Translat] \n"; return 1; }

    vector<string> Trans;
    bool Translat = false;       
    ifstream in(argv[1]);
    if (!in) { cerr << " can not open file: " << argv[1] << '\n'; }

    string word;
    while(in >> word) Trans.push_back(word);

      //for(int i=0;i<Trans.size();++i)
      //cout << Trans[i] << '\n';
     in.close();   


      /// define vector for the global nodes

                    string filemesh, matrixfile;
                    vector<Station> nodes;
                    string line;

                    while (cin >> filemesh >> matrixfile) { // start reading the part files

                    cout << "mesh: " << filemesh << " matrix: " << matrixfile << '\n';
                    in.open(filemesh.c_str());
 
                     if(!in) throw runtime_error("can not open file: " + filemesh); 
                     
                     Station st;
                     
                    
                        while(getline(in,line)){
                        string word;
                        istringstream iss(line);
                        iss >> word;
                        if( isdigit(word[0]) || word[0]=='-') break;
                        } 

                        //cout << line << '\n'; 
                        istringstream iss(line);
                        iss >> st;
                        string str_number = matrixfile+'.'+to_str<int>(st.station);
                        ///cout << str_number << '\n';
                        vector<string>::const_iterator p = find(Trans.begin(),Trans.end(),str_number);
                        (p==Trans.end()) ? st.station = 0 : st.station = p - Trans.begin() + 1;
                        //cout << st.station << '\n'; 
                        nodes.push_back(st);
 
                        while( in >> st){
                        str_number = matrixfile+'.'+to_str<int>(st.station);
                        /// cout << str_number << '\n';
                        p = find(Trans.begin(),Trans.end(),str_number);
                        (p==Trans.end()) ? st.station = 0 : st.station = p - Trans.begin() + 1;
                        //cout << st.station << '\n'; 
                        nodes.push_back(st);
                        }
                        //for(int i=0;i<nodes.size();++i)
                        //cout << nodes[i] << '\n';
                        in.close();          
 
                       }


      /// define vector for the global nodes

     if ( strcmp(argv[2],"-modal")==0 ) { cout << "modal solution\n";
     ifstream in(argv[3],ios_base::in); // file of the modal solution
        
                    while( getline(in,line) ){
                    string word;
                    istringstream iss(line);
                    iss >> word;
                    if (word=="Eigenvalues:") break;
                    }

                    //getline(in,line);
                    //cout << line << '\n';
                    vector<double> freqs;
                           while(getline(in,line) && 
                                        ( isdigit(line[0]) || (line[0]=='-') )
                                 ){ // while we read numbers
                           ///cout << line << '\n';
                           istringstream iss(line);
                           string word; 
                           double d;
                           iss >> d; //cout << d <<'\n';
                           freqs.push_back(d); 
                           }

                    ///cout << line << '\n';
                    Matrix<double> Forms;
                    in >> Forms;
                    //Forms.print();
                    ////cout << Forms.size() << '\n';
                    in.close(); 
                    int nmodes = 0; 

                    // define additional specification for the output 
                    for(int i=4; i<argc;++i){
                    if (strcmp(argv[i],"-number")==0) nmodes = atoi(argv[i+1]);
                    if (strcmp(argv[i],"-Translat")==0) Translat = true;
                    }
                    


                    (nmodes > 0) ? nmodes = nmodes : nmodes = freqs.size();
                     cout << "number of modes: " << nmodes << '\n';


                    if(Translat) { 
                                  cout << "$\n";
                                  cout << "$ Rotor CMS\n";
                                  cout << "$    Tform Translat Punchfile\n";
                                  cout << setw(10) << "R" << " " << "/free-free_rotor.pch\n";
                                  cout << "$  Station      Grid         X         M        Ip        Id\n";
                                        for(int i=0;i<nodes.size();++i)
                                        cout << setw(10) << i + 1
                                             << setw(10) << nodes[i].station
                                             << setw(10) << setprecision(4) << fixed << nodes[i].x
                                             << setw(10) << setprecision(4) << fixed << nodes[i].Mass
                                             << setw(10) << setprecision(4) << fixed << nodes[i].Ip
                                             << setw(10) << setprecision(4) << fixed << nodes[i].Id
                                             << '\n';
                                  cout << "$\n";

                                    
                                      for (int m=1;m<=nmodes;++m){
                                      vector<double> mode = Forms.column(m);
                                      cout << "$TITLE = EnDy Free-Free Rotor\n";
                                      cout << "$SUBTITLE = from file " << filemesh << '\n';
                                      cout << "$LABEL\n";
                                      cout << "$EIGENVECTOR\n";
                                      cout << "$REAL OUTPUT\n";
                                      cout << "$SUBCASE ID\n";
                                      double eigval = (2*M_PI*freqs[m-1]) * (2*M_PI*freqs[m-1]);
                                      cout << "$EIGENVALUE =" << setw(15) << scientific << setprecision(8) << eigval
                                      << setw(8) << "MODE ="
                                      << setw(6) << m
                                      << '\n';
                                      
                                            for(int k=0;k<nodes.size();++k){
                                            int node = nodes[k].station;
                                            cout << setw(10) << node
                                                 <<  setw(8) << 'G'
                                                 <<  setw(18) << scientific << setprecision(8) << mode[6*node -6 + 1 - 1] // Ux
                                                 <<  setw(18) << scientific << setprecision(8) << mode[6*node -6 + 2 - 1] // Uy
                                                 <<  setw(18) << scientific << setprecision(8) << mode[6*node -6 + 3 - 1] // Uz
                                                 << '\n'
                                                 <<  setw(18) << "-CONT-            " // begining of the next line
                                                 <<  setw(18) << scientific << setprecision(8) << mode[6*node -6 + 4 - 1] // Rx
                                                 <<  setw(18) << scientific << setprecision(8) << mode[6*node -6 + 5 - 1] // Ry
                                                 <<  setw(18) << scientific << setprecision(8) << mode[6*node -6 + 6 - 1] // Rz
                                                 << '\n';
                                            }
                                      }

                                  return 0;
                    } 


                      
                     for(int i=0;i<nmodes;++i){
                     cout << "frequency: " << freqs[i] << " Hz\n";
                     vector<double> mode = Forms.column(i+1);
                        // cout << mode.size() << '\n';
                        //for(int i=0;i<mode.size();++i)
                        // cout << mode[i] <<'\n';
 
                          for(int j=0;j<nodes.size();++j){
                           int node = nodes[j].station;
                           if (node > 0)
                           cout << fixed << setw(WIDTH) << setprecision(PRECISION) << nodes[j].x 
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << (nodes[j].Ro + nodes[j].Ri)/2
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << mode[6*node -6 + 1 - 1]
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << mode[6*node -6 + 2 - 1]
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << mode[6*node -6 + 3 - 1]
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << mode[6*node -6 + 4 - 1]
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << mode[6*node -6 + 5 - 1]
                                << fixed << setw(WIDTH) << setprecision(PRECISION) << mode[6*node -6 + 6 - 1]
                                << '\n';
                          }
                     } // cycle over the 
                      
     }


    if ( strcmp(argv[2],"-transient")==0) cout << "transient solution\n"; 

    if ( strcmp(argv[2],"-static") == 0 ){   cout << "static solution\n";

                    ifstream in(argv[3],ios_base::in); // file of the static solution
        
                    while( getline(in,line) ){
                    string word;
                    istringstream iss(line);
                    iss >> word;
                    if (word=="Static:") break;
                    }
                    //cout << line << '\n';

                    vector<double> solution;
                            while(getline(in,line) && 
                                        ( isdigit(line[0]) || (line[0]=='-') )
                                 ){ // while we read numbers
                           ///cout << line << '\n';
                           istringstream iss(line);
                           string word; 
                           double d;
                           iss >> d; //cout << d <<'\n';
                           solution.push_back(d); 
                           }
                    //for(int i=0;i<solution.size();++i)
                    //cout << solution[i] <<'\n';
                    in.close();

                     /*
                    in.open(argv[4],ios_base::in); // read file of the stiffness matrix
                    Matrix<double> K;
                    in >> K;
                    in.close();
                    //// K.print();
                     */

                           cout << setw(WIDTH) <<  " X" 
                                << setw(WIDTH) <<  " R"
                                << setw(WIDTH) << " Ux"
                                << setw(WIDTH) << " Uy"
                                << setw(WIDTH) << " Uz"
                                << setw(WIDTH) << " Rotx"
                                << setw(WIDTH) << " Roty"
                                << setw(WIDTH) << " Rotz"
                                  
                                << setw(WIDTH) << " Qx"
                                << setw(WIDTH) << " Qy"
                                << setw(WIDTH) << " Qz"
                                << setw(WIDTH) << " Mx"
                                << setw(WIDTH) << " My"
                                << setw(WIDTH) << " Mz"
                                << '\n';

                    for(int element = 1; element < nodes.size();++element){
                    int node1 = nodes[element-1].station;
                    int node2 = nodes[element].station;
                    //// cout << node1 << " " << node2 << '\n';
                    Station st1 = nodes[element - 1];
                    Station st2 = nodes[element];
                         Matrix<double> ke = FORM_LOCAL_KE_BEAM(
                          st1.F(),
                          st1.JD(),
                          st1.JY(),
                          st1.JZ(),
                          st2.x - st1.x,
                          st1.E,
                          st1.G()
                         );

                    // Matrix<double> ke(12);
                    //              for(int i=1;i<=6;++i)
                    //               for(int j=1;j<=6;++j){
                    //                       ke(i,j) =     K(6*node1-6+i,6*node1-6+j);
                    //                   ke(i,j + 6) =     K(6*node1-6+i,6*node2-6+j);
                    //                   ke(i + 6,j) =     K(6*node2-6+i,6*node1-6+j);
                    //                   ke(i + 6,j + 6) = K(6*node2-6+i,6*node2-6+j);
                    //               }
                     // set_precision(ke, WRITE_PRECISION);
                    ///cout << st1.station << " " << st2.x - st1.x << '\n';
                    //ke.print();
                    vector<double> element_solution(12); // for the beam only
                    double L = st2.x - st1.x;
                    //element_solution[2] = 0;
                    //element_solution[4] = -1.0;
                    //element_solution[8] = L;
                    ///element_solution[10] = -1.0;
                    
                          for(int dof=1;dof<=6;++dof){
                    ///          cout << "node: " << node1 << " " << 6 * node1 - 6 + dof - 1 << " " << fixed << solution[6 * node1 - 6 + dof - 1] << '\n';
                              element_solution[dof - 1] = solution[6 * node1 - 6 + dof - 1]; // define solution for the node1
                              //cout << "node: " << node2 << " " << 6 * node2 - 6 + dof - 1 << '\n';
                              element_solution[dof - 1 + 6 ] = solution[6 * node2 - 6 + dof - 1]; // define solutin for the node2
                          } 
                      
                   // for(int i=0;i<solution.size();++i)
                   // cout << solution[i] <<'\n';
                          
                    //     for(int i = 0;i < element_solution.size() ;++i)
                   ///      cout << fixed << element_solution[i] << '\n';
                       
                    element_solution = ke * element_solution;
                    

                           cout <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    nodes[element-1].x  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    (nodes[element-1].Ro + nodes[element-1].Ri)/2 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node1 -6 + 1 - 1]  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node1 -6 + 2 - 1]  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node1 -6 + 3 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node1 -6 + 4 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node1 -6 + 5 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node1 -6 + 6 - 1] 
                                  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[0] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[2] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[3] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[4] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[5] 
                                << '\n';


                           cout <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    nodes[element].x 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    (nodes[element].Ro + nodes[element].Ri)/2 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node2 -6 + 1 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node2 -6 + 2 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node2 -6 + 3 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node2 -6 + 4 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node2 -6 + 5 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    solution[6*node2 -6 + 6 - 1] 
                                  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[0] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[2] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[3] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[4] + element_solution[2] * L 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION/2) <<    element_solution[5] - element_solution[1] * L 
                                << '\n';
                          

                    }//cycle over the elements                    
   }
  

return 0;
}
