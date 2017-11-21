
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "Matrix.h"
#include "Modes.h"
#include "Read_Translat_Mesh.h"

int main(int argc, char *argv[]){

  if (argc<3) { cerr << " Stiffness eigenvectors Force file \n"; return 1; }

  //read Mass matrix eigenvectors
  // Matrix mass must be positive defiened
 
       ifstream ifs(argv[1]);

       if(!ifs) { cerr << " can not open file: " << argv[1] << '\n'; return 3; }
       vector<Mode> modes;

       cout << " read Stiffness vectors \n";
       ifs >> modes;
       ifs.close(); // compltete K modes reading


       //for(int i=0;i<modes.size();++i)
       //modes[i].print();
       //return 8;


       ifs.open(argv[2]);
       if(!ifs)  { cerr << " can not open file: " << argv[2] << '\n'; return 4; }
       vector<double> Force;
       string words;
       //cout << "Forces: " << argv[2] <<'\n';
              while(getline(ifs,words) && 
                                    !( isdigit(words[0]) || words[0] == '-' )
                    ){ cout << words << '\n'; }
              
              StringList list(words,';');
              Force.push_back( atof(list[0].c_str()) );
              
              while(getline(ifs,words) && 
                                    ( isdigit(words[0]) || words[0] == '-' )
                    ){
              list.set_string(words,';');
              ///cout << "inside: " << words << '\n';
              Force.push_back( atof(list[0].c_str()) );
              
               }
       //for(int i=0;i<Force.size();++i)
       //cout << Force[i] << '\n';

       ifs.close();

       Matrix<double> F(modes.size());
           for(int i=0;i<modes.size();++i)
               F.fill_column(modes[i].mode,i+1); 
       Matrix<double> FT = F;
       FT.transpose();      

       Matrix<double> invD(modes.size());
       /////////////cout << modes.size() << '\n';
           for(int i=0;i<modes.size();++i){
               assert(modes[i].frequency > 0);
               invD(i+1,i+1) = 1.0/ modes[i].frequency;
           }       


       //invD.print();
       //return 8;
       cout << " Static:\n";
       Force = F * invD * FT * Force;

       for(int i = 0; i< Force.size(); ++i)
       cout << fixed << setprecision(WRITE_PRECISION) << Force[i] << ";\n";

/*

      /// define vector for the global nodes

    vector<string> Trans;
           
    ifstream in(argv[3]);
    if (!in) { cerr << " can not open file: " << argv[1] << '\n'; }

    string word;
    while(in >> word) Trans.push_back(word);

      //for(int i=0;i<Trans.size();++i)
      //cout << Trans[i] << '\n';
     in.close();   



                    string filemesh, matrixfile;
                    cin >> filemesh >> matrixfile;
                    cout << "mesh: " << filemesh << " matrix: " << matrixfile << '\n';
                    
                    in.open(filemesh.c_str(),ios_base::in);
 
                     if(!in) throw runtime_error("can not open file: " + filemesh); 
                     vector<Station> nodes;
                     Station st;
                     string line;
                    
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

      /// define vector for the global nodes



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
                              element_solution[dof - 1] = Force[6 * node1 - 6 + dof - 1]; // define solution for the node1
                              //cout << "node: " << node2 << " " << 6 * node2 - 6 + dof - 1 << '\n';
                              element_solution[dof - 1 + 6 ] = Force[6 * node2 - 6 + dof - 1]; // define solutin for the node2
                          } 
                      
                   // for(int i=0;i<solution.size();++i)
                   // cout << solution[i] <<'\n';
                          
                    //     for(int i = 0;i < element_solution.size() ;++i)
                   ///      cout << fixed << element_solution[i] << '\n';
                       
                    element_solution = ke * element_solution;
                    

                           cout <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    nodes[element-1].x  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    (nodes[element-1].Ro + nodes[element-1].Ri)/2 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node1 -6 + 1 - 1]  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node1 -6 + 2 - 1]  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node1 -6 + 3 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node1 -6 + 4 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node1 -6 + 5 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node1 -6 + 6 - 1] 
                                  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[0] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[2] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[3] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[4] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[5] 
                                << '\n';


                           cout <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    nodes[element].x 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    (nodes[element].Ro + nodes[element].Ri)/2 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node2 -6 + 1 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node2 -6 + 2 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node2 -6 + 3 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node2 -6 + 4 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node2 -6 + 5 - 1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    Force[6*node2 -6 + 6 - 1] 
                                  
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[0] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[1] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[2] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[3] 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[4] + element_solution[2] * L 
                                <<   fixed << setw(WIDTH) << setprecision(PRECISION) <<    element_solution[5] - element_solution[1] * L 
                                << '\n';
                          

                    }//cycle over the elements
 *//////////////
    
return 0;
}
