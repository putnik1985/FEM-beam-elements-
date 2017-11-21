#ifndef MATRIX_H
#define MATRIX_H

#include "stl.h"
#include "Transform.h"

#define PRECISION 4
#define WIDTH 16

#define WRITE_PRECISION 10

template<typename T>
class Matrix{ // numberig starts from the 1
vector<T> v;
int dimension;

      void sin_cos(
                   double& x,
                   double& y, 
                   double& s,
                   double& c
                   )
      {

      int k = 1;
      while ( pow(4,-k)*(x*x+y*y) > 1.0 )
      ++k;

      x/= pow(2,k);
      y/= pow(2,k);

      
      c = sqrt( 0.5 * (1 + y/sqrt(x*x+y*y)) );
      s = 0.5 * ( x / sqrt(x*x+y*y) ) / c;    
      }


public:
        void balance(const double& tol = 0.000000001){ 
                      for(int i = 1; i <= dimension; ++i)
                          for(int j = i + 1; j <= dimension; ++j)
                      if (abs( get(i,j) ) < tol ) 
                         { get(i,j) = 0.0; get(j,i) = 0.0; }
        }
                      
        void clear() { v.clear(); dimension = 0;}
	Matrix(int dim = 2):
        dimension(dim)
        {v.resize(dim*dim);}

	const T& operator()(
                            const int& i, 
                            const int& j
                           ) const { 
           return v[ (i - 1) * dimension + (j-1)]; 
        } // must be sure that i <= dimension && j << dimension

	T& operator()(
                            const int& i, 
                            const int& j
                     ) { 
           return v[ (i-1) * dimension + (j-1)]; 
        } // must be sure that i <= dimension && j << dimension

        int size() const { return dimension; }

        void set_dimension(const int& n) { 

             vector<T> new_v(n*n);
             //cout << " new dimension " << n << '\n';
             int dim = min(dimension,n);

             for(int i = 1; i <= dim; ++i)
             for(int j = 1; j <= dim; ++j)
             new_v[ (i - 1) * n + (j-1)] = get(i,j);

             v.clear(); // clear current one
             dimension = n;
             v = new_v;
            
            }

      void print() const{
      for(int i=1;i<=dimension;++i){
      for(int j=1;j<=dimension;++j)   
      cout << fixed << setprecision(WRITE_PRECISION) << this->operator()(i,j) << ";";
      cout << '\n';
      }
     } 


      void print_diag(){
      for(int i=1;i<=dimension;++i){
      cout << fixed << setprecision(WRITE_PRECISION) << this->operator()(i,i) << ";";
      cout << '\n';
      }
     } 

     Matrix<T>& operator=(const Matrix<T>& m){
     //Matrix<T> *p = &m;
     if(m.size()!=dimension) 
        set_dimension(m.size());
        for(int i=1;i<=dimension;++i)
         for(int j=1;j<=dimension;++j)
             get(i,j) = m(i,j);

     return *this;
     }

     T maximum_out_diag(int& i1, int& j1){
     i1 = 1; j1 = 2;
     T d = abs(get(1,2)); // first out diagonal element
       for(int i=1;i<=dimension;++i)
        for(int j=i+1;j<=dimension;++j)
         if ( d < 
                  abs( get(i,j) )
             ){
               d = abs( get(i,j) ); 
               i1 = i;
               j1 = j; 
              }
     return d; 
     }

     T& get(const int& i, const int& j) { return this->operator()(i,j); } // must check validation of the indexes
     const T& get(const int& i, const int& j) const { return this->operator()(i,j); } // must check validation of the indexes

     Matrix<T> apply_Jacobi_rotations(const double& tolerance = 0.000001){
     Matrix<T> EVectors(dimension);
     EVectors.unit();

     int p = 0, r = 0;
     double teta =0.0;

     T d = maximum_out_diag(p,r);
         while (d > tolerance ) {
    
         if ( get(r,r) != get(p,p) ) 
         teta = 0.5*atan(2*get(p,r)/ (get(r,r) - get(p,p)));
         else{
         if ( get(p,r) > 0) teta = M_PI / 4.0 ;
         if ( get(p,r) < 0) teta = - M_PI / 4.0 ;
         if ( get(p,r) == 0.0) teta = 0.0 ;
         }
         if (teta > M_PI/4.0) teta = teta - M_PI / 2.0;
         if (teta < -M_PI/4.0) teta = teta + M_PI / 2.0;
                              

        double  x = abs( 2 * get(p,r) );
        if (teta < 0.0) x = -x;
        
        double y = abs(get(r,r) - get(p,p));
    
        double cos, sin; 
        sin_cos(x, y, sin, cos);

        //cout << " sin " << sin << '\n';
        //cout << " cos " << cos << '\n';

        double app = get(p,p);
        double arr = get(r,r);       
        double apr = get(p,r);

        
        get(p,p) = app*cos*cos - 2*cos*sin*apr + arr*sin*sin ;
        get(r,r) = app*sin*sin + 2*cos*sin*apr + arr*cos*cos ;

        
        get(p,r) = sin*cos*(app-arr)+apr*(cos*cos-sin*sin);
        get(r,p) = sin*cos*(app-arr)+apr*(cos*cos-sin*sin);
        

        //cout << " apr " << get(p,r) << '\n';
        //cout << " arp " << get(r,p) << '\n';


        //get(p,r) = 0.0 ;
        //get(r,p) = get(p,r) ;

       
        for(int i=1; i <=dimension; ++i){

            if ( (i != p) && ( i != r)){

                                   double aip = get(i,p);
                                   double air = get(i,r);

                                   get(i,p) = cos*aip - sin*air;
                                   get(i,r) = sin*aip + cos*air;

                                   get(p,i) = get(i,p);
                                   get(r,i) = get(i,r);
                                   } 
            
        double aip = EVectors(i,p);
        double air = EVectors(i,r);


        EVectors(i,p) = cos*aip - sin*air;
        EVectors(i,r) = sin*aip + cos*air;
          } // process other elements

         //$%^%$^$%^$%^$%^$%^$%
         //print();
         //#%#$%#$%#$%#$%#$%#$%#$%

           d = maximum_out_diag(p,r);
         //cout << d <<'\n';
         //cout << p << " ; " << r << '\n';
         //char c;
         //cout << " continue??? \n";
         //cin >> c;

         } // while in tolerance      
     return EVectors;
     }
  

     void fill_row(const vector<T>& v, int i) {// must be sure v.size == dimension  
                                               for(int k=1;k<=dimension;++k)
                                               get(i,k) = v[k-1];
                                              }
     void fill_column(const vector<T>& v, int i) {  // must be sure v.size == dimension 
                                               for(int k=1;k<=dimension;++k)
                                               get(k,i) = v[k-1];
                                              }

     vector<T> column(const int i) {
     vector<T> v(dimension);
       for(int k=1;k<=dimension;++k){
        v[k-1] = get(k,i);
        // cout << "column: " << v[k-1] << '\n'; 
        }
     return v;
     }


     vector<T> row(const int i) {
     vector<T> v(dimension);
       for(int k=1;k<=dimension;++k)
        v[k-1] = get(i,k);
     return v;
     }


     void random(){
     for(int i=1;i<=dimension;++i)
      for(int j=i;j<=dimension;++j){
          get(i,j) = rand();
          get(j,i) = get(i,j);
      }
     }//random filling

     void unit(){
     for(int i=1;i<=dimension;++i)
         get(i,i) = 1.0;
   
     }// unit filling


     bool symmetric() const{
      for(int i=1;i<=dimension;++i)
          for(int j=i+1;j<=dimension;++j){
          //cout << get(i,j) << " " << get(j,i) << '\n';
              if (get(i,j)!=get(j,i)) {cout << " i= " << i
                                            << " j= " << j 
                                            << " delta = " << get(i,j) - get(j,i) << '\n';
                                        return false;
                                      }
          }
     return true;
     }


     void transpose(){
     for(int i=1;i<=dimension;++i)
      for(int j=i+1;j<=dimension;++j)
      swap(get(i,j), get(j,i));
     }       

};// Matrix class 

template<typename T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B){

     if( A.size() == B.size() ){
     Matrix<T> C(A.size());
     int dim = A.size();
     for(int i=1;i<=dim;++i)
     for(int j=1;j<=dim;++j)
     C(i,j) = A(i,j) + B(i,j);
     return C;
     }
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B){

     if( A.size() == B.size() ){
     Matrix<T> C(A.size());
     int dim = A.size();
     for(int i=1;i<=dim;++i)
     for(int j=1;j<=dim;++j)
     C(i,j) = A(i,j) - B(i,j);
     return C;
     }
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B){

     if( A.size() == B.size() ){
     int dim = A.size();
     Matrix<T> C(dim);

     for(int i=1;i<=dim;++i)
     for(int j=1;j<=dim;++j)
       for(int k=1;k<=dim;++k)
       C(i,j) = C(i,j) + A(i,k) * B(k,j);
     return C;
     }
}


template<typename T>
Matrix<T> operator*(const Matrix<T>& A, const T& b){

     int dim = A.size();
     Matrix<T> C(dim);

     for(int i=1;i<=dim;++i)
     for(int j=1;j<=dim;++j)
       C(i,j) = A(i,j) * b;
     return C;
}

template<typename T>
Matrix<T> operator*(const T& b, const Matrix<T>& A){

     int dim = A.size();
     Matrix<T> C(dim);
     
     for(int i=1;i<=dim;++i)
     for(int j=1;j<=dim;++j)
        C(i,j) = A(i,j) * b;
     return C;
}

template<typename T>
ostream& operator<<(ostream& os, const Matrix<T>& M){
      int dimension = M.size();
      for(int i=1;i<=dimension;++i){
      for(int j=1;j<=dimension;++j)   
      os << M(i,j) << ";";
      os << '\n';
      }
return os;
}


template<typename T>
istream& operator>>(istream& is, Matrix<T>& M){
string line;      

      while(getline(is,line) && !isdigit(line[0]) && line[0]!='-');

      StringList list(line,';'); //assumeed that numbers are separated by ;
      //list.print(); cout << '\n';

      int dimension = list.size();
      M.set_dimension(dimension);
      int k =1;
      M.fill_row(list.numbers(),k); // T must be a double
      

     
      while( getline(is,line) && 
               ( 
                 isdigit(line[0]) || 
                 line[0]=='-' 
               )
            ){ 
      ++k;
      list.set_string(line,';');
      M.fill_row(list.numbers(),k); // T must be a double
      }
      
return is;
}


template<typename T>
vector<T> operator*(const Matrix<T>& A, const vector<T>& v){

     int dim = A.size();
     vector<T> c;
     ////cout << "dimension: " << dim << '\n';
     for(int i=1;i<=dim;++i){
     double s =0;

         for(int j=1;j<=dim;++j)
           s+= A(i,j) * v[j-1];

     c.push_back(s);  
    }
     return c;
}


template<typename T>
string to_str(const T& a){
stringstream ss;
ss << a;
return ss.str();
}

#endif
