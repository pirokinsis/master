#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "lyapunov_funcs.hpp"
#include <vector>
#include <fstream>
#include <tuple>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix2d;
using Eigen::Matrix4d;
using Eigen::Vector4d;
using Eigen::Map;
using Eigen::RealSchur;


int main()
{

int n=5;

//MatrixXd A = MatrixXd::Random(n,n);

// for (int i=0;i<n;i++){
//     for (int j=i+1;j<n;j++)
//     {
//         A(i,j)=A(j,i)+2;
//     }
// }

MatrixXd B = MatrixXd::Random(n,n);
VectorXd C = VectorXd::Random(n,1).array().abs();


for(int i=0;i<n;i++) B(i,i)-=1;

MatrixXd i(5,5);

/* i << 1, 0,0,1, 0,
     0,1,1, 0, 0,
     1,0,1, 0, 0,
      0,1,0, 0, 1,
       1,0,0, 0, 1;*/
i=MatrixXd::Constant(n,n,1);

MatrixXd Cmat= C.asDiagonal();

B=B.cwiseProduct(i);
MatrixXd A = SolveLypunov(B,Cmat,n);


MatrixXd B3 = MatrixXd::Random(n,n);
VectorXd C3 = VectorXd::Random(n,1).array().abs();
for(int i=0;i<n;i++) B3(i,i)-=3;

MatrixXd B2;
VectorXd C2;



RealSchur<MatrixXd> schur(B);
MatrixXd T = schur.matrixT();
//cout<<T;



B3=B3.cwiseProduct(i);
cout<<A<<endl;



int converge_flag=0;
for (double j=0.0001;j<10;j*=pow(100000,0.01))
{
    converge_flag=0;
    int count=0;
    while( (converge_flag==0) )
    {
        count++;

        tie (B2, C2,converge_flag)=ComputeGradient(A,B3,C3,n,0.5,0,j,0.3,i,0.0001,100);
        B3=B2;C3=C2;
    }
    cout<<j<<endl;
    cout<<B2.norm()<<endl;
    Cmat=C2.asDiagonal();
    cout<< (SolveLypunov(B2,Cmat,n)-A).norm()<<endl;
    cout<<"number of steps: "<<count<<endl;
    cout<<i<<endl;
}
//  cout<<B<<endl<<C<<endl;
//  cout<<"SOLUTION:_______________________;"<<endl;
//  cout<<B3<<endl<<C3<<endl<<endl;
//cout<<i<<endl;
//  ofstream myfile;
//   myfile.open ("example.txt");
//   myfile <<n<<endl;
//   myfile<<B3;
//   myfile.close();
}
