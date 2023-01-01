#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "lyapunov_funcs.hpp"

using namespace std;
using namespace Eigen;


int main()
{
Matrix3d a=Matrix3d::Random();
a=Matrix3d::Random();
MatrixXd b(3,3);
double lambda=0.2;
b<<1,1,1,
1,1,1,
1,1,1;


cout<<a<<endl;
cout<<a.cwiseAbs().sum()<<"|"<<a.diagonal().cwiseAbs().sum()<<"|"<<a.cwiseAbs().sum()-a.diagonal().cwiseAbs().sum()<<endl;
cout<<SoftThreshold(a,lambda,3,b)<<endl;
cout<<b<<endl;

a=Matrix3d::Random();

cout<<a<<endl;
cout<<SoftThreshold(a,lambda,3,b)<<endl;
cout<<b<<endl;


}


