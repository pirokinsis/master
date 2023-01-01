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
a=Matrix3d::Random();
MatrixXd b(3,3);
double lambda=0.6;
b<<1,1,1,
1,1,1,
1,1,1;
cout<<a<<endl;
MatrixXd res_mat;
double res_num;
tie(res_num,res_mat)=PairwiseLassoErrorSoftThreshold(a,lambda,3,b);

cout<<res_num<<endl<<res_mat<<endl;
cout<<b<<endl;


a=Matrix3d::Random();
tie(res_num,res_mat)=PairwiseLassoErrorSoftThreshold(a,lambda,3,b);

cout<<res_num<<endl<<res_mat<<endl;
cout<<b<<endl;

}


