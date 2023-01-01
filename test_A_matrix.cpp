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
 int n=4;
 MatrixXd active_vars(n,n);
 active_vars<<1,0,0,0,
 1,1,1,1,
 1,1,1,1,
 0,0,0,1;
 MatrixXd sigma(n,n);
 sigma=MatrixXd::Random(n,n).array().abs();

 sigma=(sigma+sigma.transpose()).eval();

MatrixXd solution=Generate_A_Matrix(active_vars,sigma,n).transpose();

Eigen::ColPivHouseholderQR<MatrixXd> decomp(solution);

cout<<active_vars.sum()<<" "<<decomp.rank()<<endl;


}