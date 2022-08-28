#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

int main()
{
    Vector3d X;
    X<< 1, 2, 3;
    Matrix3d m;
    m << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;
    VectorXd y= VectorXd::Constant(10, 1, 1.0);
    cout<<VectorXd::Constant(10, 1, 1.0)<<endl;

    X=0*X;
    cout<<(X.array() != 0.0).any()<<endl;
}
