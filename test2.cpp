#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

Matrix2d Solve2dSyl(Matrix2d R,Matrix2d X)
{
    Matrix4d KroSum;

    KroSum.block<2,2>(0,0)=R+R(0,0)*Matrix2d::Identity(2,2);
    KroSum.block<2,2>(2,2)=R+R(1,1)*Matrix2d::Identity(2,2);
    KroSum.block<2,2>(2,0)=R(1,0)*Matrix2d::Identity(2,2);
    KroSum.block<2,2>(0,2)=R(0,1)*Matrix2d::Identity(2,2);
    Vector4d X_vec(Map<Vector4d>(X.data(), 4));
    Vector4d Y_vec=KroSum.inverse()*X_vec;
    return (Map<Matrix2d>(Y_vec.data()));

}

int main()
{
Matrix2d a=Matrix2d::Random();
Matrix2d b=Matrix2d::Random();
Matrix2d Y=Solve2dSyl(a,b);
cout<<a<<endl<<b<<endl<< a*Y+Y*a.transpose();
}


