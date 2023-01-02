#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
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

MatrixXd Generate_A_Matrix(MatrixXd active_vars, MatrixXd Sigma, int n,bool diagonal_c=false)
{
    MatrixXd output;
    if (!diagonal_c)  output= MatrixXd::Zero(n * (n + 1) / 2, active_vars.sum());
    else output = MatrixXd::Zero(n * (n + 1) / 2, (int)active_vars.sum()+n-1);
    int col = 0;
    int row = 0;
    for (int k = 0; k < n; k++)
        for (int l = k; l < n; l++)
        {
            col = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    if (active_vars(j,i) != 0)
                    {
                        if ((j == l) && (j == k))
                            output(row, col) = 2 * Sigma(j, i);
                        else if (j == l)
                            output(row, col) = Sigma(k, i);
                        else if (j == k)
                            output(row, col) = Sigma(l, i);
                        col++;
                    }
            row++;
        }
    if (diagonal_c)
    {
        row=0;
        col=0;
        for (int k = 0; k < n; k++)
            for (int l = k; l < n; l++)
            {
                if ((k==l) && (k!=0))
                {

                    output(row,(int)active_vars.sum()+col)=-1;
                    col++;
                }
                row++;
            }
        return output;
    }
}

MatrixXd matrix_from_bitmask(unsigned bitmask,int n)
{
    MatrixXd output = MatrixXd::Identity(n,n);

    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
        {
            if (i!=j)
            {
                if (bitmask&1==1)output(i,j)=1;
                bitmask=bitmask>>1;
            }
        }
    return output;
}


int main()
{
int n =4;
MatrixXd sigma(n,n);
sigma=MatrixXd::Random(n,n).array().abs();
int cnt =0;
for (int i =0;i<(1<<(n*(n-1)));i++)
{
    cout<<i<<endl;
    MatrixXd active_vars=matrix_from_bitmask(i,n);
    MatrixXd solution=Generate_A_Matrix(active_vars,sigma,n,true);
    Eigen::ColPivHouseholderQR<MatrixXd> decomp(solution);
    if (active_vars.sum()==decomp.rank())cnt++;
}
cout<<cnt<<"/"<<(1<<(n*(n-1)))<<endl;
// int n=5;
// MatrixXd active_vars(n,n);
// active_vars<<
// 1,0,1,
// 1,1,0,
// 0,1,1,1,1,
// 1,1,1,1,0;

//
//
//MatrixXd solution=Generate_A_Matrix(active_vars,sigma,n);




}
