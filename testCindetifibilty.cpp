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

vector<MatrixXd> all_c_identifible(int n)
    {

        vector<MatrixXd> identifiable;
        MatrixXd sigma(n,n);
        int cnt =0;
        for (int i =0;i<(1<<(n*n-n));i++)
        {
                    
            sigma=MatrixXd::Random(n,n).array().abs();
            MatrixXd active_vars=matrix_from_bitmask(i,n);
            MatrixXd solution=Generate_A_Matrix(active_vars,sigma,n,true);
            Eigen::ColPivHouseholderQR<MatrixXd> decomp(solution);
            if ((active_vars.rowwise().sum().minCoeff()>1) && (active_vars.colwise().sum().minCoeff()>1))
            if (((active_vars.sum()+n-1)==decomp.rank()) && ((active_vars.sum()==((n*n-n+2)/2))))
                {   
                    cout<<active_vars;
                    cout<<active_vars.rowwise().sum().minCoeff()<<" "<<active_vars.colwise().sum().minCoeff()<<endl;
                    identifiable.push_back(active_vars);
                    //cout<<active_vars<<endl;
                }
        }
        return identifiable;
    }

int main()
{
    int n =5;
    MatrixXd i=MatrixXd::Zero(n,n);
    MatrixXd B; 
    string m;
    VectorXd C;
    MatrixXd Cmat;

    //vector<MatrixXd> identifible= all_c_identifible(n);
    //cout<<identifible.size()<<endl;

    do {
        
        B=MatrixXd::Random(n,n).array().abs();
        C= VectorXd::Random(n,1).array().abs();


        Cmat= C.asDiagonal();

        //i=identifible[rand() % identifible.size()];
        i<<1,1,0,0,0,0,1,0,0,1,1,0,1,0,0,0,0,1,1,0,1,0,0,1,1;
        i.diagonal()=VectorXd::Zero(n);

        B=B.cwiseProduct(i);
        B.diagonal()-=B.rowwise().sum();
        B.diagonal()-=VectorXd::Random(n,1).array().abs().matrix();
    }
    while (abs(B.determinant())<0.1);

    cout<<B.determinant()<<endl;

    MatrixXd A = SolveLypunov(B,-Cmat,n);
    MatrixXd B_lasso = A.inverse().eval().array()*-0.5;
    VectorXd C_lasso = VectorXd::Random(n,1).array().abs();
    
    MatrixXd active_vars=((B.array().abs()>0).cast<double>());
    B_lasso=B_lasso.array()*active_vars.array();
    B_lasso.diagonal()-=B.rowwise().sum();


    cout<<"original:"<<endl;
    cout<<A<<endl;
    cout<<"first_guess:"<<endl;
     Cmat= C_lasso.asDiagonal();
    cout<<SolveLypunov(B_lasso,Cmat,n)<<endl;


    int converge_flag=0;
    tie(B_lasso,C_lasso,converge_flag)=FitLyapunov(A,B_lasso,C_lasso,n,0.5,0,0,0.001,active_vars,0.0000001,100,true,10000);

    cout<<converge_flag<<endl;

     B=B.array()/C[0];
    C=C.array()/C[0];

    B_lasso=B_lasso.array()/C_lasso[0];
    C_lasso=C_lasso.array()/C_lasso[0];

    cout<< B<< endl<<B_lasso<<endl;

    cout<<C<<endl<<C_lasso<<endl;
    

    cout<<A;
    MatrixXd AS= Generate_A_Matrix(active_vars,A,n,true);
    MatrixXd bvec= MatrixXd::Zero(active_vars.sum()+n-1,1);
    bvec(0,0)=1;

    Eigen::ColPivHouseholderQR<MatrixXd> decomp(AS);
    

    cout<<AS<<endl;
    cout<<"det:"<<AS.determinant()<<endl;
    cout<<"rank:"<<decomp.rank()<<endl;

    cout<<AS.inverse()<<endl;
    cout<<bvec<<endl;
    cout<<AS.inverse()*bvec;

}
