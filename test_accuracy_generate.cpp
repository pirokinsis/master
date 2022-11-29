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
srand((unsigned int) 1000);
ofstream myfile("data.txt");
for (int k=0;k<1000;k++ )
{

    cout<<"RUN NUMBER: "<<k<<endl;
    int n=rand()%10+5;
    double p= static_cast <double> (rand()) / static_cast <double> (RAND_MAX);;
    //MatrixXd A = MatrixXd::Random(n,n);

    // for (int i=0;i<n;i++){
    //     for (int j=i+1;j<n;j++)
    //     {
    //         A(i,j)=A(j,i)+2;
    //     }
    // }

    MatrixXd B = MatrixXd::Random(n,n).array().abs();
    VectorXd C = VectorXd::Random(n,1).array().abs();
    MatrixXd Cmat= C.asDiagonal();
    MatrixXd i=MatrixXd::Random(n,n);


    do {
    i=(i.array().abs()>p).cast<double>();

    i.diagonal()=VectorXd::Zero(n);

    B=B.cwiseProduct(i);
    B.diagonal()-=B.rowwise().sum();
    B.diagonal()-=VectorXd::Random(n,1).array().abs().matrix();
    }
    while (B.determinant()==0);

    cout<<B.determinant()<<endl;

    MatrixXd A = SolveLypunov(B,-Cmat,n);

    MatrixXd B_lasso = A.inverse().eval().array()*-0.5;
    VectorXd C_lasso = VectorXd::Random(n,1).array().abs();

    MatrixXd B_group = A.inverse().eval().array()*-0.5;
    VectorXd C_group = VectorXd::Random(n,1).array().abs();


    //B3=B3.cwiseProduct(i);

    MatrixXd active_vars_lasso=MatrixXd::Constant(n,n,1);
    MatrixXd active_vars_group=MatrixXd::Constant(n,n,1);
    MatrixXd active_vars_prev;
    int converge_flag=0;

    double min_dif = n;
    double min_diff_accurracy;
    double min_pair = n;
    double min_pair_accurracy;
    double i_pairs = ((i.array() + i.transpose().array()) > 0).cast<double>().sum()/2.0;

    for (double j=0.0001;j<1;j*=pow(10000,0.01))
    {
        converge_flag=0;
        int count=0;
        while( (converge_flag==0) and (count<=100) )
        {
            count++;

            tie (B_lasso, C_lasso,converge_flag)=ComputeGradient(A,B_lasso,C_lasso,n,0.5,j,0,0.3,active_vars_lasso,0.0001,100);
        }

        if (active_vars_lasso.sum()==n) break;

        if(active_vars_prev.sum()!=active_vars_lasso.sum())
        {
            cout<<i<<endl;
            cout<<active_vars_lasso<<endl;
            cout<<"prediction accuracy:"<<accuracy(i,active_vars_lasso,n)<<endl;
            cout<<"prediction pairwise accuracy:"<<pairwise_accuracy(i,active_vars_lasso,n)<<endl;
            double act_pairs=(((active_vars_lasso.array()+active_vars_lasso.transpose().array())>0).cast<double>().sum()-n)/2.0;
            if ( abs(i.sum()-active_vars_lasso.sum()+n)<min_dif)
                {
                    min_dif=abs(i.sum()-active_vars_lasso.sum()+n);
                    min_diff_accurracy=accuracy(i,active_vars_lasso,n);
                }
        
            if(abs(i_pairs-act_pairs)<min_pair)
                {
                    min_pair=abs(i_pairs-act_pairs);
                    min_pair_accurracy=pairwise_accuracy(i,active_vars_lasso,n);
                }
                
        }
        active_vars_prev=active_vars_lasso;
    }
    double min_pair_g = n;
    double min_pair_accurracy_g;
    cout<<"Starting group lasso:"<<endl;
    for (double j=0.0001;j<1;j*=pow(10000,0.01))
    {
        converge_flag=0;
        int count=0;
        while( (converge_flag==0) and (count<=100) )
        {
            count++;

            tie (B_group, C_group,converge_flag)=ComputeGradient(A,B_group,C_group,n,0.5,0,j,0.3,active_vars_group,0.0001,100);
        }

        if (active_vars_group.sum()==n) break;

        if(active_vars_prev.sum()!=active_vars_group.sum())
        {
            cout<<i<<endl;
            cout<<active_vars_group<<endl;
            cout<<"prediction pairwise accuracy:"<<pairwise_accuracy(i,active_vars_group,n)<<endl;
            double act_pairs=(((active_vars_group.array()+active_vars_group.transpose().array())>0).cast<double>().sum()-n)/2.0;
            if(abs(i_pairs-act_pairs)<min_pair_g)
                {
                    min_pair_g=abs(i_pairs-act_pairs);
                    min_pair_accurracy_g=pairwise_accuracy(i,active_vars_group,n);
                }
                
        }
        active_vars_prev=active_vars_group;
    }
    cout<<"smallest diffrence: "<<min_dif<<" accuracy: "<<min_diff_accurracy<<endl;
    cout<<"smallest pairs diffrence: "<<min_pair<<" accuracy: "<<min_pair_accurracy<<endl;
    cout<<"GROUP smallest pairs diffrence: "<<min_pair_g<<" accuracy: "<<min_pair_accurracy_g<<endl;
    myfile<<min_diff_accurracy<<" "<<min_pair_accurracy<<" "<<min_pair_accurracy_g<<endl;

}
myfile.close();
}
