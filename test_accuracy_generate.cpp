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
ofstream myfile_error("failed_to_converge.txt");
for (int k=0;k<10000;k++ )
{
    
    cout<<"RUN NUMBER: "<<k<<endl;

    int n;
    double p;
    MatrixXd i;
    MatrixXd B; 
    string m;
    VectorXd C;
    MatrixXd Cmat;
    do {
        n = rand() % 20 + 5;
        p= static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        
        B=MatrixXd::Random(n,n).array().abs();
        C= VectorXd::Random(n,1).array().abs();
        Cmat= C.asDiagonal();

        i=MatrixXd::Random(n,n);
        i=(i.array().abs()>p).cast<double>();

        i.diagonal()=VectorXd::Zero(n);

        B=B.cwiseProduct(i);
        B.diagonal()-=B.rowwise().sum();
        B.diagonal()-=VectorXd::Random(n,1).array().abs().matrix();
    }
    while (abs(B.determinant())<0.1);

    cout<<B.determinant()<<endl;

    MatrixXd A = SolveLypunov(B,-Cmat,n);

    bool group[2] ={false,true};
    bool method[2]={true,false};
    for (bool meth: method)
        for (bool gr : group)
            {

                MatrixXd B_lasso = A.inverse().eval().array()*-0.5;
                VectorXd C_lasso = VectorXd::Random(n,1).array().abs();
                MatrixXd active_vars_lasso=MatrixXd::Constant(n,n,1);
                MatrixXd active_vars_prev;
                int converge_flag=0;
                int tp1,fp1,tn1,fn1,tp2,fp2,tn2,fn2;

                tie(B_lasso,C_lasso,converge_flag)=ComputeGradient(A,B_lasso,C_lasso,n,0.5,0,0,0.001,active_vars_lasso,0.0001,100,meth);
                
                
                double max=B_lasso.cwiseAbs().sum();
                


                for (double j=0.0001;j<max;j*=pow(max*10000,0.01))
                {
                    active_vars_lasso=MatrixXd::Constant(n,n,1);
                    tie (B_lasso, C_lasso,converge_flag)=FitLyapunov(A,B_lasso,C_lasso,n,0.5,j*(!gr),j*gr,0.001,active_vars_lasso,0.001,1000,meth,100);
                    if (converge_flag==1)
                        { 
                            myfile_error<<k<<endl;
                            break;
                        }        
                    active_vars_lasso=(B_lasso.array()!=0.0).cast<double>();

                    tie(tp1,tn1,fp1,fn1)=accuracy(i,active_vars_lasso,n);
                    tie(tp2,tn2,fp2,fn2)=pairwise_accuracy(i,active_vars_lasso,n);

                    if (meth) m="varadno";
                    else m ="fitch";

                    if (gr)
                        {
                         myfile<<k<<";"<<p<<";"<<n<<";"<<"group;"<<m<<";"<<tp2<<";"<<tn2<<";"<<fp2<<";"<<fn2<<endl;
                        }
                    else
                        {
                        myfile<<k<<";"<<p<<";"<<n<<";"<<"normal;"<<m<<";"<<tp1<<";"<<tn1<<";"<<fp1<<";"<<fn1<<endl;
                        myfile<<k<<";"<<p<<";"<<n<<";"<<"normalg;"<<m<<";"<<tp2<<";"<<tn2<<";"<<fp2<<";"<<fn2<<endl;

                        }
                    if (active_vars_lasso.sum()==n) break;
                        if (active_vars_lasso.sum()!=active_vars_prev.sum())
                        {
                        for (int k=0;k<n;k++)
                            {

                                cout<<i.row(k)<<"   ";
                                cout<<active_vars_lasso.row(k)<<endl;
                            }
                        cout<<"TP: "<<tp1<<" TN: "<<tn1<<" FP: "<<fp1<<" FN: "<<fn1<<endl;
                        cout<<"TP: "<<tp2<<" TN: "<<tn2<<" FP: "<<fp2<<" FN: "<<fn2<<endl;
                        }
                    active_vars_prev=active_vars_lasso;
                }
            }

}
myfile.close();
}
