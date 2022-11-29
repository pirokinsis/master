#include "lyapunov_funcs.hpp"

double accuracy(MatrixXd A,MatrixXd B,int n)
    {
        MatrixXi cmp1 = (A.array()>0).cast<int>();
        MatrixXi cmp2 = (B.array()>0).cast<int>();
        cmp1.diagonal()=VectorXi::Zero(n).matrix();
        cmp2.diagonal()=VectorXi::Zero(n).matrix();

        return ((cmp1.array()*cmp2.array()).sum())/((double)cmp2.sum());
    }
double pairwise_accuracy(MatrixXd A,MatrixXd B,int n)
    {
        return accuracy(A.array().abs()+A.transpose().array().abs(),B.array().abs()+B.transpose().array().abs(),n);
    }

MatrixXd SoftThreshold(MatrixXd B, double lambda, int n, MatrixXd &active_vars)
{

        MatrixXd minus_l=(B.array().abs()-lambda);

        minus_l.diagonal()=B.diagonal().cwiseAbs();

        active_vars = (active_vars.array()*(minus_l.array()>0).cast<double>()).matrix();
        minus_l=minus_l.array()*active_vars.array()*B.array().sign();



        return minus_l;
        
}

tuple<double,MatrixXd> PairwiseLassoErrorSoftThreshold(MatrixXd B, double lambda, int n, MatrixXd &active_vars)
{
    MatrixXd B_squared=B.cwiseProduct(B).eval();


    B_squared+=B_squared.transpose().eval();

    B_squared=B_squared.cwiseSqrt();
    double result=(B_squared.sum()-B_squared.diagonal().sum())/2*lambda;
    B_squared+=((B_squared.array()==0).cast<double>()*Eigen::Infinity).matrix(); 

    B_squared=B.cwiseAbs()-lambda*B.cwiseAbs().cwiseQuotient(B_squared);
    B_squared.diagonal()=B.diagonal().cwiseAbs();
    
    active_vars=active_vars.array()*(B_squared.array()>0).cast<double>();

    B_squared=B_squared.array()*active_vars.array();
    B_squared=B_squared.array()*B.array().sign();

    return make_tuple(result,B_squared);
}



MatrixXd SoftThreshold(MatrixXd B, double lambda, int n)
{
        MatrixXd active= MatrixXd::Constant(n,n,1);
        return SoftThreshold(B, lambda,n,active);

}

Matrix2d Solve2dSyl(Matrix2d R1,Matrix2d R2,Matrix2d X)
{
    Matrix4d KroSum;

    KroSum.block<2,2>(0,0)=R1+R2(0,0)*Matrix2d::Identity(2,2);
    KroSum.block<2,2>(2,2)=R1+R2(1,1)*Matrix2d::Identity(2,2);
    KroSum.block<2,2>(2,0)=R2(1,0)*Matrix2d::Identity(2,2);
    KroSum.block<2,2>(0,2)=R2(0,1)*Matrix2d::Identity(2,2);
    Vector4d X_vec(Map<Vector4d>(X.data(), 4));
    Vector4d Y_vec=KroSum.inverse()*X_vec;
    return (Map<Matrix2d>(Y_vec.data()));

}

MatrixXd SolveSides(MatrixXd D,MatrixXd U, MatrixXd R_22, vector<int> starts,vector<int> sizes,int last)
{
    MatrixXd X=MatrixXd::Zero(starts[last],sizes[last]);
    for (int i=last-1;i>=0;i--)
    {

        MatrixXd Dj= D.block(starts[i],0,sizes[i],sizes[last]);

        for (int j=i+1;j<last;j++)
            Dj-= U.block(starts[i],starts[j],sizes[i],sizes[j])*X.block(starts[j],0,sizes[j],sizes[last]);

        if ((sizes[last]==2) && (sizes[i]==2)){

            X.block(starts[i],0,sizes[i],sizes[last])=
                Solve2dSyl(U.block(starts[i],starts[i],sizes[i],sizes[i]),R_22,Dj);
        }
        else if (sizes[i]==2)
        {

            X.block(starts[i],0,sizes[i],sizes[last])=
                (U.block(starts[i],starts[i],sizes[i],sizes[i])+R_22(0,0)*Matrix2d::Identity(2,2)).inverse()*Dj;
        }
        else if (sizes[last]==2)
        {
            X.block(starts[i],0,sizes[i],sizes[last])=
                ((U(starts[i],starts[i])*Matrix2d::Identity(2,2) + R_22).inverse()*Dj.transpose()).transpose();
        }

        else
        {

            X(starts[i],0)=Dj(0,0)/(U(starts[i],starts[i])+R_22(0,0));
        }

    }
    return X;

}


MatrixXd SolveLypunov (MatrixXd B,MatrixXd C,int n)
{
    vector <int> sizes;
    vector <int> starts;
    MatrixXd Y = MatrixXd::Zero(n,n);
    MatrixXd Y_22;

    RealSchur<MatrixXd> schur(B);
    MatrixXd U = schur.matrixU();

    MatrixXd oldC= U.transpose()*C*U;

    C=U.transpose()*C*U;
    MatrixXd T = schur.matrixT();

    for (int i=0;i<n;)
    {
        starts.push_back(i);
        if((i==n-1) || (T(i+1,i)==0)){
            sizes.push_back(1);
            i+=1;
        }
        else {
            sizes.push_back(2);
            i+=2;
        }
    }
    for (int i=starts.size()-1;i>=0;i--)
    {

        MatrixXd R_22= T.block(starts[i], starts[i], sizes[i], sizes[i]);
        MatrixXd R_12= T.block(0, starts[i], starts[i], sizes[i]);
        MatrixXd C_22=C.block(starts[i], starts[i], sizes[i], sizes[i]);
        MatrixXd C_12=C.block(0, starts[i], starts[i], sizes[i]);
        MatrixXd C_21=C.block(starts[i], 0, sizes[i],starts[i]);

        if (sizes[i]==2){
            Y_22= Solve2dSyl(R_22,R_22,C_22);
                }
        else {
            Y_22=C_22*(2*R_22).inverse();

            }
        MatrixXd C_12_hat=C_12-R_12*Y_22;
        MatrixXd C_21_hat=C_21.transpose()-R_12*Y_22.transpose();
        Y.block(starts[i],starts[i],sizes[i],sizes[i])=Y_22;

        if (i==0) break;

        MatrixXd Y_12= SolveSides(C_12_hat,T,R_22,starts,sizes,i);
        MatrixXd Y_21= SolveSides(C_21_hat,T,R_22,starts,sizes,i);

        Y.block(0,starts[i],starts[i],sizes[i])=Y_12;
        Y.block(starts[i],0,sizes[i],starts[i])=Y_21.transpose();



        C.block(0,0,starts[i],starts[i])-=R_12*Y_21.transpose()+Y_12*R_12.transpose();

    }
    return U*Y*U.transpose();
}


tuple<MatrixXd,VectorXd,int> ComputeGradient(MatrixXd Sigma, MatrixXd B, VectorXd C,int n, double alpha,double lambda,double lambda2, double kappa,MatrixXd &active_vars,double min_change, int max_steps)
{


    if (active_vars.cols()==0)
        active_vars=MatrixXd::Constant(n,n,1);

    //cout<<active_vars<<endl;
    MatrixXd Cmat= C.asDiagonal();
    MatrixXd S=SolveLypunov(B,-Cmat,n);


    double pairwise_error;
    MatrixXd pairwise_gradient;

    double f= (Sigma-S).norm();
    f=f*f;
    f+=lambda*(B.cwiseAbs().sum()-B.diagonal().cwiseAbs().sum());
    if (lambda2>0)
        {
            tie(pairwise_error,pairwise_gradient)=PairwiseLassoErrorSoftThreshold(B,lambda2,n,active_vars);
            f+=pairwise_error;
        }
    

    MatrixXd delta= (Sigma-S);
    MatrixXd D=SolveLypunov(B.transpose(),delta,n);


    MatrixXd grad_B= (2*S*D).cwiseProduct(active_vars);
    VectorXd grad_C=2*D.diagonal();//+2*kappa*(C-VectorXd::Constant(n, 1, 1.0));

//   cout<<grad_B<<endl;
//    cout<<grad_C<<endl;

    double step = 1;
    double step_B= 1;
    double step_C = 1;
    MatrixXd B_old=B;
    VectorXd C_old=C;

    bool bad_step=true;

    int count2=0;
    while(bad_step)
    {   

        count2++;
        B=B_old-step*step_B*grad_B;


        if (lambda>0)
            B=SoftThreshold(B,step*step_B*lambda,n,active_vars);
        if (lambda2>0)
            tie(pairwise_error,B)=PairwiseLassoErrorSoftThreshold(B,step*step_B*lambda2,n,active_vars);



        bad_step=false;
        RealSchur<MatrixXd> schur(B);
        MatrixXd T = schur.matrixT();

        for (int i=0;i<n;)
            {
                if ( (i==n-1) || (T(i+1,i)==0) )
                    {
                        if (T(i,i)>0) bad_step=true;
                        i+=1;
                    }
                else
                {
                        if (T(i,i)+T(i+1,i+1)>0) bad_step=true;
                        i+=2;
                }
            }
        step_B*=alpha;
    }
    C=0*C;
    while ((C.array()<=0).any())
        {
            C=C_old-step*step_C*grad_C;
            step_C*=alpha;
        }
    bad_step=true;
    double f_old=f;
    int counter=0;
    while (bad_step)
    {
        counter++;

        if (counter==max_steps){
            return make_tuple(B,C,1);
        }    
        bad_step=false;
        Cmat= C.asDiagonal();
        S=SolveLypunov(B,-Cmat,n);




        f= (Sigma-S).norm();
        f=f*f;
        f+=lambda*(B.cwiseAbs().sum()-B.diagonal().cwiseAbs().sum());
        if (lambda2>0)
            {
                tie(pairwise_error,pairwise_gradient)=PairwiseLassoErrorSoftThreshold(B,lambda2,n,active_vars);
                f+=pairwise_error;
            }
    
    
        

        if ((f_old-f)<min_change)
            {
                 step*= alpha; 
                 C=C_old-step*step_C*grad_C;
                 B=B_old-step*step_B*grad_B;
                 if (lambda>0)
                    B=SoftThreshold(B,step*step_B*lambda,n,active_vars);
                 if (lambda2>0)
                    tie(pairwise_error,B)=PairwiseLassoErrorSoftThreshold(B,step*step_B*lambda2,n,active_vars);
                 bad_step=true;
            }
        


    }
    // cout<<"f old:"<< f_old<<"   B norm"<<B.norm()<<"  C norm"<<C.norm()<<endl;
    // cout<< " S norm:" <<S.norm()<<"  C_grad_norm:"<<grad_C.norm()<<"  B_grad_norm"<<grad_B.norm()<<endl;
    // cout<<"f new:"<<f<<" number of steps "<<counter<<" step_b counter:"<<count2<<endl;
    return make_tuple(B,C,0);
}