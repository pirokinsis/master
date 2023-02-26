#include "lyapunov_funcs.hpp"

tuple<int,int,int,int> accuracy(MatrixXd Y, MatrixXd X, int n)
{
    MatrixXi cmp1 = (Y.array() != 0).cast<int>();
    MatrixXi cmp2 = (X.array() != 0).cast<int>();
    cmp1.diagonal() = VectorXi::Zero(n).matrix();
    cmp2.diagonal() = VectorXi::Zero(n).matrix();

    int tp = (((cmp1.array() == cmp2.array()).cast<int>()) * cmp2.array()).sum();
    int tn = (((cmp1.array() == cmp2.array()).cast<int>()) * ((cmp2.array()) == 0).cast<int>()).sum();
    int fp = (((cmp1.array() != cmp2.array()).cast<int>()) * cmp2.array()).sum();
    int fn = (((cmp1.array() != cmp2.array()).cast<int>()) * ((cmp2.array()) == 0).cast<int>()).sum();

    return make_tuple(tp,(tn-n),fp,fn);
}
tuple<int,int,int,int>  pairwise_accuracy(MatrixXd A, MatrixXd B, int n)
{
    return accuracy(A.array().abs() + A.transpose().array().abs(), B.array().abs() + B.transpose().array().abs(), n);
}

MatrixXd SoftThreshold(MatrixXd B, double lambda, int n, MatrixXd active_vars)
{

    MatrixXd minus_l = (B.array().abs() - lambda);

    minus_l.diagonal() = B.diagonal().cwiseAbs();

    active_vars = active_vars.array() * (minus_l.array() > 0).cast<double>();
    minus_l = minus_l.array() * B.array().sign() * active_vars.array();

    return minus_l;
}

tuple<double, MatrixXd> PairwiseLassoErrorSoftThreshold(MatrixXd B, double lambda, int n, MatrixXd active_vars)
{
    MatrixXd B_squared = B.cwiseProduct(B).eval();

    B_squared += B_squared.transpose().eval();

    B_squared = B_squared.cwiseSqrt();
    double result = (B_squared.sum() - B_squared.diagonal().sum()) / 2 * lambda;
    B_squared += ((B_squared.array() == 0).cast<double>() * Eigen::Infinity).matrix();

    B_squared = B.cwiseAbs() - lambda * B.cwiseAbs().cwiseQuotient(B_squared);
    B_squared.diagonal() = B.diagonal().cwiseAbs();

    active_vars = active_vars.array() * (B_squared.array() > 0).cast<double>();

    B_squared = B_squared.array() * active_vars.array();
    B_squared = B_squared.array() * B.array().sign();

    return make_tuple(result, B_squared);
}

MatrixXd SoftThreshold(MatrixXd B, double lambda, int n)
{
    MatrixXd active = MatrixXd::Constant(n, n, 1);
    return SoftThreshold(B, lambda, n, active);
}

Matrix2d Solve2dSyl(Matrix2d R1, Matrix2d R2, Matrix2d X)
{
    Matrix4d KroSum;

    KroSum.block<2, 2>(0, 0) = R1 + R2(0, 0) * Matrix2d::Identity(2, 2);
    KroSum.block<2, 2>(2, 2) = R1 + R2(1, 1) * Matrix2d::Identity(2, 2);
    KroSum.block<2, 2>(2, 0) = R2(1, 0) * Matrix2d::Identity(2, 2);
    KroSum.block<2, 2>(0, 2) = R2(0, 1) * Matrix2d::Identity(2, 2);
    Vector4d X_vec(Map<Vector4d>(X.data(), 4));
    Vector4d Y_vec = KroSum.inverse() * X_vec;
    return (Map<Matrix2d>(Y_vec.data()));
}

MatrixXd SolveSides(MatrixXd D, MatrixXd U, MatrixXd R_22, vector<int> starts, vector<int> sizes, int last)
{
    MatrixXd X = MatrixXd::Zero(starts[last], sizes[last]);
    for (int i = last - 1; i >= 0; i--)
    {

        MatrixXd Dj = D.block(starts[i], 0, sizes[i], sizes[last]);

        for (int j = i + 1; j < last; j++)
            Dj -= U.block(starts[i], starts[j], sizes[i], sizes[j]) * X.block(starts[j], 0, sizes[j], sizes[last]);

        if ((sizes[last] == 2) && (sizes[i] == 2))
        {

            X.block(starts[i], 0, sizes[i], sizes[last]) =
                Solve2dSyl(U.block(starts[i], starts[i], sizes[i], sizes[i]), R_22, Dj);
        }
        else if (sizes[i] == 2)
        {

            X.block(starts[i], 0, sizes[i], sizes[last]) =
                (U.block(starts[i], starts[i], sizes[i], sizes[i]) + R_22(0, 0) * Matrix2d::Identity(2, 2)).inverse() * Dj;
        }
        else if (sizes[last] == 2)
        {
            X.block(starts[i], 0, sizes[i], sizes[last]) =
                ((U(starts[i], starts[i]) * Matrix2d::Identity(2, 2) + R_22).inverse() * Dj.transpose()).transpose();
        }

        else
        {

            X(starts[i], 0) = Dj(0, 0) / (U(starts[i], starts[i]) + R_22(0, 0));
        }
    }
    return X;
}

MatrixXd SolveLypunov(MatrixXd B, MatrixXd C, int n)
{
    vector<int> sizes;
    vector<int> starts;
    MatrixXd Y = MatrixXd::Zero(n, n);
    MatrixXd Y_22;

    RealSchur<MatrixXd> schur(B);
    MatrixXd U = schur.matrixU();

    MatrixXd oldC = U.transpose() * C * U;

    C = U.transpose() * C * U;
    MatrixXd T = schur.matrixT();

    for (int i = 0; i < n;)
    {
        starts.push_back(i);
        if ((i == n - 1) || (T(i + 1, i) == 0))
        {
            sizes.push_back(1);
            i += 1;
        }
        else
        {
            sizes.push_back(2);
            i += 2;
        }
    }
    for (int i = starts.size() - 1; i >= 0; i--)
    {

        MatrixXd R_22 = T.block(starts[i], starts[i], sizes[i], sizes[i]);
        MatrixXd R_12 = T.block(0, starts[i], starts[i], sizes[i]);
        MatrixXd C_22 = C.block(starts[i], starts[i], sizes[i], sizes[i]);
        MatrixXd C_12 = C.block(0, starts[i], starts[i], sizes[i]);
        MatrixXd C_21 = C.block(starts[i], 0, sizes[i], starts[i]);

        if (sizes[i] == 2)
        {
            Y_22 = Solve2dSyl(R_22, R_22, C_22);
        }
        else
        {
            Y_22 = C_22 * (2 * R_22).inverse();
        }
        MatrixXd C_12_hat = C_12 - R_12 * Y_22;
        MatrixXd C_21_hat = C_21.transpose() - R_12 * Y_22.transpose();
        Y.block(starts[i], starts[i], sizes[i], sizes[i]) = Y_22;

        if (i == 0)
            break;

        MatrixXd Y_12 = SolveSides(C_12_hat, T, R_22, starts, sizes, i);
        MatrixXd Y_21 = SolveSides(C_21_hat, T, R_22, starts, sizes, i);

        Y.block(0, starts[i], starts[i], sizes[i]) = Y_12;
        Y.block(starts[i], 0, sizes[i], starts[i]) = Y_21.transpose();

        C.block(0, 0, starts[i], starts[i]) -= R_12 * Y_21.transpose() + Y_12 * R_12.transpose();
    }
    return U * Y * U.transpose();
}

bool check_instablity(MatrixXd B, int n)
{
    RealSchur<MatrixXd> schur(B);
    MatrixXd T = schur.matrixT();

    bool bad_step = false;

    for (int i = 0; i < n;)
    {
        if ((i == n - 1) || (T(i + 1, i) == 0))
        {
            if (T(i, i) > 0)
                bad_step = true;
            i += 1;
        }
        else
        {
            if (T(i, i) + T(i + 1, i + 1) > 0)
                bad_step = true;
            i += 2;
        }
    }

    return bad_step;
}

tuple<double, MatrixXd, VectorXd> ForbeniusLossGradient(MatrixXd Sigma, MatrixXd B, MatrixXd C, int n, double lambda, double lambda2, double kappa, MatrixXd active_vars)
{
    MatrixXd Cmat = C.asDiagonal();
    MatrixXd S = SolveLypunov(B, -Cmat, n);
    MatrixXd delta = (Sigma - S);
    MatrixXd D = SolveLypunov(B.transpose(), delta, n);

    double f = (Sigma - S).norm();
    f = f * f;
    double pairwise_error;
    MatrixXd pairwise_gradient;

    f += lambda * (B.cwiseAbs().sum() - B.diagonal().cwiseAbs().sum());
    if (lambda2 > 0)
    {
        tie(pairwise_error, pairwise_gradient) = PairwiseLassoErrorSoftThreshold(B, lambda2, n, active_vars);
        f += pairwise_error;
    }

    f += kappa * (C - VectorXd::Constant(n, 1, 1.0)).squaredNorm();

    MatrixXd grad_B = (2 * S * D).cwiseProduct(active_vars);
    VectorXd grad_C = 2 * D.diagonal() + 2 * kappa * (C - VectorXd::Constant(n, 1, 1.0));

    return make_tuple(f, grad_B, grad_C);
}

tuple<double, MatrixXd, VectorXd> LassoGradient(MatrixXd Sigma, MatrixXd B, VectorXd C, int n, double lambda, double lambda2, double kappa, MatrixXd active_vars)
{
    MatrixXd Cmat = C.asDiagonal();
    double f = (B * Sigma + Sigma * B.transpose() + Cmat).squaredNorm();
    double pairwise_error;
    MatrixXd pairwise_gradient;
    f += lambda * (B.cwiseAbs().sum() - B.diagonal().cwiseAbs().sum());

    if (lambda2 > 0)
    {

        tie(pairwise_error, pairwise_gradient) = PairwiseLassoErrorSoftThreshold(B, lambda2, n, active_vars);

        f += pairwise_error;
    }

    MatrixXd grad_B = 2 * ((B * Sigma + Sigma * B.transpose()) + Cmat ) * Sigma;
    VectorXd grad_C =   ((B * Sigma + Sigma * B.transpose()).diagonal() + C) + 2 * kappa * (C - VectorXd::Constant(n, 1, 1.0));

    return make_tuple(f, grad_B, grad_C);
}

tuple<MatrixXd, VectorXd, int> ComputeGradient(MatrixXd Sigma, MatrixXd B, VectorXd C, int n, double alpha, double lambda, double lambda2, double kappa, MatrixXd active_vars, double min_change, int max_steps, bool fitch_or_verando)
{

    if (active_vars.cols() == 0)
        active_vars = MatrixXd::Constant(n, n, 1);



    double f;
    double pairwise_error;
    MatrixXd grad_B;
    VectorXd grad_C;
    MatrixXd temp1;
    VectorXd temp2;

    if (fitch_or_verando)
        tie(f, grad_B, grad_C) = ForbeniusLossGradient(Sigma, B, C, n, lambda, lambda2, kappa, active_vars);
    else
        tie(f, grad_B, grad_C) = LassoGradient(Sigma, B, C, n, lambda, lambda2, kappa, active_vars);
    
    grad_B=grad_B.cwiseProduct(active_vars);
    
    double step = 1;
    if (!fitch_or_verando) step=0.01;
    double step_B = 1;
    double step_C = 1;
    MatrixXd B_old = B;
    VectorXd C_old = C;

    bool bad_step = true;

    int counter = 0;
    while (bad_step && (counter < max_steps))
    {
        counter++;

        B = B_old - step * step_B * grad_B;
        if (lambda > 0)
            B = SoftThreshold(B, step * step_B * lambda, n, active_vars);
        if (lambda2 > 0)
            tie(pairwise_error, B) = PairwiseLassoErrorSoftThreshold(B, step * step_B * lambda2, n, active_vars);

        bad_step = check_instablity(B, n);
        step_B *= alpha;
    }

    C = 0 * C;
    while ((C.array() <= 0).any())
    {
        C = C_old - step * step_C * grad_C;
        step_C *= alpha;
    }
    bad_step = true;
    double f_old = f;

    while (bad_step)
    {
        counter++;

        if (counter > max_steps)
        {
            cout<<f_old<<" "<<f<<endl;
            cout<<B_old<<endl;
            cout<<B<<endl;
            cout<<grad_B<<endl;
            return make_tuple(B, C, 1);
        }

        bad_step = false;

        if (fitch_or_verando)
            tie(f, temp1, temp2) = ForbeniusLossGradient(Sigma, B, C, n, lambda, lambda2, kappa, active_vars);
        else
            tie(f, temp1, temp2) = LassoGradient(Sigma, B, C, n, lambda, lambda2, kappa, active_vars);

        if ((f_old - f) < 0)
        {
            step *= alpha;
            C = C_old - step * step_C * grad_C;
            B = B_old - step * step_B * grad_B;
            if (lambda > 0)
                B = SoftThreshold(B, step * step_B * lambda, n, active_vars);
            if (lambda2 > 0)
                tie(pairwise_error, B) = PairwiseLassoErrorSoftThreshold(B, step * step_B * lambda2, n, active_vars);
            bad_step = true;
        }
        if ((f_old - f) < min_change)
            return make_tuple(B, C, 2);
    }

    cout<<"f old:"<< f_old<<" l1: "<<lambda*B.norm()<<" Lambda: "<<lambda<<"   B norm: "<<B.norm()<<"  C norm"<<C.norm()<<endl;
    cout<<"B old norm: "<<B_old.norm()<<endl;
    cout<<"f new:"<<f<<" number of steps "<<counter<<" step_b counter:"<<counter<<endl;
    cout<<"#############################################################################"<<endl;
    return make_tuple(B, C, 0);
}

tuple<MatrixXd, VectorXd, int> FitLyapunov(MatrixXd Sigma, MatrixXd B, VectorXd C, int n, double alpha, double lambda, double lambda2, double kappa, MatrixXd active_vars, double min_change, int max_steps, bool fitch_or_verando, int max_decent)
{
    int converge_flag = 0;
    int count = 0;
    MatrixXd B_new;
    VectorXd C_new;

    while ((converge_flag == 0) and (count <= max_decent))
    {
        count++;

        tie(B_new, C_new, converge_flag) = ComputeGradient(Sigma, B, C, n, alpha, lambda, lambda2, kappa, active_vars, min_change, max_steps, fitch_or_verando);
        B = B_new;
        C = C_new;
    }
    return make_tuple(B_new, C_new, converge_flag);
}

MatrixXd Generate_A_Matrix(MatrixXd active_vars, MatrixXd Sigma, int n,bool diagonal_c)
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

    }
    return output;
}