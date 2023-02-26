#ifndef PLAYER_H 
#define PLAYER_H

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
using Eigen::ArrayXXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;

Matrix2d Solve2dSyl(Matrix2d R1,Matrix2d R2,Matrix2d X);
MatrixXd SolveSides(MatrixXd D,MatrixXd U, MatrixXd R_22, vector<int> starts,vector<int> sizes,int last);
MatrixXd SolveLypunov (MatrixXd B,MatrixXd C,int n);
MatrixXd SoftThreshold(MatrixXd B, double lambda, int n, MatrixXd active_vars);
tuple<double,MatrixXd> PairwiseLassoErrorSoftThreshold(MatrixXd B, double lambda, int n, MatrixXd active_vars);
MatrixXd SoftThreshold(MatrixXd B, double lambda, int n);
tuple<MatrixXd,VectorXd,int> ComputeGradient(MatrixXd Sigma, MatrixXd B, VectorXd C,int n, double alpha,double lambda,double lambda2,double kappa, MatrixXd active_vars ,double min_change=1e-5, int max_steps=100,bool fitch_or_varando=true);
tuple<MatrixXd,VectorXd,int> FitLyapunov(MatrixXd Sigma, MatrixXd B, VectorXd C,int n, double alpha,double lambda,double lambda2, double kappa,MatrixXd active_vars,double min_change, int max_steps, bool fitch_or_verand,int max_decent);
MatrixXd Generate_A_Matrix(MatrixXd active_vars, MatrixXd Sigma, int n,bool diagonal_c=false);
tuple<int,int,int,int>  accuracy(MatrixXd A,MatrixXd B,int n);
tuple<int,int,int,int>  pairwise_accuracy(MatrixXd A,MatrixXd B,int n);



#endif