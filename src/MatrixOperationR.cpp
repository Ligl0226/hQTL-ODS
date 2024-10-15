// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace std;
using namespace Eigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Transpose;
using Rcpp::wrap;
using Rcpp::Named;
using Eigen::Lower;
using Eigen::LLT;
using Eigen::HouseholderQR;
using Eigen::SVDBase;

// transpose of matrix A
// [[Rcpp::export]]
SEXP trans_tX_Cpp(const Eigen::Map<Eigen::MatrixXd> A){
    Eigen::MatrixXd tA = A.transpose();
    return Rcpp::wrap(tA);
}

// product of two maritx A and B
// [[Rcpp::export]]
SEXP prod_XY_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
	Eigen::MatrixXd C = A * B;
	return Rcpp::wrap(C);
}

// product of 3 maritx A, B and C
// [[Rcpp::export]]
SEXP prod_XYZ_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B, const Eigen::Map<Eigen::MatrixXd> C){
	Eigen::MatrixXd M = A * B * C;
	return Rcpp::wrap(M);
}

// product of 5 maritx A, B, C, D and E
// [[Rcpp::export]]
SEXP prod_ABCDE_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B, const Eigen::Map<Eigen::MatrixXd> C, const Eigen::Map<Eigen::MatrixXd> D, const Eigen::Map<Eigen::MatrixXd> E){
	Eigen::MatrixXd M = A * B * C * D * E;
	return Rcpp::wrap(M);
}

// product of 6 maritx A, B, C, D, E and F
// [[Rcpp::export]]
SEXP prod_ABCDEF_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B, const Eigen::Map<Eigen::MatrixXd> C, const Eigen::Map<Eigen::MatrixXd> D, const Eigen::Map<Eigen::MatrixXd> E, const Eigen::Map<Eigen::MatrixXd> F){
	Eigen::MatrixXd M = A * B * C * D * E * F;
	return Rcpp::wrap(M);
}

// the variance of vector x
// [[Rcpp::export]]
SEXP vecvar_x_Cpp(const Eigen::Map<Eigen::VectorXd> x){
	VectorXd centered = x.array() - x.mean();
    double result = (centered.array() * centered.array()).sum() / (x.size() - 1);
	return Rcpp::wrap(result);
}

// the variance of each column of matrix X
// [[Rcpp::export]]
SEXP colvar_X_Cpp(const Eigen::Map<Eigen::MatrixXd> X){
	MatrixXd centered = X.rowwise() - X.colwise().mean();
    VectorXd colvar = (centered.array() * centered.array()).colwise().sum() / (X.rows() - 1);
	return Rcpp::wrap(colvar);
}

// pairwise cor between matrix X and Y (based on the Rcode from YongJiang)
// [[Rcpp::export]]
SEXP pairwisecor_XY_Cpp(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Y){
	MatrixXd X1 = X.rowwise() - X.colwise().mean();
    MatrixXd Y1 = Y.rowwise() - Y.colwise().mean();
	VectorXd sdx1 = (X1.array() * X1.array()).colwise().mean().sqrt();
    VectorXd sdy1 = (Y1.array() * Y1.array()).colwise().mean().sqrt();
	VectorXd z_numerator = (X1.array() * Y1.array()).colwise().mean();
    VectorXd z_denominator = sdx1.array() * sdy1.array();
	VectorXd z = z_numerator.array() / z_denominator.array();
	return Rcpp::wrap(z);
}

// the wise-element product of matrix X and Y
// [[Rcpp::export]]
SEXP wiseprod_XY_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
	Eigen::MatrixXd C = A.array() * B.array();
	return Rcpp::wrap(C);
}

// the wise-element product of matrix X and X
// [[Rcpp::export]]
SEXP wiseprod_XX_Cpp(const Eigen::Map<Eigen::MatrixXd> A){
	Eigen::MatrixXd C = A.array() * A.array();
	return Rcpp::wrap(C);
}

// the wise-element product of matrix X and 1/Y
// [[Rcpp::export]]
SEXP wiseprod_XYinv_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
	Eigen::MatrixXd C = A.array() / B.array();
	return Rcpp::wrap(C);
}

// colsum for the wise-element product of X and Y
// [[Rcpp::export]]
SEXP wiseprod_colsum_XY_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
	Eigen::VectorXd C = (A.array() * B.array()).colwise().sum().array();
	return Rcpp::wrap(C);
}

// colsum for the wise-element product of X and X
// [[Rcpp::export]]
SEXP wiseprod_colsum_X2_Cpp(const Eigen::Map<Eigen::MatrixXd> A){
	Eigen::VectorXd C = (A.array() * A.array()).colwise().sum().array();
	return Rcpp::wrap(C);
}

// triple rank among Matrix X,Y and Z (based on the Rcode from YongJiang)
// [[Rcpp::export]]
SEXP mattriplerank_XYZ_Cpp(const Eigen::Map<Eigen::MatrixXd> X,
						   const Eigen::Map<Eigen::MatrixXd> Y,
						   const Eigen::Map<Eigen::MatrixXd> Z){
	ArrayXd z(X.cols());
    for (int i = 0; i < X.cols(); i++){
        MatrixXd tempMet(X.rows(),3);
        tempMet << X.col(i), Y.col(i), Z.col(i);
        JacobiSVD<MatrixXd> svd(tempMet.jacobiSvd(ComputeThinU | ComputeThinV));
        z(i) = (svd.singularValues().array().abs().array() > 1e-8).count();
    }
	return Rcpp::wrap(z);
}

// the product of matrix t(X) and Y
// [[Rcpp::export]]
SEXP crossprod_tXY_Cpp(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
	Eigen::MatrixXd C = A.adjoint() * B;
	return Rcpp::wrap(C);
}

// the produnt of matrix t(X) and X
// [[Rcpp::export]]
SEXP crossprod_tXX_Cpp(const Eigen::Map<Eigen::MatrixXd> A){
	int ncol = A.cols();
	Eigen::MatrixXd tAA(Eigen::MatrixXd(ncol, ncol).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()));
	return Rcpp::wrap(tAA);
}

// the product of matrix X and t(X)
// [[Rcpp::export]]
SEXP crossprod_XtX_Cpp(const Eigen::Map<Eigen::MatrixXd> A){
	int nrow = A.rows();
    Eigen::MatrixXd AtA(Eigen::MatrixXd(nrow, nrow).setZero().selfadjointView<Lower>().rankUpdate(A));
	return Rcpp::wrap(AtA);
}

// the product of matrix X and t(X) with n_cores
// [[Rcpp::export]]
SEXP crossprod_XtX_ncores_Cpp(const Eigen::Map<Eigen::MatrixXd> A, int n_cores){
	Eigen::setNbThreads(n_cores);
	int nrow = A.rows();
    Eigen::MatrixXd AtA(Eigen::MatrixXd(nrow, nrow).setZero().selfadjointView<Lower>().rankUpdate(A));
	return Rcpp::wrap(AtA);
}

// wise product of vector x and matrix Y
// [[Rcpp::export]]
SEXP wiseprod_VxMY_Cpp(const Eigen::Map<Eigen::VectorXd> x, const Eigen::Map<Eigen::MatrixXd> Y){
	Eigen::MatrixXd Z;
	Z = Y.array().colwise() * x.array();
	return Rcpp::wrap(Z);
}