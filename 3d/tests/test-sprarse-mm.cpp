/*
 * A test to check the sparse matrix multiplication for the MICCG solver
 */
#include "conjugate_gradient.hpp"
#include <Eigen/Sparse>	//used for the matrix A

cg::SparseMat from_eigen(const Eigen::SparseMatrix<double, Eigen::RowMajor> A_eigen) {
	// construct a sparse matrix from an eigen sparse matrix
	const unsigned int n_values = A_eigen.nonZeros();
	const unsigned int n_row_idx = A_eigen.rows()+1;
	cg::SparseMat A(n_values, n_row_idx);

	std::copy(A_eigen.valuePtr(), A_eigen.valuePtr()+n_values, A.values);
	std::copy(A_eigen.innerIndexPtr(), A_eigen.innerIndexPtr()+n_values, A.col_idx);
	std::copy(A_eigen.outerIndexPtr(), A_eigen.outerIndexPtr()+n_row_idx, A.row_idx);
	return A;
}


Eigen::SparseMatrix<double, Eigen::RowMajor> to_eigen(const cg::SparseMat A, int n_cols) {
	// construct a sparse matrix from an eigen sparse matrix
	const unsigned int n_values = A.v;
	const unsigned int n_row_idx = A.r;

	Eigen::SparseMatrix<double, Eigen::RowMajor> A_eigen(n_row_idx-1, n_cols);

	std::copy(A.values, A.values+n_values, A_eigen.valuePtr());
	std::copy(A.col_idx, A.col_idx+n_values,  A_eigen.innerIndexPtr());
	std::copy(A.col_idx, A.col_idx+n_row_idx, A_eigen.outerIndexPtr());
	return A_eigen;
}


int main() {
	Eigen::SparseMatrix<double, Eigen::RowMajor> A_eigen(4,4);
	A_eigen.setIdentity();
	cg::SparseMat A = from_eigen(A_eigen);

	Eigen::VectorXd b_eigen;
	b_eigen << 1, 2, 3, 4;

	const double b[] = {1, 2, 3, 4};

	Eigen::VectorXd result_eigen = Eigen::VectorXd(A_eigen.dot(b_eigen));

	double result_test[4];
	cg::ICConjugateGradientSolver::Mat_mult(&A, b, result_test)
}
