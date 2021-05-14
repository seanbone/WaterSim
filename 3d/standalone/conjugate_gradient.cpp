#include "conjugate_gradient.hpp"

// standalone file for development of conjugate gradient solver


int main() {
	// build a test matrix A
	constexpr int nnz = 20;
	constexpr int n_rows = 7;
	double values[nnz];
	int col_idx[nnz];
	int row_idx[n_rows];

	// build a rhs
	double rhs[n_rows];

	// instaciate the solver
	cg::ICConjugateGradientSolver solver;

	// solve the sparse LSE
	double solution[n_rows];

	solver.solve(values, row_idx, col_idx, rhs, solution);

}
