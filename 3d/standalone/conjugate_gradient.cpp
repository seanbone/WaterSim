#include "conjugate_gradient.hpp"

// standalone file for development of conjugate gradient solver


int main() {
	// build a test matrix A
	constexpr int nnz = 20;
	constexpr int n_rows = 7;
	unsigned max_steps = 100;

	cg::SparseMat A ( 20,7);
    cg::SparseMat Pre ( 20,7);


    // build a rhs
	double rhs[n_rows];

	// instaciate the solver
	cg::ICConjugateGradientSolver solver(Pre, max_steps);

	// solve the sparse LSE
	double solution[n_rows];

	//solver.solve(A.values, A.row_idx, A.col_idx, rhs, solution);
    solver.solve(&A, rhs, solution);


}
