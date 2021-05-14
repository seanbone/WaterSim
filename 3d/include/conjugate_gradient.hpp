namespace cg {
	// incomplete cholesky conjugate gradient
	class ICConjugateGradientSolver {
public:
		void solve(
				const double* A_vals, const int* A_row_idx, const int* A_col_idx,
				const double* rhs, double* sol) {
			// for the matrix A, defined by (A_vals, A_row_idx, A_col_idx)
			// and the right hand side rhs, compute the solution using
			// a conjugate gradient solver with incomplete cholesky preconditioning.
			
		}
	};
}
