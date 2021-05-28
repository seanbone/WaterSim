#ifndef __conjugate_gradient_hpp_
#define __conjugate_gradient_hpp_

// NOTE included new header
#include<algorithm>
#include<cmath>
#include "Mac3d.h"
// incomplete cholesky conjugate gradient
struct SparseMat {
	SparseMat(unsigned a, unsigned b);
	~SparseMat();

	// copy constructor
	// SparseMat(const SparseMat A);

	// values, rows
	const unsigned v,r;
	double *values;
	unsigned *col_idx;
	unsigned *row_idx;
};

class ICConjugateGradientSolver {
	const Mac3d& grid;
	const unsigned n_cells_x, n_cells_y, n_cells_z;

	const unsigned stride_x = 1;
	const unsigned stride_y = n_cells_x;
	const unsigned stride_z = n_cells_x * n_cells_y;

	const double tau;
	public:
	// number of rows aka. len of rhs aka. len of res, guess vector ect.
	const unsigned num_cells;

	protected:
	// intermediate vector for pressure solution
	double *q;

	// residual vector
	double *r;

	// auxiliary vector
	double *z;

	// search vector
	double *s;

	// diagonal of preconditioner matrix
	double *precon_diag;

	// diagonal of matrix A
	double* A_diag;

	// current step and max steps
	unsigned step, max_steps;

	// threshhold
	const double thresh = 1e-9;

	public:
	ICConjugateGradientSolver();
	~ICConjugateGradientSolver();
	ICConjugateGradientSolver(unsigned max_steps, const Mac3d& grid);

	void computePreconDiag();
	void applyPreconditioner(const double *r, double *z) const;
	void applyA(const double *s, double *z) const;

    void solve(const double* rhs, double* p);
};




#endif
