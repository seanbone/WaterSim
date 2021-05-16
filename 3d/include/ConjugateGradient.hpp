#ifndef __conjugate_gradient_hpp_
#define __conjugate_gradient_hpp_

// NOTE included new header
#include<algorithm>
#include<cmath>
#include "Mac3d.h"

namespace cg {
	// incomplete cholesky conjugate gradient
    struct SparseMat {
        SparseMat() = default;
        SparseMat(unsigned a, unsigned b);
        ~SparseMat();

        // copy constructor
        // SparseMat(const cg::SparseMat A);

        // values, rows
        const unsigned v,r;
        double *values;
        unsigned *col_idx;
        unsigned *row_idx;
    };
	class ICConjugateGradientSolver {

	    // preconditioner Matrix

		const Mac3d& grid;
		const unsigned n_cells_x, n_cells_y, n_cells_z;

		const double tau;

	    // number of rows aka. len of rhs aka. len of res, guess vector ect.
        const unsigned num_cells;

        // initial guess
	    double *p;

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
	    const double thresh = 1.e-9;

	    // unsure where rho comes from
	    const double rho = 1.;
	    const double rho_inv = 1.;

        public:
        ICConjugateGradientSolver();
        ~ICConjugateGradientSolver();
        ICConjugateGradientSolver(unsigned max_steps, const Mac3d& grid);

		void computePreconDiag();
		void applyPreconditioner(const double *r, double *z);
		void applyA(const double *z, double *s);

        // not sure where to put this
        //res[i] = a[i] * b + c[i]
        inline void sca_product(double* a, double b, double* c,  const unsigned n, double *res);

        // res[] += a[] * b
        inline void sca_add_product(const double * a, const double b, const unsigned n, double *res);

        inline double dot_product(double * a, double * b, const unsigned n);

        static void Mat_mult(SparseMat *M, const double *a, double *res);

        void solve(const double* rhs, double* p);
	};

}
//stuff here can be moved to extra files
cg::SparseMat::SparseMat(unsigned a, unsigned b): v(a), r(b){
    values = new double [v];
    col_idx = new unsigned [v];
    row_idx = new unsigned [r];
}
cg::SparseMat::~SparseMat(){
    delete [] values;
    delete [] col_idx;
    delete [] row_idx;
}
/*
cg::SparseMat::SparseMat(SparseMat &M): v(M.v), r(M.r){
    double values[v];
    int col_idx[v];
    int row_idx[r];
};
*/
cg::ICConjugateGradientSolver::ICConjugateGradientSolver(unsigned max_steps, const Mac3d& grid)
	:
		grid{grid},
		n_cells_x{grid.get_num_cells_x()}, n_cells_y{grid.get_num_cells_y()}, n_cells_z{grid.get_num_cells_z()},
		tau{0.97},
		num_cells{n_cells_x * n_cells_y * n_cells_z},
		max_steps(max_steps)
{
    step = 0;
    p = new double [num_cells];
    r = new double [num_cells];
    z = new double [num_cells];
    s = new double [num_cells];
    precon_diag = new double [num_cells];
}

cg::ICConjugateGradientSolver::~ICConjugateGradientSolver() {
    delete [] p;
    delete [] r;
    delete [] z;
    delete [] s;
    delete [] precon_diag;

}


inline double cg::ICConjugateGradientSolver::dot_product(double *a, double *b, const unsigned int n) {
    double tmp = 0;
    for(unsigned i = 0 ; i < n; ++i ) {
        tmp += a[i]*b[i];
    }
    return tmp;
}

inline void cg::ICConjugateGradientSolver::Mat_mult(SparseMat *M, const double *a, double *res) {
    const unsigned v = M->v;
    const unsigned r = M->r;
    double tmp = 0;
    // we assume first val of row dx is 0
    // and last val is past the end val by convention
    unsigned *curr_row_idx = M->row_idx;
    unsigned *curr_col_idx = M->col_idx;

    for(unsigned j = 0; j < v;){
        tmp = 0;
        while(j < *(curr_row_idx+1)){
            tmp += M->values[j] * a[*curr_col_idx];
            ++j;
            ++curr_col_idx;
        }
        res[*curr_row_idx] = tmp;
    }
}

inline void cg::ICConjugateGradientSolver::sca_add_product(const double *a, const double b, const unsigned int n, double *res) {
    for(unsigned i = 0 ; i < n; ++i ) {
        res [i] += a[i] * b;
    }
}
inline void cg::ICConjugateGradientSolver::sca_product(double *a, double b,double *c, const unsigned int n, double *res) {
    for(unsigned i = 0 ; i < n; ++i ) {
        res[i] = a[i] * b + c[i];
    }
}

#endif
