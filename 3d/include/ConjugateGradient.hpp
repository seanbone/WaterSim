#ifndef __conjugate_gradient_hpp_
#define __conjugate_gradient_hpp_

// NOTE included new header
#include<algorithm>
#include<cmath>
#include "Mac3d.h"

namespace cg {
	// incomplete cholesky conjugate gradient
    struct SparseMat {
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

		const unsigned stride_x = 1;
		const unsigned stride_y = n_cells_x;
		const unsigned stride_z = n_cells_x * n_cells_y;

		const double tau;

		public:
	    // number of rows aka. len of rhs aka. len of res, guess vector ect.
        const unsigned num_cells;

        // pressure
	    double *p;
		
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

	    // unsure where rho comes from

        public:
        ICConjugateGradientSolver();
        ~ICConjugateGradientSolver();
        ICConjugateGradientSolver(unsigned max_steps, const Mac3d& grid);

		void computePreconDiag();
		void applyPreconditioner(const double *r, double *z) const;
		void applyA(const double *s, double *z) const;

        // not sure where to put this
        //res[i] = a[i] * b + c[i]
        inline void sca_product(const double* a, double b, double* c,  const unsigned n, double *res);

        // res[] += a[] * b
        inline void sca_add_product(const double * a, const double b, const unsigned n, double *res);

        inline double dot_product(const double * a, const double * b, const unsigned n);

        static void Mat_mult(SparseMat *M, const double *a, double *res);

        void solve(const double* rhs, double* p);
	};

}
/*
cg::SparseMat::SparseMat(SparseMat &M): v(M.v), r(M.r){
    double values[v];
    int col_idx[v];
    int row_idx[r];
};
*/

inline double cg::ICConjugateGradientSolver::dot_product(const double *a, const double *b, const unsigned int n) {
    double tmp = 0;
    for(unsigned i = 0 ; i < n; ++i ) {
        tmp += a[i]*b[i];
    }
    return tmp;
}

inline void cg::ICConjugateGradientSolver::Mat_mult(SparseMat *M, const double *a, double *res) {
    const unsigned v = M->v;
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
        ++curr_row_idx;
    }
}

inline void cg::ICConjugateGradientSolver::sca_add_product(const double *a, const double b, const unsigned int n, double *res) {
    for(unsigned i = 0 ; i < n; ++i ) {
        res [i] += a[i] * b;
    }
}
inline void cg::ICConjugateGradientSolver::sca_product(const double *a, double b,double *c, const unsigned int n, double *res) {
    for(unsigned i = 0 ; i < n; ++i ) {
        res[i] = a[i] * b + c[i];
    }
}

void print_array_head(const double* array, std::string prefix="", unsigned number=20);
#endif
