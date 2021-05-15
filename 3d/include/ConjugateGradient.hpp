#ifndef __conjugate_gradient_hpp_
#define __conjugate_gradient_hpp_

// NOTE included new header
#include<algorithm>

namespace cg {
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
	// incomplete cholesky conjugate gradient
	class ICConjugateGradientSolver {

	    // preconditioner Matrix
	    SparseMat M;

	    // number of rows aka. len of rhs aka. len of res, guess vector ect.
        unsigned num_rows;

        // initial guess
	    double *p;

	    // residual vector
	    double *r;

	    // auxiliary vector
	    double *z;

	    // search vector
	    double *s;

	    // sigmas, alpha and beta
	    double sig, sig_new, alpha, beta;

        // current step and max steps
	    unsigned step, max_steps;

	    // threshhold
	    const double thresh = 1.e-9;

	    // unsure where rho comes from
	    const double rho = 1.;
	    const double rho_inv = 1.;

	    // result of dotprod of s and z
	    double dots;




        public:
        ICConjugateGradientSolver();
        ~ICConjugateGradientSolver();
        ICConjugateGradientSolver(SparseMat Precon, unsigned max_steps);

        // not sure where to put this
        //res[i] = a[i] * b + c[i]
        inline void sca_product(double* a, double b, double* c,  const unsigned n, double *res);

        // res[] += a[] * b
        inline void sca_add_product(double * a, double b, const unsigned n, double *res);

        inline void dot_product(double * a, double * b, const unsigned n, double *res);

        static void Mat_mult(SparseMat *M, const double *a, double *res);

        void solve( cg::SparseMat* A, const double* rhs, double* sol);
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
cg::ICConjugateGradientSolver::ICConjugateGradientSolver(cg::SparseMat Precon, unsigned max_steps)
        :M(Precon), max_steps(max_steps) {
    // the -1 is beacuse we assume a past the end entry for r
    num_rows = M.r - 1;
    step = 0;
    p = new double [num_rows];
    r = new double [num_rows];
    z = new double [num_rows];
    s = new double [num_rows];


}
cg::ICConjugateGradientSolver::~ICConjugateGradientSolver() {
    delete [] p;
    delete [] r;
    delete [] z;
    delete [] s;


}

void cg::ICConjugateGradientSolver::solve( cg::SparseMat* A, const double* rhs, double* sol) {
    // initial phase
    std::copy(rhs,rhs+num_rows,r);
    // unsure about pass by ref
    Mat_mult(&M,r,z);
    std::copy(z,z+num_rows,s);
    dot_product(z,r,num_rows, &sig);


    while( step < max_steps){
        Mat_mult(A,s,z);
        dot_product(z,s,num_rows,&dots);
        alpha = rho/dots;
        sca_add_product(s,alpha,num_rows,p);
        // r [] -= alph * z[]
        sca_add_product(z,-alpha,num_rows,r);

        //check if exit cond is met : inf norm < then some thresh
        if( *std::max_element(r,r+num_rows) < thresh)
            sol = p; return;
        Mat_mult(&M,r,z);
        dot_product(z,r,num_rows,&sig_new);

        beta = sig_new * rho_inv;
        //Bug potential: aliasing
        sca_product(s,beta,z,num_rows,s);
    }
    sol = p;
}

inline void cg::ICConjugateGradientSolver::dot_product(double *a, double *b, const unsigned int n, double *res) {
    double tmp = 0;
    for(unsigned i = 0 ; i < n; ++i ) {
        tmp += a[i]*b[i];
    }
    *res= tmp;
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

inline void cg::ICConjugateGradientSolver::sca_add_product(double *a, double b, const unsigned int n, double *res) {
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
