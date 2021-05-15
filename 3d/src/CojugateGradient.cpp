#include "ConjugateGradient.hpp"

using namespace cg;

// apply the preconditioner (L L^T)^-1 by solving Lq = d and Lp = q
void ICConjugateGradientSolver::applyPreconditioner(const double *r, double *z) {
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++) {
			}
		}
	}
}

// apply the matrix A
// element order: k-j-i
void ICConjugateGradientSolver::applyA(const double *z, double *s) {
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++) {
			}
		}
	}
}

void cg::ICConjugateGradientSolver::solve(const double* rhs, double* p) {
    // initialize initial guess and residual
    std::fill(p,p+num_rows,0);
    std::copy(rhs,rhs+num_rows,r);

	applyPreconditioner(r, z);
    std::copy(z,z+num_rows,s);

    double sigma = dot_product(z, r, num_rows);


    for(unsigned step = 0; step < max_steps; step++){
		applyA(z, s);

        double dots = dot_product(z,s,num_rows);
        double alpha = rho/dots;

        sca_add_product(s,alpha,num_rows,p);
        sca_add_product(z,-alpha,num_rows,r);

        //check if exit cond is met : inf norm < then some thresh
		{
			bool thresh_exceeded = false;
			for (double* el = r; el < r+num_rows; el++) {
				if (std::abs(*el) > thresh) thresh_exceeded = true;
			}
			if (not thresh_exceeded) return;
		}

		applyPreconditioner(r, z);
		sigma = dot_product(z, r, num_rows);
        double beta = sigma * rho_inv;
        //Bug potential: aliasing
        sca_product(s, beta, z, num_rows, s);
    }
}
