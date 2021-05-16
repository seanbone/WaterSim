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

void ICConjugateGradientSolver::computePreconDiag() {
    unsigned cellidx = 0;
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++, cellidx++) {
				double e = A_diag[cellidx]
					- std::pow(-1 * precon_diag[cellidx-1], 2)
					- std::pow(-1 * precon_diag[cellidx-n_cells_x], 2)
					- std::pow(-1 * precon_diag[cellidx-n_cells_x * n_cells_y], 2);
				precon_diag[cellidx] = 1 / std::sqrt(e + 1e-30);
			}
		}
	}
}

// apply the matrix A
// element order: k-j-i
void ICConjugateGradientSolver::applyA(const double *z, double *s) {
    unsigned cellidx = 0;
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++, cellidx++) {
                // Copy diagonal entry
                auto& diag_e = grid.get_a_diag()[i + j*n_cells_x + n_cells_x*n_cells_y*k];
                // triplets.push_back(diag_e);
				s[diag_e.row()] += diag_e.value() * z[diag_e.col()];

                // Compute off-diagonal entries
                if (grid.is_fluid(i, j, k)){
                    // x-adjacent cells
                    if (i+1 < n_cells_x && grid.is_fluid(i+1, j, k)){

                        // Compute (i+1,j,k)
                        // triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx+1, -1));
						s[cellidx] += -1 * z[cellidx+1];

                        // Use symmetry to avoid computing (i-1,j,k) separately
                        // triplets.push_back(Mac3d::Triplet_t(cellidx+1, cellidx, -1));
						s[cellidx+1] += -1 * z[cellidx];
                    }

                    // y-adjacent cells
                    if (j+1 < n_cells_y && grid.is_fluid(i, j+1, k)){

                        // Compute (i,j+1,k)
                        // triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + n_cells_x, -1));
						s[cellidx] += -1 * z[cellidx+n_cells_x];

                        // Use symmetry to avoid computing (i,j-1,k) separately
                        // triplets.push_back(Mac3d::Triplet_t(cellidx + n_cells_x, cellidx, -1));
						s[cellidx+n_cells_x] += -1 * z[cellidx];
                    }

                    // z-adjacent cells
                    if (k+1 < n_cells_z && grid.is_fluid(i, j, k+1)){

                        // Compute (i,j,k+1)
                        // triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + n_cells_x*n_cells_y, -1));
						s[cellidx] += -1 * z[cellidx+n_cells_x*n_cells_y];

                        // Use symmetry to avoid computing (i,j-1) separately
                        // triplets.push_back(Mac3d::Triplet_t(cellidx + n_cells_x*n_cells_y, cellidx, -1));
						s[cellidx+n_cells_x*n_cells_y] += -1 * z[cellidx];
                    }
				}
			}
		}
	}
}

void cg::ICConjugateGradientSolver::solve(const double* rhs, double* p) {
    // initialize initial guess and residual
    std::fill(p,p+num_cells,0);
    std::copy(rhs,rhs+num_cells,r);

	applyPreconditioner(r, z);
    std::copy(z,z+num_cells,s);

    double sigma = dot_product(z, r, num_cells);

	// load A_diag

    for(unsigned step = 0; step < max_steps; step++){
		applyA(z, s);

        double dots = dot_product(z,s,num_cells);
        double alpha = rho/dots;

        sca_add_product(s,alpha,num_cells,p);
        sca_add_product(z,-alpha,num_cells,r);

        //check if exit cond is met : inf norm < then some thresh
		{
			bool thresh_exceeded = false;
			for (double* el = r; el < r+num_cells; el++) {
				if (std::abs(*el) > thresh) thresh_exceeded = true;
			}
			if (not thresh_exceeded) return;
		}

		applyPreconditioner(r, z);
		sigma = dot_product(z, r, num_cells);
        double beta = sigma * rho_inv;
        //Bug potential: aliasing
        sca_product(s, beta, z, num_cells, s);
    }
}
