#include "ConjugateGradient.hpp"
#include <cassert>

using namespace cg;

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
    q = new double [num_cells];
    r = new double [num_cells];
    z = new double [num_cells];
    s = new double [num_cells];
    precon_diag = new double [num_cells];
    A_diag = new double [num_cells];

	// load A_diag
    unsigned cellidx = 0;
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++, cellidx++) {
                auto& diag_e = grid.get_a_diag()[i + j*n_cells_x + n_cells_x*n_cells_y*k];
				A_diag[cellidx] = diag_e.value();
			}
		}
	}

}

// apply the preconditioner (L L^T)^-1 by solving Lq = d and Lp = q
void ICConjugateGradientSolver::applyPreconditioner(const double *r, double *z) {
    unsigned cellidx = stride_x + stride_y + stride_z;
	// element 1,1,1
	for (unsigned k = 1; k < n_cells_z; k++, cellidx += stride_z) {
		for (unsigned j = 1; j < n_cells_y; j++, cellidx += stride_y) {
			for (unsigned i = 1; i < n_cells_x; i++, cellidx++) {
				double t = r[cellidx]
					- (-1 * precon_diag[cellidx - stride_x] * q[cellidx - stride_x])
					- (-1 * precon_diag[cellidx - stride_y] * q[cellidx - stride_y])
					- (-1 * precon_diag[cellidx - stride_z] * q[cellidx - stride_z]);
				q[cellidx] = t * precon_diag[cellidx];
			}
		}
	}
	for (unsigned k = n_cells_z-1; k >= 1; k--, cellidx -= stride_z) {
		for (unsigned j = n_cells_y-1; j >= 1; j--, cellidx -= stride_y) {
			for (unsigned i = n_cells_x-1; i >= 1; i--, cellidx--) {
				double t = q[cellidx]
					- (-1 * precon_diag[cellidx] * z[cellidx + stride_x])
					- (-1 * precon_diag[cellidx] * z[cellidx + stride_y])
					- (-1 * precon_diag[cellidx] * z[cellidx + stride_z]);
				z[cellidx] = t * precon_diag[cellidx];
			}
		}
	}
    assert (cellidx == stride_x + stride_y + stride_z);
}

void ICConjugateGradientSolver::computePreconDiag() {
    unsigned cellidx = stride_x + stride_y + stride_z;
	// Q: where is precon_diag[i=0|j=0|k=0] initialized?
	for (unsigned k = 1; k < n_cells_z; k++, cellidx += stride_z) {
		for (unsigned j = 1; j < n_cells_y; j++, cellidx += stride_y) {
			for (unsigned i = 1; i < n_cells_x; i++, cellidx++) {
				const double e = A_diag[cellidx]
					- std::pow(-1 * precon_diag[cellidx-stride_x], 2)
					- std::pow(-1 * precon_diag[cellidx-stride_y], 2)
					- std::pow(-1 * precon_diag[cellidx-stride_z], 2);
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
						s[cellidx] += -1 * z[cellidx+stride_x];

                        // Use symmetry to avoid computing (i-1,j,k) separately
                        // triplets.push_back(Mac3d::Triplet_t(cellidx+1, cellidx, -1));
						s[cellidx+stride_x] += -1 * z[cellidx];
                    }

                    // y-adjacent cells
                    if (j+1 < n_cells_y && grid.is_fluid(i, j+1, k)){

                        // Compute (i,j+1,k)
                        // triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + n_cells_x, -1));
						s[cellidx] += -1 * z[cellidx+stride_y];

                        // Use symmetry to avoid computing (i,j-1,k) separately
                        // triplets.push_back(Mac3d::Triplet_t(cellidx + n_cells_x, cellidx, -1));
						s[cellidx+stride_y] += -1 * z[cellidx];
                    }

                    // z-adjacent cells
                    if (k+1 < n_cells_z && grid.is_fluid(i, j, k+1)){

                        // Compute (i,j,k+1)
                        // triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + n_cells_x*n_cells_y, -1));
						s[cellidx] += -1 * z[cellidx+stride_z];

                        // Use symmetry to avoid computing (i,j-1) separately
                        // triplets.push_back(Mac3d::Triplet_t(cellidx + n_cells_x*n_cells_y, cellidx, -1));
						s[cellidx+stride_z] += -1 * z[cellidx];
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

	// double sigma = dot_product(z, r, num_cells);

	// precompute the diagonal of the preconditioner
	computePreconDiag();

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
		double sigma = dot_product(z, r, num_cells);
        double beta = sigma * rho_inv;
        //Bug potential: aliasing
        sca_product(s, beta, z, num_cells, s);
    }
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


cg::ICConjugateGradientSolver::~ICConjugateGradientSolver() {
    delete [] p;
    delete [] q;
    delete [] r;
    delete [] z;
    delete [] s;
    delete [] precon_diag;
    delete [] A_diag;

}
