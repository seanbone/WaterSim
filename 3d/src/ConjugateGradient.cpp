#include "ConjugateGradient.hpp"
#include <cassert>
#include <cmath>
#include <iostream>


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
                auto& diag_e = grid.get_a_diag()[i + j*stride_y + k*stride_z];
				A_diag[cellidx] = diag_e.value();
			}
		}
	}

}

// apply the preconditioner (L L^T)^-1 by solving Lq = d and Lp = q
void ICConjugateGradientSolver::applyPreconditioner(const double *r, double *z) const {
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++) {
				const unsigned cellidx = i + j*stride_y + k*stride_z;
				if (not grid.is_fluid(i, j, k)) {
					q[cellidx] = 0;
					continue;
				}
				double t = r[cellidx];
				if (i > 0 && grid.is_fluid(i-1, j, k)) t += precon_diag[cellidx - stride_x] * q[cellidx - stride_x];
				if (j > 0 && grid.is_fluid(i, j-1, k)) t += precon_diag[cellidx - stride_y] * q[cellidx - stride_y];
				if (k > 0 && grid.is_fluid(i, j, k-1)) t += precon_diag[cellidx - stride_z] * q[cellidx - stride_z];
				q[cellidx] = t * precon_diag[cellidx];
			}
		}
	}
	for (int k = n_cells_z-1; k >= 0; k--) {
		for (int j = n_cells_y-1; j >= 0; j--) {
			for (int i = n_cells_x-1; i >= 0; i--) {
				const unsigned cellidx = i + j*stride_y + k*stride_z;
				if (not grid.is_fluid(i, j, k)) {
					z[cellidx] = 0;
					continue;
				}
				double t = q[cellidx];
				if (i + 1 < (int) n_cells_x && grid.is_fluid(i+1, j, k)) t += precon_diag[cellidx] * z[cellidx + stride_x];
				if (j + 1 < (int) n_cells_y && grid.is_fluid(i, j+1, k)) t += precon_diag[cellidx] * z[cellidx + stride_y];
				if (k + 1 < (int) n_cells_z && grid.is_fluid(i, j, k+1)) t += precon_diag[cellidx] * z[cellidx + stride_z];
				z[cellidx] = t * precon_diag[cellidx];
			}
		}
	}
}

void ICConjugateGradientSolver::computePreconDiag() {
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++) {
				const unsigned cellidx = i + j*stride_y + k*stride_z;
				if (not grid.is_fluid(i, j, k)) {
					precon_diag[cellidx] = 0;
					continue;
				}
				double e = A_diag[cellidx];
				if (i > 0 && grid.is_fluid(i-1, j, k)) e -= std::pow(-1 * precon_diag[cellidx-stride_x], 2);
				if (j > 0 && grid.is_fluid(i, j-1, k)) e -= std::pow(-1 * precon_diag[cellidx-stride_y], 2);
				if (k > 0 && grid.is_fluid(i, j, k-1)) e -= std::pow(-1 * precon_diag[cellidx-stride_z], 2);
				precon_diag[cellidx] = 1 / std::sqrt(e + 1e-30);
			}
		}
	}
}

// apply the matrix A: y <- A b
// element order: k-j-i
void ICConjugateGradientSolver::applyA(const double *b, double *y) const{
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++) {
				const unsigned cellidx = i + j*stride_y + k*stride_z;
				y[cellidx] = 0;
			}
		}
	}
	for (unsigned k = 0; k < n_cells_z; k++) {
		for (unsigned j = 0; j < n_cells_y; j++) {
			for (unsigned i = 0; i < n_cells_x; i++) {
				const unsigned cellidx = i + j*stride_y + k*stride_z;
                auto& diag_e = grid.get_a_diag()[cellidx];
				y[diag_e.row()] += diag_e.value() * b[diag_e.col()];

                // Compute off-diagonal entries
                if (grid.is_fluid(i, j, k)){
                    // x-adjacent cells
                    if (i+1 < n_cells_x && grid.is_fluid(i+1, j, k)){
						y[cellidx] 			+= -1 * b[cellidx+stride_x];
						y[cellidx+stride_x] += -1 * b[cellidx];
                    }

                    // y-adjacent cells
                    if (j+1 < n_cells_y && grid.is_fluid(i, j+1, k)){
						y[cellidx] 			+= -1 * b[cellidx+stride_y];
						y[cellidx+stride_y] += -1 * b[cellidx];
                    }

                    // z-adjacent cells
                    if (k+1 < n_cells_z && grid.is_fluid(i, j, k+1)){
						y[cellidx] 			+= -1 * b[cellidx+stride_z];
						y[cellidx+stride_z] += -1 * b[cellidx];
                    }
				}
			}
		}
	}
}

void checknan(const double* array, int len, std::string array_name="array") {
	int count_left = 20;
	for (int i = 0; i<len; i++) {
		if(std::isnan(array[i])) {
			std::cout
				<< "WARNING: NaN found in " << array_name <<
				" at position " << i << "!" << std::endl;
			if (--count_left == 0)
				std::cout << "[further occurrences ignored]" << std::endl;
		}
	}
}

void print_array_head(const double* array, std::string prefix, unsigned number) {
	std::cout << prefix;
	for (unsigned i=0; i < number; i++) std::cout << array[i] << ' ';
	std::cout << std::endl;
}

void cg::ICConjugateGradientSolver::solve(const double* rhs, double* p) {
    // initialize initial guess and residual
	print_array_head(rhs,  "RHS: ");

	// p = 0
    std::fill(p,p+num_cells,0);
	// r = d
    std::copy(rhs,rhs+num_cells,r);
	// catch zero rhs early
	{
		double max_abs_val = 0;
		for (double* el = r; el < r+num_cells; el++) {
			const double abs_val = std::abs(*el);
			if (max_abs_val < abs_val) max_abs_val = abs_val;
		}
		// std::cout << "Max abs val: " << max_abs_val << std::endl;
		if (max_abs_val < thresh) {
			std::cout << "Number of steps: " << step + 1 << std::endl;
			return;
		}
	}

	computePreconDiag();
	// s = M⁻¹ r
	applyPreconditioner(r, s);

	// precompute the diagonal of the preconditioner
	print_array_head(precon_diag,  "PreconDiag: ");

	// ρ = <r,s>
	double rho = dot_product(r,s,num_cells);

    for(unsigned step = 0; step < max_steps; step++){
		// std::cout << "\nNew Round! (" << step << ")\n";
		// print_array_head(p,  "p: ");
		// print_array_head(r,  "r: ");
		// z = A s
		applyA(s, z);
		//checknan(s, num_cells, "s");
		// print_array_head(z,  "z: ");

		// α = ρ / <s,z>
        const double dots = dot_product(z,s,num_cells);
		const double alpha = rho / dots;

        sca_add_product(s,alpha,num_cells,p);
        sca_add_product(z,-alpha,num_cells,r);

        //check if exit cond is met : inf norm < then some thresh
		{
			double max_abs_val = 0;
			for (double* el = r; el < r+num_cells; el++) {
				const double abs_val = std::abs(*el);
				if (max_abs_val < abs_val) max_abs_val = abs_val;
			}
			// std::cout << "Max abs val: " << max_abs_val << std::endl;
			if (max_abs_val < thresh) {
				std::cout << "Number of steps: " << step + 1 << std::endl;
				return;
			}
			else if (step + 1 == max_steps) {
				std::cout << "WARNING: Failed to find a solution in " << step + 1 << " steps (|r|=" << max_abs_val <<")!" << std::endl;
			}
		}

		// z = M⁻¹ r
		applyPreconditioner(r, z);
		const double rho_new = dot_product(z, r, num_cells);
		const double beta = rho_new / rho;
		rho = rho_new;
        //Bug potential: aliasing
        sca_product(s, beta, z, num_cells, s);
    }
}

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
    //delete [] r;
    delete [] z;
    delete [] s;
    delete [] precon_diag;
    delete [] A_diag;

}
