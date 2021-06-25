/*
 * A test to check the correctness of the pressure solver
 */
#include <iostream>
#include <cassert>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers> // solve sparse systems

#include "includes/watersim-test-common.h"
#include "Mac3d.h"
#include "FLIP.h"
#include "NcReader.h"
#include "ConjugateGradient.hpp"
#include "tsc_x86.hpp"


void compute_pressure_matrix(Mac3d* MACGrid_, Eigen::SparseMatrix<double>& A_) {
    // Compute matrix for pressure solve and store in A_
    // See eq. (4.19) and (4.24) in SIGGRAPH notes

    // Vector of triplets to construct pressure matrix with
    std::vector< Mac3d::Triplet_t > triplets;

    // Get total number of cells on each axis
    unsigned nx = MACGrid_->get_num_cells_x();
    unsigned ny = MACGrid_->get_num_cells_y();
    unsigned nz = MACGrid_->get_num_cells_z();

    // Index of the current grid-cell [0, nx*ny*nz[ (in the matrix)
    unsigned cellidx = 0;

    // Iterate over all grid-cells
    for (unsigned k = 0; k < nz; ++k){
        for (unsigned j = 0; j < ny; ++j){
            for (unsigned i = 0; i < nx; ++i, ++cellidx){

                unsigned index = i + j*nx + nx*ny*k;

                // Copy diagonal entry
                triplets.push_back(Mac3d::Triplet_t(index, index, MACGrid_->A_diag_val[index]));

                // Compute off-diagonal entries
                if (MACGrid_->is_fluid(i, j, k)){

                    // x-adjacent cells
                    if (i+1 < nx && MACGrid_->is_fluid(i+1, j, k)){

                        // Compute (i+1,j,k)
                        triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx+1, -1));

                        // Use symmetry to avoid computing (i-1,j,k) separately
                        triplets.push_back(Mac3d::Triplet_t(cellidx+1, cellidx, -1));
                    }

                    // y-adjacent cells
                    if (j+1 < ny && MACGrid_->is_fluid(i, j+1, k)){

                        // Compute (i,j+1,k)
                        triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + nx, -1));

                        // Use symmetry to avoid computing (i,j-1,k) separately
                        triplets.push_back(Mac3d::Triplet_t(cellidx + nx, cellidx, -1));
                    }

                    // z-adjacent cells
                    if (k+1 < nz && MACGrid_->is_fluid(i, j, k+1)){

                        // Compute (i,j,k+1)
                        triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + nx*ny, -1));

                        // Use symmetry to avoid computing (i,j-1) separately
                        triplets.push_back(Mac3d::Triplet_t(cellidx + nx*ny, cellidx, -1));
                    }
                } // if is_fluid(i,j,k)
            }
        } // Outer for
    }
    // Set A_ to zero. A_ could be resized only at the start of the sim
    A_.resize(nx*ny*nz, nx*ny*nz);
    A_.setZero();
    A_.setFromTriplets(triplets.begin(), triplets.end());
}


void compute_pressure_rhs(Mac3d* MACGrid_, Eigen::VectorXd& d_, const double fluid_density_, const double dt) {
    // Compute right-hand side of the pressure equations and store in d_
    // See eq. (4.19) and (4.24) in SIGGRAPH notes
    // Note: u_{solid} = 0

    // Get total number of cells on each axis
    unsigned nx = MACGrid_->get_num_cells_x();
    unsigned ny = MACGrid_->get_num_cells_y();
    unsigned nz = MACGrid_->get_num_cells_z();

    // Alias for MAC Grid
    auto& g = MACGrid_;

    // Set d_ to zero. d_ could be resized only at the start of the sim
    d_.resize(nx*ny*nz);
    d_.setZero();

    // Index of the current grid-cell [0, nx*ny*nz[ (in the matrix)
    unsigned cellidx = 0;

    // Iterate over all grid-cells
    for (unsigned k = 0; k < nz; ++k) {
        for (unsigned j = 0; j < ny; ++j) {
            for (unsigned i = 0; i < nx; ++i, ++cellidx) {
                if (g->is_fluid(i,j,k)) {

                    // Apply the formulas of the SIGGRAPH notes
                    double d_ij = -(g->get_u(i+1,j,k) - g->get_u(i,j,k));
                    d_ij -= g->get_v(i,j+1,k) - g->get_v(i,j,k);
                    d_ij -= g->get_w(i,j,k+1) - g->get_w(i,j,k);

                    // Note: u_{solid} = 0

                    // Check each adjacent cell. If solid, alter term as in (4.24)
                    // Consider cells outside of the boundary as solid

                    // (i+1, j, k)
                    if ((i < (nx-1) && g->is_solid(i+1,j,k)) || i == nx-1) {
                        d_ij += g->get_u(i+1,j,k);
                    }

                    // (i-1, j, k)
                    if ((i > 0 && g->is_solid(i-1,j,k)) || i == 0) {
                        d_ij += g->get_u(i,j,k);
                    }

                    // (i, j+1, k)
                    if ((j < (ny-1) && g->is_solid(i,j+1,k)) || j == ny-1) {
                        d_ij += g->get_v(i,j+1,k);
                    }

                    // (i, j-1, k)
                    if ((j > 0 && g->is_solid(i,j-1,k)) || j == 0) {
                        d_ij += g->get_v(i,j,k);
                    }

                    // (i, j, k+1)
                    if ((k < (nz-1) && g->is_solid(i,j,k+1)) || k == nz-1) {
                        d_ij += g->get_w(i,j,k+1);
                    }

                    // (i, j, k-1)
                    if ((k > 0 && g->is_solid(i,j,k-1)) || k == 0) {
                        d_ij += g->get_w(i,j,k);
                    }

                    d_(cellidx) = fluid_density_ * g->get_cell_sizex() * d_ij / dt;

                } else { // if is_fluid(i,j,k)

                    // Set the entry to zero if the current cell is a
                    // fluid-cell
                    d_(cellidx) = 0;
                }
            }
        }
    }
}


int main(){

	/* We need to perform the following:
	 * instanciate the grid
	 * compute pressure matrix
	 * compute pressure rhs
	 * instanciate solver
	 * run solver
	 */


	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");

	// initialize the FLIP instance which will also give us the MACGrid
	std::cout << "Initializing NcReader" << std::endl;
	NcReader nc_reader("ref-340.nc", "config.json");

	nc_reader.toFlipStructures();

	// initialize the MACGrid
	Mac3d& grid = *(nc_reader.MACGrid);

	// initialize the A_ matrix
	Eigen::SparseMatrix<double> A_;

	// initialize the d_ vector
	Eigen::VectorXd d_;
	
	// initialize the dt value
	double dt = nc_reader.cfg.getTimeStep();
	double fluid_density = nc_reader.cfg.getDensity();
	
	// first step: get the old solver to run on some data
	// ###############  Reference Solver  #################
    // Compute A matrix
	// note that A_ is only written
    compute_pressure_matrix(&grid, A_);
	std::vector<Eigen::Triplet<double>> a_diag = grid.get_a_diag();
    // Compute rhs d
	// note that d_ is only written
    compute_pressure_rhs(&grid, d_, fluid_density, dt);

    // Solve for p: Ap = d (MICCG(0))
    using namespace Eigen;
    using solver_t = ConjugateGradient< SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double> >;


	tsctimer.start_timing("particle_to_grid");
	tsctimer.stop_timing("particle_to_grid", true, "");

	tsctimer.start_timing("eigen solver");
	std::cout << "Solving using the reference solver..." << std::endl;
    solver_t solver;
    solver.setMaxIterations(100);
    solver.compute(A_);
    VectorXd p = solver.solve(d_);
	tsctimer.stop_timing("eigen solver", true, "");

	std::cout << "Vector p after reference solver: ";
	for(int i = 0; i < p.size() % 40; i++) std::cout << p(i) << " ";
	std::cout << std::endl;


	// second step: get the new solver to run on the same data
	// ###############  Optimized Solver  #################
	//

	// get a raw array of the right-hand side
	double* rhs = d_.data();

	ICConjugateGradientSolver cg_solver(100, grid);
	unsigned num_cells = cg_solver.num_cells;
	std::vector<double> p_vec(num_cells);

	// should really be the same size!
	assert (num_cells == d_.size());
	
	tsctimer.start_timing("own pcg solver");
	std::cout << "Solving using the new solver..." << std::endl;
	cg_solver.solve(rhs, p_vec.data());
	tsctimer.stop_timing("own pcg solver", true, "");

	std::cout << "Vector p after optimized solver: ";
	for(int i = 0; i < p.size() % 40; i++) std::cout << p_vec[i] << " ";
	std::cout << std::endl;

}
