/**
 * This file contains the implementation for the apply_pressure_correction FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "ConjugateGradient.hpp"
#include "tsc_x86.hpp"


/*** PRESSURE SOLVING ***/
void FLIP::apply_pressure_correction(const double dt) {

    // Compute & apply pressure gradients to field

    // Compute rhs d
    compute_pressure_rhs(dt);

    // Solve for p: Ap = d (MICCG(0))
	// work directly on grid pressure array
	cg_solver.solve(d_, MACGrid_->ppressure_);

    // Apply pressure gradients to velocity field
    //     -> see SIGGRAPH §4
    apply_pressure_gradients(dt);
}

void FLIP::compute_pressure_rhs(const double dt) {

    // Compute right-hand side of the pressure equations and store in d_
    // See eq. (4.19) and (4.24) in SIGGRAPH notes
    // Note: u_{solid} = 0

    // Get total number of cells on each axis
    unsigned nx = MACGrid_->get_num_cells_x();
    unsigned ny = MACGrid_->get_num_cells_y();
    unsigned nz = MACGrid_->get_num_cells_z();

    // Alias for MAC Grid
    auto& g = MACGrid_;

    // Set d_ to zero
	std::fill(d_, d_+nx*ny*nz, 0);

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

                    d_[cellidx] = fluid_density_ * g->get_cell_sizex() * d_ij / dt;

                } else { // if is_fluid(i,j,k)

                    // Set the entry to zero if the current cell is a
                    // fluid-cell
                    d_[cellidx] = 0;
                }
            }
        }
    }
}


void FLIP::apply_pressure_gradients(const double dt) {

    // Apply pressure gradients to velocity field

    // Get total number of cells on each axis
    unsigned nx = MACGrid_->get_num_cells_x();
    unsigned ny = MACGrid_->get_num_cells_y();
    unsigned nz = MACGrid_->get_num_cells_z();

    // Alias for MAC Grid
    auto& g = MACGrid_;

    // Get the length of an edge
    double dx = g->get_cell_sizex();

    // Iterate over all grid-cells
    for (unsigned k = 0; k < nz; ++k) {
        for (unsigned j = 0; j < ny; ++j) {
            for (unsigned i = 0; i < nx; ++i) {

                // Update grid-velocities with new velocities induced by
                // pressures
                if (i != 0) {
                    // get_u(i,j,k) = u_{ (i-1/2, j, k) }
                    // See SIGGRAPH eq. (4.6)
                    double du = (g->get_pressure(i,j,k) - g->get_pressure(i-1,j,k));
                    du *= (dt/(dx*fluid_density_));

                    g->set_u(i,j,k, g->get_u(i,j,k) - du);
                }

                if (j != 0) {

                    // get_v(i,j,k) = v_{ (i, j-1/2, k) }
                    // See SIGGRAPH eq. (4.7)
                    double dv = (g->get_pressure(i,j,k) - g->get_pressure(i,j-1,k));
                    dv *= (dt/(dx*fluid_density_));

                    g->set_v(i,j,k, g->get_v(i,j,k) - dv);
                }

                if (k != 0) {

                    // get_w(i,j,k) = w_{ (i, j, k-1/2) }
                    // See SIGGRAPH eq. (4.8)

                    double dw = (g->get_pressure(i,j,k) - g->get_pressure(i,j,k-1));
                    dw *= (dt/(dx*fluid_density_));

                    g->set_w(i,j,k, g->get_w(i,j,k) - dw);
                }
            }
        }
    }
}
