/**
 * This file contains the implementation for the apply_boundary_conditions FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/*** BOUNDARY CONDITIONS ***/
void FLIP::apply_boundary_conditions() {

    // Get total number of cells on each axis
    unsigned nx = MACGrid_->get_num_cells_x();
    unsigned ny = MACGrid_->get_num_cells_y();
    unsigned nz = MACGrid_->get_num_cells_z();

    // Enforce boundary conditions for outer (system) boundaries
    for (unsigned k = 0; k < nz; k++) {
        for (unsigned i = 0; i < nx; i++) {
            MACGrid_->set_v(i, 0, k, 0);
            MACGrid_->set_v(i, ny, k, 0);
        }
    }

    for (unsigned k = 0; k < nz; k++) {
        for (unsigned j = 0; j < ny; j++) {
            MACGrid_->set_u(0, j, k, 0);
            MACGrid_->set_u(nx, j, k, 0);
        }
    }

    for (unsigned j = 0; j < ny; j++) {
        for (unsigned i = 0; i < nx; i++) {
            MACGrid_->set_w(i, j, 0, 0);
            MACGrid_->set_w(i, j, nz, 0);
        }
    }
}
