/**
 * This file contains the implementation for the apply_forces FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/*** APPLY EXTERNAL FORCES ***/
void FLIP::apply_forces(const double dt) {

    // Compute&apply external forces (gravity, vorticity confinement, ...)
    // Apply them to the velocity field via forward euler
    // Only worry about gravity for now

    // Alias for MAC Grid
    auto& g = MACGrid_;

    // Get total number of cells on each axis
    const unsigned N = g->get_num_cells_x();
    const unsigned M = g->get_num_cells_y();
    const unsigned L = g->get_num_cells_z();

    // Iterate over cells & update: dv = dt*g
    for(unsigned k = 0; k < L; ++k){
        for (unsigned j = 0; j <= M; ++j){
            for (unsigned i = 0; i < N; ++i){
                g->set_v(i, j, k, g->get_v(i,j,k) - dt*gravity_mag_);
            }
        }
    }
}


void FLIP::explode(const double dt, const unsigned long step, const int x, const int y, const int z, const double r, const double value){

    // Apply external forces to simulate meteorite crash
    const double force = value;

    // Set the slope of the fall
    const double slope_x = 1;
    const double slope_y = 1;
    const double slope_z = 1;

    // Move the forces based on the current step
    const double sty = slope_y*step;
    const double stx = slope_x*sty;
    const double stz = slope_z*sty;

    // Get total number of cells on each axis
    const int nx = MACGrid_->get_num_cells_x();
    const int ny = MACGrid_->get_num_cells_y();
    const int nz = MACGrid_->get_num_cells_z();

    // Indices of the grid-cell at the center of the "meteorite"
    const int x_center = x + slope_x*nx/2 - stx;
    const int y_center = y + slope_y*ny/2 - sty;
    const int z_center = z + slope_z*nz/2 - stz;

    // Set the radius of the "meteorite"
    const double radius = r;

    // Iterate over all grid-velocities. If they are within the
    // radius r from the center, update the velocties with the given
    // force value (the force always points outwards)
    for( int k = 0; k < nz; ++k ){
        for( int j = 0; j < ny; ++j ){
            for( int i = 0; i < nx; ++i ){

                if (std::abs(i-x_center) <= radius and std::abs(j-y_center) <= radius and std::abs(k-z_center) <= radius){

                    if (i-x_center == 0 or std::signbit(i-x_center)){
                        MACGrid_->set_u(i, j, k, MACGrid_->get_u(i, j, k) - dt*force*(2.+(i-x_center)/radius));
                    } else {
                        MACGrid_->set_u(i, j, k, MACGrid_->get_u(i, j, k) + dt*force*(2.-(i-x_center)/radius));
                    }

                    if (j-y_center == 0 or std::signbit(j-y_center)){
                        MACGrid_->set_v(i, j, k, MACGrid_->get_v(i, j, k) - dt*force*(2.+(j-y_center)/radius));
                    } else {
                        MACGrid_->set_v(i, j, k, MACGrid_->get_v(i, j, k) + dt*force*(2.-(j-y_center)/radius));
                    }

                    if (k-z_center == 0 or std::signbit(k-z_center)){
                        MACGrid_->set_w(i, j, k, MACGrid_->get_w(i, j, k) - dt*force*(2.+(k-z_center)/radius));
                    } else {
                        MACGrid_->set_w(i, j, k, MACGrid_->get_w(i, j, k) + dt*force*(2.-(k-z_center)/radius));
                    }
                }
            }
        }
    }
}
