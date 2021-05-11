/**
 * This file contains the implementation for the grid_to_particle FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"


/*** UPDATE PARTICLE VELOCITIES ***/
void FLIP::grid_to_particle(){

    // FLIP grid to particle transfer
    //  -> See slides Fluids II, FLIP_explained.pdf

    // FLIP: alpha = 0.
    // PIC: alpha = 1.
    double alpha = alpha_;

    // Get total number of cells on each axis
    int nx = MACGrid_->get_num_cells_x();
    int ny = MACGrid_->get_num_cells_y();
    int nz = MACGrid_->get_num_cells_z();

    // Iterate over all particles
    for(unsigned i = 0; i < num_particles_; ++i){

        // Store the initial positions and velocities of the particles
        Eigen::Vector3d initial_position = (particlesOLD_ + i)->get_position();
        Eigen::Vector3d initial_velocity = (particlesOLD_ + i)->get_velocity();

        // Get the index of the grid-cell containing the current
        // particle
        auto initial_idx = MACGrid_->index_from_coord(initial_position(0),
                                                      initial_position(1),
                                                      initial_position(2));

        // Initialization of variables
        Eigen::Vector3d interp_u_star;
        interp_u_star.setZero();
        Eigen::Vector3d interp_u_n1;
        interp_u_n1.setZero();
        Eigen::Vector3d u_update;
        u_update.setZero();

        // Aliases for initial_position components
        double x = initial_position[0];
        double y = initial_position[1];
        double z = initial_position[2];

        // With u*, v* and w* we can do the interpolation by means of
        // interp(u*, x_p)
        // With the new u, v and w we can do the interpolation by means
        // of interp(u_n1, x_p)

        // Update the horizontal velocity (trilinear interpolation)
        interp_u_star[0] = MACGrid_->get_interp_u(x,y,z,true);
        interp_u_n1[0] = MACGrid_->get_interp_u(x,y,z);

        // Update the vertical velocity (trilinear interpolation)
        interp_u_star[1] = MACGrid_->get_interp_v(x,y,z,true);
        interp_u_n1[1] = MACGrid_->get_interp_v(x,y,z);

        // Update the outgoing velocity (trilinear interpolation)
        interp_u_star[2] = MACGrid_->get_interp_w(x,y,z,true);
        interp_u_n1[2] = MACGrid_->get_interp_w(x,y,z);

        // Blend PIC and FLIP, use double the amount of PIC on boundary
        if (initial_idx(0) == 0 or initial_idx(0) == nx-1
            or initial_idx(1) == 0 or initial_idx(1) == ny-1
            or initial_idx(2) == 0 or initial_idx(2) == nz-1){

            u_update = initial_velocity*(1 - std::min(1., 2*alpha)) + interp_u_n1 - interp_u_star*(1 - std::min(1., 2*alpha));

        } else {

            u_update = initial_velocity*(1 - alpha) + interp_u_n1 - interp_u_star*(1 - alpha);

        }

        // Finally, update the velocities of the particles
        (particlesOLD_ + i)->set_velocity(u_update);
    }
}

