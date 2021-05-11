/**
 * This file contains the implementation for the advance_particles FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/** COMPUTE ADVECTION TIMESTEP BASED ON CFL CONDITIONS */
double FLIP::compute_timestep( const double dt ){

    // New timestep that satisfies CFL condition
    double dt_new;

    // Particle velocity
    Eigen::Vector3d vel;

    // Get the maximal particle velocity components
    double u_max = 0;
    double v_max = 0;
    double w_max = 0;
    for( unsigned int n = 0; n < num_particles_; ++n ){
        vel = (particles_ + n)->get_velocity();
        if ( std::abs(vel(0)) > std::abs(u_max) ){
            u_max = vel(0);
        }
        if ( std::abs(vel(1)) > std::abs(v_max) ){
            v_max = vel(1);
        }
        if ( std::abs(vel(2)) > std::abs(w_max) ){
            w_max = vel(2);
        }
    }

    // Check if the fastest particles travel a distance larger than the
    // length of an edge of a cell
    if ( u_max == 0 ){
        dt_new = dt;
    } else {
        dt_new = std::abs(MACGrid_->get_cell_sizex()/u_max);
        if ( dt_new > dt){
            dt_new = dt;
        }
    }

    if ( v_max == 0 ){
        dt_new = dt;
    } else {
        double tmp = std::abs(MACGrid_->get_cell_sizey()/v_max);
        if ( tmp < dt_new){
            dt_new = tmp;
        }
        if ( dt_new > dt ){
            dt_new = dt;
        }
    }

    if ( w_max == 0 ){
        dt_new = dt;
    } else {
        double tmp = std::abs(MACGrid_->get_cell_sizez()/w_max);
        if ( tmp < dt_new){
            dt_new = tmp;
        }
        if ( dt_new > dt ){
            dt_new = dt;
        }
    }

    return dt_new;
}


/** PARTICLE ADVECTION */
void FLIP::advance_particles(const double dt, const unsigned long step) {

    // Iterate over all particles
    for(unsigned n = 0; n < num_particles_; ++n){

        // Get current position and velocity of the particle
        Eigen::Vector3d pos_curr = (particles_ + n)->get_position();
        Eigen::Vector3d vel = (particles_ + n)->get_velocity();

        // Coordinates of the future location of the particle
        Eigen::Vector3d pos_next;

        // Euler estimate
        Eigen::Vector3d pos_half = pos_curr + 0.5*dt*vel;

        // Get the size of the grid (in meters)
        double size_x = MACGrid_->get_grid_size()(0);
        double size_y = MACGrid_->get_grid_size()(1);
        double size_z = MACGrid_->get_grid_size()(2);

        // Get the length of the edges of a grid-cell (in meters)
        double cell_sizex = MACGrid_->get_cell_sizex();
        double cell_sizey = MACGrid_->get_cell_sizey();
        double cell_sizez = MACGrid_->get_cell_sizez();

        // Coordinates of intermediate position, computed with Euler
        double x_half = pos_half(0);
        double y_half = pos_half(1);
        double z_half = pos_half(2);

        // Check if pos_half is out of the grid
        if ((x_half <= -0.5*cell_sizex) or (x_half >= size_x - 0.5*cell_sizex)
            or (y_half <= -0.5*cell_sizey) or (y_half >= size_y - 0.5*cell_sizey)
            or (z_half <= -0.5*cell_sizez) or (z_half >= size_z - 0.5*cell_sizez)){
            continue;
        }

        // RK2
        pos_next(0) = pos_curr(0) + dt*MACGrid_->get_interp_u(x_half, y_half, z_half);
        pos_next(1) = pos_curr(1) + dt*MACGrid_->get_interp_v(x_half, y_half, z_half);
        pos_next(2) = pos_curr(2) + dt*MACGrid_->get_interp_w(x_half, y_half, z_half);

        // Aliases for pos_next coordinates
        double x = pos_next(0);
        double y = pos_next(1);
        double z = pos_next(2);

        // Check if the particle exits the grid
        if (x <= -0.5*cell_sizex) {
            pos_next(0) = 0.;
        }

        if (x >= size_x - 0.5*cell_sizex) {
            pos_next(0) = size_x - cell_sizex;
        }

        if (y <= -0.5*cell_sizey) {
            pos_next(1) = 0.;
        }

        if (y >= size_y - 0.5*cell_sizey) {
            pos_next(1) = size_y - cell_sizey;
        }

        if (z <= -0.5*cell_sizez) {
            pos_next(2) = 0.;
        }

        if (z >= size_z - 0.5*cell_sizez) {
            pos_next(2) = size_z - cell_sizez;
        }

        // Check if the particle enters in a solid

        // Get the indices of the grid-cells containing the particle at
        // the current time and in the future
        auto prev_indices = MACGrid_->index_from_coord(pos_curr(0), pos_curr(1), pos_curr(2));
        auto new_indices = MACGrid_->index_from_coord(pos_next(0), pos_next(1), pos_curr(2));

        // Get the length of the edges of a grid-cell (in meters)
        // [shorter names]
        double sx = MACGrid_->get_cell_sizex();
        double sy = MACGrid_->get_cell_sizey();
        double sz = MACGrid_->get_cell_sizez();

        // Shift a particle if it would exit the system
        // (should not happen)
        if (MACGrid_->is_solid(new_indices(0), new_indices(1), new_indices(2))) {

            if (prev_indices(0)  > new_indices(0))
                pos_next(0) = (prev_indices(0) - 0.25) * sx;
            else if (prev_indices(0) < new_indices(0))
                pos_next(0) = (prev_indices(0) + 0.25) * sx;

            if (prev_indices(1) > new_indices(1))
                pos_next(1) = (prev_indices(1) - 0.25) * sy;
            else if (prev_indices(1) < new_indices(1))
                pos_next(1) = (prev_indices(1) + 0.25) * sy;

            if (prev_indices(2) > new_indices(2))
                pos_next(2) = (prev_indices(2) - 0.25) * sz;
            else if (prev_indices(2) < new_indices(2))
                pos_next(2) = (prev_indices(2) + 0.25) * sz;
        }

        // Update the position of the current particle
        (particles_ + n)->set_position(pos_next);
    }
}

