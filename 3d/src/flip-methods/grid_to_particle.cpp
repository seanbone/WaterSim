/**
 * This file contains the implementation for the grid_to_particle FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"


/*** UPDATE PARTICLE VELOCITIES ***/
void FLIP::grid_to_particle(){

	// FLIP grid to particle transfer
	//  -> See slides Fluids II, FLIP_explained.pdf

	const double coeff_internal = 1. - alpha_;
	const double coeff_boundary = 1. - std::min(1., 2.*alpha_);

	// Get total number of cells on each axis
	const Mac3d::cellIdx_t nx = MACGrid_->get_num_cells_x();
	const Mac3d::cellIdx_t ny = MACGrid_->get_num_cells_y();
	const Mac3d::cellIdx_t nz = MACGrid_->get_num_cells_z();

	// Position and velocity of the current particle
	double x_particle;
	double y_particle;
	double z_particle;
	double u_particle;
	double v_particle;
	double w_particle;

	double interp_u_star;
	double interp_v_star;
	double interp_w_star;

	double interp_u_n1;
	double interp_v_n1;
	double interp_w_n1;

	double u_update;
	double v_update;
	double w_update;

	// Indices of the cell containing the current particle
	Mac3d::cellIdx_t cell_idx_x;
	Mac3d::cellIdx_t cell_idx_y;
	Mac3d::cellIdx_t cell_idx_z;

	// Iterate over all particles
	for( Particles::particleIdx_t n = 0; n < num_particles_; ++n ){

		// Store the initial positions and velocities of the particles
		x_particle = particles_.x[n];
		y_particle = particles_.y[n];
		z_particle = particles_.z[n];
		u_particle = particles_.u[n];
		v_particle = particles_.v[n];
		w_particle = particles_.w[n];

		// Get updated velocities by interpolation
		std::tie(interp_u_n1, interp_u_star) = MACGrid_->grid_interpolate<Mac3d::GRID_U, Mac3d::INTERPOLATE_BOTH>(x_particle, y_particle, z_particle);
		std::tie(interp_v_n1, interp_v_star) = MACGrid_->grid_interpolate<Mac3d::GRID_V, Mac3d::INTERPOLATE_BOTH>(x_particle, y_particle, z_particle);
		std::tie(interp_w_n1, interp_w_star) = MACGrid_->grid_interpolate<Mac3d::GRID_W, Mac3d::INTERPOLATE_BOTH>(x_particle, y_particle, z_particle);

		// Get the index of the grid-cell containing the current
		// particle
		particles_.get_cell_index(n, cell_idx_x, cell_idx_y, cell_idx_z);

		// Blend PIC and FLIP, use double the amount of PIC on boundary
		if( cell_idx_x == 0 or cell_idx_x == nx-1 or 
			cell_idx_y == 0 or cell_idx_y == ny-1 or 
			cell_idx_z == 0 or cell_idx_z == nz-1 )
		{

			u_update = interp_u_n1 + (u_particle - interp_u_star) * coeff_boundary;
			v_update = interp_v_n1 + (v_particle - interp_v_star) * coeff_boundary;
			w_update = interp_w_n1 + (w_particle - interp_w_star) * coeff_boundary;

		} else {

			u_update = interp_u_n1 + (u_particle - interp_u_star) * coeff_internal;
			v_update = interp_v_n1 + (v_particle - interp_v_star) * coeff_internal;
			w_update = interp_w_n1 + (w_particle - interp_w_star) * coeff_internal;

		}

		// Finally, update the velocities of the particles
		particles_.u[n] = u_update;
		particles_.v[n] = v_update;
		particles_.w[n] = w_update;
	}
}