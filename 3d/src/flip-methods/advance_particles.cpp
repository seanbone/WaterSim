/**
 * This file contains the implementation for the advance_particles FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/** COMPUTE ADVECTION TIMESTEP BASED ON CFL CONDITIONS */
double FLIP::compute_timestep( const double dt ){

	// Particle velocity
	double u_particle;
	double v_particle;
	double w_particle;

	// Maximum particle velocity
	double u_max = 0.;
	double v_max = 0.;
	double w_max = 0.;

	// Get the maximal particle velocity components
	for( Particles::particleIdx_t n = 0; n < num_particles_; ++n ){

		u_particle = particles_.u[n];
		v_particle = particles_.v[n];
		w_particle = particles_.w[n];
		
		if( std::abs(u_particle) > u_max ) u_max = std::abs(u_particle);
		if( std::abs(v_particle) > v_max ) v_max = std::abs(u_particle);
		if( std::abs(w_particle) > w_max ) w_max = std::abs(w_particle);
	}

	// Check if the fastest particles travel a distance larger than the
	// length of an edge of a cell

	// New timestep that satisfies CFL condition
	double tmp;
	double dt_new = dt;

	if( u_max != 0. ){
		tmp = MACGrid_->cell_sizex_/u_max;
		if( tmp < dt ) dt_new = tmp;
	}

	if( v_max != 0. ){
		tmp = MACGrid_->cell_sizey_/v_max;
		if( tmp < dt_new ) dt_new = tmp;
	}

	if( w_max != 0. ){
		tmp = MACGrid_->cell_sizez_/w_max;
		if( tmp < dt_new ) dt_new = tmp;
	}

	return dt_new;
}


/** PARTICLE ADVECTION */
void FLIP::advance_particles(const double dt) {

	// Grid dimensions
	const Mac3d::cellIdx_t nx = MACGrid_->N_;
	const Mac3d::cellIdx_t ny = MACGrid_->M_;

	// Half timestep
	const double dt_half = 0.5 * dt;

	// Get the lower and upper bounds of the grid (in meters)
	const double x_lower_bound = -0.5 * MACGrid_->cell_sizex_;
	const double y_lower_bound = -0.5 * MACGrid_->cell_sizey_;
	const double z_lower_bound = -0.5 * MACGrid_->cell_sizez_;
	const double x_upper_bound = MACGrid_->sizex_ + x_lower_bound;
	const double y_upper_bound = MACGrid_->sizey_ + y_lower_bound;
	const double z_upper_bound = MACGrid_->sizez_ + z_lower_bound;

	const double x_last_center = MACGrid_->sizex_ - MACGrid_->cell_sizex_;
	const double y_last_center = MACGrid_->sizey_ - MACGrid_->cell_sizey_;
	const double z_last_center = MACGrid_->sizez_ - MACGrid_->cell_sizez_;

	// Position and velocity of the current particle
	double x_particle;
	double y_particle;
	double z_particle;

	// Particle coordinates after half timestep (computed with Euler)
	double x_half;
	double y_half;
	double z_half;

	// Coordinates of the future location of the particle
	double x_next;
	double y_next;
	double z_next;

	// Cell indices of the particle (before and after advection)
	Mac3d::cellIdx_t curr_cell_idx_x;
	Mac3d::cellIdx_t curr_cell_idx_y;
	Mac3d::cellIdx_t curr_cell_idx_z;
	Mac3d::cellIdx_t next_cell_idx_x;
	Mac3d::cellIdx_t next_cell_idx_y;
	Mac3d::cellIdx_t next_cell_idx_z;

	// Iterate over all particles
	for( Particles::particleIdx_t n = 0; n < num_particles_; ++n ){

		// Get current position and velocity of the particle
		x_particle = particles_.x[n];
		y_particle = particles_.y[n];
		z_particle = particles_.z[n];

		// Euler estimate
		x_half = x_particle + dt_half*particles_.u[n];
		y_half = y_particle + dt_half*particles_.v[n];
		z_half = z_particle + dt_half*particles_.w[n];

		// Check if the particle is out of the grid after the euler step
		if( x_half <= x_lower_bound or x_half >= x_upper_bound or 
			y_half <= y_lower_bound or y_half >= y_upper_bound or 
			z_half <= z_lower_bound or z_half >= z_upper_bound )
		{
			continue;
		}

		// RK2
		x_next = x_particle + dt*MACGrid_->grid_interpolate<Mac3d::GRID_U>(x_half, y_half, z_half);
		y_next = y_particle + dt*MACGrid_->grid_interpolate<Mac3d::GRID_V>(x_half, y_half, z_half);
		z_next = z_particle + dt*MACGrid_->grid_interpolate<Mac3d::GRID_W>(x_half, y_half, z_half);

		// Check if the particle exits the grid
		if(x_next <= x_lower_bound) x_next = 0.;
		if(y_next <= y_lower_bound) y_next = 0.;
		if(z_next <= z_lower_bound) z_next = 0.;
		if(x_next >= x_upper_bound) x_next = x_last_center; 
		if(y_next >= y_upper_bound) y_next = y_last_center;
		if(z_next >= z_upper_bound) z_next = z_last_center;

		// Check if the particle enters in a solid

		// Get the indices of the grid-cells containing the particle at
		// the current time and in the future
		//particles_.get_cell_index(n, curr_cell_idx_x, curr_cell_idx_y, curr_cell_idx_z);
		//MACGrid_->index_from_coord(x_next, y_next, z_next, next_cell_idx_x, next_cell_idx_y, next_cell_idx_z);

		//// Shift a particle if it would exit the system
		//// (should not happen)
		//if( MACGrid_->psolid_[next_cell_idx_x + nx*next_cell_idx_y + nx*ny*next_cell_idx_z] ) {

		//	if     ( curr_cell_idx_x > next_cell_idx_x ) x_next = (curr_cell_idx_x - 0.25) * MACGrid_->cell_sizex_;
		//	else if( curr_cell_idx_x < next_cell_idx_x ) x_next = (curr_cell_idx_x + 0.25) * MACGrid_->cell_sizex_;

		//	if     ( curr_cell_idx_y > next_cell_idx_y ) y_next = (curr_cell_idx_y - 0.25) * MACGrid_->cell_sizey_;
		//	else if( curr_cell_idx_y < next_cell_idx_y ) y_next = (curr_cell_idx_y + 0.25) * MACGrid_->cell_sizey_;

		//	if     ( curr_cell_idx_z > next_cell_idx_z ) z_next = (curr_cell_idx_z - 0.25) * MACGrid_->cell_sizez_;
		//	else if( curr_cell_idx_z < next_cell_idx_z ) z_next = (curr_cell_idx_z + 0.25) * MACGrid_->cell_sizez_;
		//}

		// Update the position of the current particle
		particles_.x[n] = x_next;
		particles_.y[n] = y_next;
		particles_.z[n] = z_next;
	}
}

