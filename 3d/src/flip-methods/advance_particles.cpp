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
		vel = (particlesOLD_ + n)->get_velocity();
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

	// Half timestep
	double dt_half = 0.5 * dt;

	// Get the lower and upper bounds of the grid (in meters)
	double x_lower_bound = -0.5 * MACGrid_->cell_sizex_;
	double y_lower_bound = -0.5 * MACGrid_->cell_sizey_;
	double z_lower_bound = -0.5 * MACGrid_->cell_sizez_;
	double x_upper_bound = MACGrid_->sizex_ + x_lower_bound;
	double y_upper_bound = MACGrid_->sizey_ + y_lower_bound;
	double z_upper_bound = MACGrid_->sizez_ + z_lower_bound;

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
		x_next = x_particle + dt*MACGrid_->get_interp_u(x_half, y_half, z_half);
		y_next = y_particle + dt*MACGrid_->get_interp_v(x_half, y_half, z_half);
		z_next = z_particle + dt*MACGrid_->get_interp_w(x_half, y_half, z_half);

		// Check if the particle exits the grid
		if(x_next <= x_lower_bound) x_next = 0.;
		if(y_next <= y_lower_bound) y_next = 0.;
		if(z_next <= z_lower_bound) z_next = 0.;
		if(x_next >= x_upper_bound) x_next = MACGrid_->sizex_ - MACGrid_->cell_sizex_; 
		if(y_next >= y_upper_bound) y_next = MACGrid_->sizey_ - MACGrid_->cell_sizey_;
		if(z_next >= z_upper_bound) z_next = MACGrid_->sizez_ - MACGrid_->cell_sizez_;

		// Check if the particle enters in a solid

		// Get the indices of the grid-cells containing the particle at
		// the current time and in the future
		particles_.get_cell_index(n, curr_cell_idx_x, curr_cell_idx_y, curr_cell_idx_z);
		MACGrid_->index_from_coord(x_next, y_next, z_next, next_cell_idx_x, next_cell_idx_y, next_cell_idx_z);

		// Shift a particle if it would exit the system
		// (should not happen)
		if( MACGrid_->is_solid(next_cell_idx_x, next_cell_idx_y, next_cell_idx_z) ) {

			if     ( curr_cell_idx_x > next_cell_idx_x ) x_next = (curr_cell_idx_x - 0.25) * MACGrid_->cell_sizex_;
			else if( curr_cell_idx_x < next_cell_idx_x ) x_next = (curr_cell_idx_x + 0.25) * MACGrid_->cell_sizex_;

			if     ( curr_cell_idx_y > next_cell_idx_y ) y_next = (curr_cell_idx_y - 0.25) * MACGrid_->cell_sizey_;
			else if( curr_cell_idx_y < next_cell_idx_y ) y_next = (curr_cell_idx_y + 0.25) * MACGrid_->cell_sizey_;

			if     ( curr_cell_idx_z > next_cell_idx_z ) z_next = (curr_cell_idx_z - 0.25) * MACGrid_->cell_sizez_;
			else if( curr_cell_idx_z < next_cell_idx_z ) z_next = (curr_cell_idx_z + 0.25) * MACGrid_->cell_sizez_;
		}

		// Update the position of the current particle
		particles_.x[n] = x_next;
		particles_.y[n] = y_next;
		particles_.z[n] = z_next;
	}
}

