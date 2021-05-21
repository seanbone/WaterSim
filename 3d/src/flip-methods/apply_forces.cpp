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

	// Get total number of faces on the y-axis
	const Mac3d::globalCellIdx_t n_vfaces = MACGrid_->N_ * MACGrid_->L_ * (MACGrid_->M_+1);

	const double dv = -dt * gravity_mag_;

	// Iterate over cells & update: dv = dt*g
	for(Mac3d::globalCellIdx_t idx = 0; idx < n_vfaces; ++idx){

		MACGrid_->pv_[idx] += dv;
	}
}


void FLIP::explode( const double dt, 
					const unsigned long step, 
					const Mac3d::cellIdx_t x, 
					const Mac3d::cellIdx_t y, 
					const Mac3d::cellIdx_t z, 
					const double radius, 
					const double force )
{

	// Set the slope of the fall
	const double slope_x = 1.;
	const double slope_y = 1.;
	const double slope_z = 1.;

	// Move the forces based on the current step
	const double sty = slope_y * step;
	const double stx = slope_x * sty;
	const double stz = slope_z * sty;

	// Get total number of cells on each axis
	const Mac3d::cellIdx_t nx = MACGrid_->N_;
	const Mac3d::cellIdx_t ny = MACGrid_->M_;
	const Mac3d::cellIdx_t nz = MACGrid_->L_;

	// Indices of the grid-cell at the center of the "meteorite"
	const Mac3d::cellIdx_t x_center = x + 0.5*slope_x*nx - stx;
	const Mac3d::cellIdx_t y_center = y + 0.5*slope_y*ny - sty;
	const Mac3d::cellIdx_t z_center = z + 0.5*slope_z*nz - stz;

	const double dvel      = dt*force;
	const double dvel2     = 2.*dvel;
	const double dvel_rinv = dvel/radius;
	const double du_center = x_center * dvel_rinv;
	const double dv_center = y_center * dvel_rinv;
	const double dw_center = z_center * dvel_rinv;

	// Iterate over all grid-velocities. If they are within the
	// radius r from the center, update the velocties with the given
	// force value (the force always points outwards)
	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		if( std::abs(k-z_center) <= radius ){
			for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
				if( std::abs(j-y_center) <= radius ){
					for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
						if(std::abs(i-x_center) <= radius){

							if( i-x_center == 0 or std::signbit(i-x_center) ){
								MACGrid_->pu_[i + (nx+1)*(j + ny*k)] -= dvel2 + i*dvel_rinv - du_center;
							} else {
								MACGrid_->pu_[i + (nx+1)*(j + ny*k)] += dvel2 - i*dvel_rinv + du_center;
							}

							if( j-y_center == 0 or std::signbit(j-y_center) ){
								MACGrid_->pv_[i + nx*(j + (ny+1)*k)] -= dvel2 + j*dvel_rinv - dv_center;
							} else {
								MACGrid_->pv_[i + nx*(j + (ny+1)*k)] += dvel2 - j*dvel_rinv + dv_center;
							}

							if( k-z_center == 0 or std::signbit(k-z_center) ){
								MACGrid_->pw_[i + nx*(j + ny*k)] -= dvel2 + k*dvel_rinv - dw_center;
							} else {
								MACGrid_->pw_[i + nx*(j + ny*k)] += dvel2 - k*dvel_rinv + dw_center;
							}
						}
					}
				}
			}
		}
	}

}
