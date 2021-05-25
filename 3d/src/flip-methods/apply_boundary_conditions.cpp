/**
 * This file contains the implementation for the apply_boundary_conditions FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/*** BOUNDARY CONDITIONS ***/
void FLIP::apply_boundary_conditions() {

	// Get total number of cells on each axis
	const Mac3d::cellIdx_t nx = MACGrid_->N_;
	const Mac3d::cellIdx_t ny = MACGrid_->M_;
	const Mac3d::cellIdx_t nz = MACGrid_->L_;

	// Enforce boundary conditions for outer (system) boundaries
	for(Mac3d::cellIdx_t k = 0; k < nz; ++k){
		for(Mac3d::cellIdx_t j = 0; j < ny; ++j){
			MACGrid_->pu_[     (nx+1) * (j + ny * k)] = 0.;
			MACGrid_->pu_[nx + (nx+1) * (j + ny * k)] = 0.;
		}
	}

	for(Mac3d::cellIdx_t k = 0; k < nz; ++k){
		for(Mac3d::cellIdx_t i = 0; i < nx; ++i){
			MACGrid_->pv_[i + nx *       (ny+1) * k ] = 0.;
			MACGrid_->pv_[i + nx * (ny + (ny+1) * k)] = 0.;
		}
	}

	for(Mac3d::cellIdx_t j = 0; j < ny; ++j){
		for(Mac3d::cellIdx_t i = 0; i < nx; ++i){
			MACGrid_->pw_[i + nx *  j           ] = 0.;
			MACGrid_->pw_[i + nx * (j + ny * nz)] = 0.;
		}
	}
}
