/**
 * This file declares the data structure housing the FLIP particles.
 */

#ifndef WATERSIM_PARTICLES_H
#define WATERSIM_PARTICLES_H

#include "Mac3d.h"
#include "SimConfig.h"

struct Particles {
	/**
	 * Index type for particles.
	 */
	using particleIdx_t = unsigned int;

	// Positions
	double *x = nullptr;
	double *y = nullptr;
	double *z = nullptr;

	// Velocities
	double *u = nullptr;
	double *v = nullptr;
	double *w = nullptr;

	Particles(Particles::particleIdx_t nParticles, const Mac3d &macGrid);

	~Particles();

	/**
	 * Get the indices of the cell the particle is in
	 */
	inline void get_cell_index(particleIdx_t particleIdx, Mac3d::cellIdx_t &cellIdxX, Mac3d::cellIdx_t &cellIdxY,
	                           Mac3d::cellIdx_t &cellIdxZ) const {
		cellIdxX = (Mac3d::cellIdx_t) (x[particleIdx] * rcell_size_x_ + 0.5);
		cellIdxY = (Mac3d::cellIdx_t) (y[particleIdx] * rcell_size_y_ + 0.5);
		cellIdxZ = (Mac3d::cellIdx_t) (z[particleIdx] * rcell_size_z_ + 0.5);
	}

	/**
	 * Returns the number of particles.
	 */
	inline particleIdx_t get_num_particles() const { return num_particles_; }

private:

	//! 1 / cell_size_xyz
	double rcell_size_x_;
	double rcell_size_y_;
	double rcell_size_z_;

	//! Number of particles
	particleIdx_t num_particles_;
};


#endif //WATERSIM_PARTICLES_H
