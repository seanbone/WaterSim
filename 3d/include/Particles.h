/**
 * This file declares the data structure housing the FLIP particles.
 */

#ifndef WATERSIM_PARTICLES_H
#define WATERSIM_PARTICLES_H

#include "Mac3d.h"
#include "SimConfig.h"

struct Particles {
	using particleIdx_t = unsigned int;

	// Positions
	double *x = nullptr;
	double *y = nullptr;
	double *z = nullptr;

	// Velocities
	double *u = nullptr;
	double *v = nullptr;
	double *w = nullptr;

	Particles(Particles::particleIdx_t nParticles, const SimConfig &cfg, const Mac3d &macGrid);

	~Particles();

	/**
	 * Get the indices of the cell the particle is in
	 */
	inline void get_cell_index(particleIdx_t particleIdx, Mac3d::cellIdx_t &cellIdxX, Mac3d::cellIdx_t &cellIdxY,
	                           Mac3d::cellIdx_t &cellIdxZ) const;

	/** Returns the number of particles.
	 */
	particleIdx_t get_num_particles() const;

private:

	//! 1 / cell_size_xyz
	double rcell_size_x_;
	double rcell_size_y_;
	double rcell_size_z_;

	//! Number of particles
	particleIdx_t num_particles_;
};


#endif //WATERSIM_PARTICLES_H
