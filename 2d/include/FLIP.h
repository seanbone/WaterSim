#ifndef FLIP_H
#define FLIP_H

#include "Particle.h"
#include "Mac2d.h"
#include <Eigen/Sparse>
#include <algorithm> // std::copy
#include <Eigen/IterativeLinearSolvers> // solve sparse systems

class FLIP {
	public:
	
	using SparseMat_t = Eigen::SparseMatrix<double>;
	using Triplet_t = Eigen::Triplet<double>;

	/** Constructor
	 * Params:
	 * - particles is a list of length num_particles containing all particles
	 *    to be used for the simulation
	 * - num_particles is the size of particles list
	 * - MACGrid is an initialized MAC grid
	 */
	FLIP(Particle* particles, const unsigned num_particles, Mac2d* MACGrid);
	
	/** Perform one FLIP step
	 * Params:
	 * - dt is the amount of time to advance the sim by (may be subsampled)
	 * - time is the current time of the simulation
	 * - step is the number of steps performed
	 */
	void step_FLIP(const double dt, const double time, const unsigned long step);
	
	private:
	
	// Array of all particles
	Particle* particles_;
	const unsigned num_particles_;
	// TODO: pass as sim parameter
	const double fluid_density_ = 1.;
	
	// MAC Grid
	Mac2d* MACGrid_;
	
	// Pressure-Matrix
	SparseMat_t A_;
	// RHS of pressure LSE
	Eigen::VectorXd d_;

	// Compute velocity field by particle-to-grid transfer
	//   and extrapolating into air region
	void compute_velocity_field();

	// Apply external forces to velocities on grid
	void apply_forces();

	// Compute and apply pressures on velocity field
	void do_pressures(const double dt);
	void compute_pressure_matrix();
	void compute_pressure_rhs(const double dt);
	void apply_pressure_gradients(const double dt);

	// transfer grid velocities to particles
	void grid_to_particle();

	// Apply forward euler or RK2 to update particle positions
	void advance_particles();
};

#endif
