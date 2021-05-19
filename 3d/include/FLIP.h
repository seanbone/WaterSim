#ifndef FLIP_H
#define FLIP_H

#include "Mac3d.h"
#include "SimConfig.h"
#include "NcWriter.h"
#include "Particles.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm> // std::copy
#include <Eigen/IterativeLinearSolvers> // solve sparse systems
#include <chrono>

class FLIP {
public:
	
	using SparseMat_t = Eigen::SparseMatrix<double>;
	using Triplet_t = Eigen::Triplet<double>;

	/** Constructor
	 * Params:
	 * - particles is a list of length num_particles containing all 
	 * 	 particles
	 *   to be used for the simulation
	 * - num_particles is the size of particles list
	 * - MACGrid is an initialized MAC grid
	 * - density is the density of the fluid
	 * - gravity is the value of a force directed downwards since the 
	 * 	 beginning of the simulation
	 * - alpha is the percentage of PIC method mixed with the FLIP 
	 *   method
	 */
	// TODO: remove old particle structure and num_particles params
	FLIP(Particles& particles, Mac3d* MACGrid, const SimConfig& cfg);

    /** Destructor
	 */	
	~FLIP();

	/** Perform one FLIP step
	 * Params:
	 * - dt is the amount of time to advance the simulation by (may be 
	 * 	 subsampled)
	 * - step is the number of steps performed
	 */
	void step_FLIP(double dt, unsigned long step);

	/** Compute velocity field by particle-to-grid transfer
	 *  and extrapolate into air region
	 */
	void particle_to_grid();

	/** Apply external forces to velocities on grid
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 */
	void apply_forces(const double dt);

	/** Simulate explosion of the meteorite
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 * - step is the number of steps performed
	 * - x, y and z are the indices of the grid-cell on which the impact
	 *	 will take place
	 * - r is the radius (in number of cells) of the meteorite
	 * - value is the amount of force applied by the meteorite on the 
	 * 	 grid-cells
	 */
	void explode(const double dt, const unsigned long step, const int x, const int y, const int z, const double r, const double value);

	/** Enforce boundary conditions for grid & solid boundaries
	 */
	void apply_boundary_conditions();

	/** Compute and apply pressures on velocity field
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 */
	void apply_pressure_correction(const double dt);

	/** Transfer grid velocities to particles
	 */
	void grid_to_particle();

	/** Compute timestep to satisfy CFL condition
	 * Params:
	 * - dt is the amount of time to advance the simulation by, and to 
	 * 	 subsample
	 */
	double compute_timestep( const double dt );

	/** Apply RK2 to update particle positions
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 * - step is the number of steps performed
	 */
	void advance_particles(const double dt, const unsigned long step);

private:

	// The configuration
	const SimConfig cfg_;

	// NcWriter object to write reference data
	NcWriter* ncWriter_;
	
	// FLIP particles
	Particles& particles_;
	
	// Total number of particles
	const unsigned num_particles_;
	
	// Density of the fluid to be simulated
	const double fluid_density_;
	
	// Value of the gravitational force
	const double gravity_mag_;

	// FLIP: alpha = 0.
	// PIC: alpha = 1.
	const double alpha_;
	
	// Pointer to MAC Grid
	Mac3d* MACGrid_;
	
	// Pressure-Matrix
	SparseMat_t A_;
	
	// RHS of pressure LSE
	Eigen::VectorXd d_;

	void compute_pressure_matrix();
	void compute_pressure_rhs(const double dt);
	void apply_pressure_gradients(const double dt);
};

#endif
