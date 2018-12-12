#ifndef FLIP_H
#define FLIP_H

#include "Particle.h"
#include "Mac2d.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
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
	FLIP(Particle* particles, const unsigned num_particles, Mac2d* MACGrid,
		 const double density, const double gravity, const double alpha);
	
	/** Perform one FLIP step
	 * Params:
	 * - dt is the amount of time to advance the sim by (may be subsampled)
	 * - time is the current time of the simulation
	 * - step is the number of steps performed
	 */
	void step_FLIP(const double dt, const unsigned long step);
	
	private:
	
	// Array of all particles
	Particle* particles_;
	const unsigned num_particles_;
	const double fluid_density_;
	const double gravity_mag_;

	// FLIP: alpha = 0.
	// PIC: alpha = 1.
	const double alpha_;
	
	// MAC Grid
	Mac2d* MACGrid_;
	
	// Pressure-Matrix
	SparseMat_t A_;
	// RHS of pressure LSE
	Eigen::VectorXd d_;
	
	// Compute timestep to satisfy CFL condition
	double compute_timestep( const double dt );

	// Compute velocity field by particle-to-grid transfer
	//   and extrapolating into air region
	void compute_velocity_field();
	
	// Check the distance between a particle and a point on the MACGrid
	bool check_threshold( const Eigen::Vector3d& particle_coord, 
						  const Eigen::Vector3d& grid_coord, 
						  const double h );
	
	// Compute the weight using the SPH Kernels multiplied by
	// the norm ||x_p - x_uij||
	double compute_weight( const Eigen::Vector3d& particle_coord, 
						const Eigen::Vector3d& grid_coord, 
						const double h );
	
	// Accumulate velocities and weights for u					
	void accumulate_u( const Eigen::Vector3d& pos,
					   const Eigen::Vector3d& vel,
					   const Eigen::Vector3d& grid_coord,
					   const double h,
					   const int i,
					   const int j );
	
	// Accumulate velocities and weights for v
	void accumulate_v( const Eigen::Vector3d& pos,
					   const Eigen::Vector3d& vel,
					   const Eigen::Vector3d& grid_coord,
					   const double h,
					   const int i,
					   const int j );
					   
	// Normalize accumulated horizontal velocities
	void normalize_accumulated_u( bool* const visited_u );
	
	// Normalize accumulated vertical velocities
	void normalize_accumulated_v( bool* const visited_v );
	
	// Extrapolate horizontal velocities into air cells
	void extrapolate_u( const bool* const visited_u );
	
	// Extrapolate vertical velocities into air cells
	void extrapolate_v( const bool* const visited_v );
	
	// Apply external forces to velocities on grid
	void apply_forces(const double dt);

	// Enforce boundary conditions for grid & solid boundaries
	void apply_boundary_conditions();

	// Compute and apply pressures on velocity field
	void do_pressures(const double dt);
	void compute_pressure_matrix();
	void compute_pressure_rhs(const double dt);
	void apply_pressure_gradients(const double dt);

	// transfer grid velocities to particles
	void grid_to_particle();

	// Apply forward euler or RK2 to update particle positions
	void advance_particles(const double dt, const unsigned step);
	
	void explode(const double dt, const unsigned long step, const double value);
};

#endif
