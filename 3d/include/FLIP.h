#ifndef FLIP_H
#define FLIP_H

#include "Particle.h"
#include "Mac3d.h"
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
	FLIP(Particle* particles, const unsigned num_particles, Mac3d* MACGrid,
		 const double density, const double gravity, const double alpha);
	
	/** Perform one FLIP step
	 * Params:
	 * - dt is the amount of time to advance the simulation by (may be 
	 * 	 subsampled)
	 * - step is the number of steps performed
	 */
	void step_FLIP(const double dt, const unsigned long step);
	
	private:
	
	// Array of all particles
	Particle* particles_;
	
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
	
	/** Compute timestep to satisfy CFL condition
	 * Params:
	 * - dt is the amount of time to advance the simulation by, and to 
	 * 	 subsample
	 */
	double compute_timestep( const double dt );

	/** Compute velocity field by particle-to-grid transfer
	 *  and extrapolate into air region
	 */
	void compute_velocity_field();
	
	/** Check the distance between a particle and a point on the MACGrid
	 * Params:
	 * - particle_coord is an Eigen vector containing the coordinates of 
	 * 	 a particle
	 * - grid_coord is an Eigen vector containing the coordinates of 
	 * 	 the center of a face on which a grid-velocity is saved
	 * - h is the distance threshold to be satisfied
	 */
	bool check_threshold( const Eigen::Vector3d& particle_coord, 
						  const Eigen::Vector3d& grid_coord, 
						  const double h );
	
	/** Compute the weight using the SPH Kernels multiplied by the norm 
	 * 	||x_p - x_uij||
	 * Params:
	 * - particle_coord is an Eigen vector containing the coordinates of 
	 * 	 a particle
	 * - grid_coord is an Eigen vector containing the coordinates of 
	 * 	 the center of a face on which a grid-velocity is saved
	 * - h is the distance threshold to be satisfied
	 */
	double compute_weight( const Eigen::Vector3d& particle_coord, 
						const Eigen::Vector3d& grid_coord, 
						const double h );
	
	/** Accumulate velocities and weights for u	
	 * Params:
	 * - pos is an Eigen vector containing the coordinates of a particle
	 * - vel is an Eigen vector containing the velocity components of a 
	 * 	 particle
	 * - grid_coord is an Eigen vector containing the coordinates of 
	 * 	 the center of a face on which a grid-velocity is saved
	 * - h is the distance threshold to be satisfied by check_threshold
	 * - i, j and k are the indices corresponding to a grid-cell
	 */		
	void accumulate_u( const Eigen::Vector3d& pos,
					   const Eigen::Vector3d& vel,
					   const Eigen::Vector3d& grid_coord,
					   const double h,
					   const int i,
					   const int j,
					   const int k );
	
	/** Accumulate velocities and weights for v	
	 * Params:
	 * - pos is an Eigen vector containing the coordinates of a particle
	 * - vel is an Eigen vector containing the velocity components of a 
	 * 	 particle
	 * - grid_coord is an Eigen vector containing the coordinates of 
	 * 	 the center of a face on which a grid-velocity is saved
	 * - h is the distance threshold to be satisfied by check_threshold
	 * - i, j and k are the indices corresponding to a grid-cell
	 */
	void accumulate_v( const Eigen::Vector3d& pos,
					   const Eigen::Vector3d& vel,
					   const Eigen::Vector3d& grid_coord,
					   const double h,
					   const int i,
					   const int j, 
					   const int k );
					   
	/** Accumulate velocities and weights for w	
	 * Params:
	 * - pos is an Eigen vector containing the coordinates of a particle
	 * - vel is an Eigen vector containing the velocity components of a 
	 * 	 particle
	 * - grid_coord is an Eigen vector containing the coordinates of 
	 * 	 the center of a face on which a grid-velocity is saved
	 * - h is the distance threshold to be satisfied by check_threshold
	 * - i, j and k are the indices corresponding to a grid-cell
	 */
	void accumulate_w( const Eigen::Vector3d& pos,
					   const Eigen::Vector3d& vel,
					   const Eigen::Vector3d& grid_coord,
					   const double h,
					   const int i,
					   const int j, 
					   const int k );
					   
	/** Normalize accumulated horizontal velocities
	 * Params:
	 * - visited_u is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from compute_velocity_field
	 */
	void normalize_accumulated_u( bool* const visited_u );
	
	/** Normalize accumulated vertical velocities
	 * Params:
	 * - visited_v is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from compute_velocity_field
	 */
	void normalize_accumulated_v( bool* const visited_v );
	
	/** Normalize accumulated outgoing velocities
	 * Params:
	 * - visited_w is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from compute_velocity_field
	 */
	void normalize_accumulated_w( bool* const visited_w );
	
	/** Extrapolate horizontal velocities into air cells
	 * Params:
	 * - visited_u is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from compute_velocity_field
	 */
	void extrapolate_u( const bool* const visited_u );
	
	/** Extrapolate vertical velocities into air cells
	 * Params:
	 * - visited_v is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from compute_velocity_field
	 */
	void extrapolate_v( const bool* const visited_v );
	
	/** Extrapolate outgoing velocities into air cells
	 * * Params:
	 * - visited_w is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from compute_velocity_field
	 */
	void extrapolate_w( const bool* const visited_w );
	
	/** Apply external forces to velocities on grid
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 */
	void apply_forces(const double dt);

	/** Enforce boundary conditions for grid & solid boundaries
	 */
	void apply_boundary_conditions();

	/** Compute and apply pressures on velocity field
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 */
	void do_pressures(const double dt);
	void compute_pressure_matrix();
	void compute_pressure_rhs(const double dt);
	void apply_pressure_gradients(const double dt);

	/** Transfer grid velocities to particles
	 */
	void grid_to_particle();

	/** Apply RK2 to update particle positions
	 * Params:
	 * - dt is the amount of time to advance the simulation by
	 * - step is the number of steps performed
	 */
	void advance_particles(const double dt, const unsigned long step);
	
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
};

#endif
