#ifndef FLIP_H
#define FLIP_H

#include "Mac3d.h"
#include "SimConfig.h"
#include "Particles.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm> // std::copy
#include <Eigen/IterativeLinearSolvers> // solve sparse systems
#include <chrono>
#include <vector>
#include "ConjugateGradient.hpp"
#include <immintrin.h>

#ifdef WRITE_REFERENCE
#include "NcWriter.h"
#endif

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
	 * - radius is the radius (in number of cells) of the meteorite
	 * - force is the amount of force applied by the meteorite on the 
	 * 	 grid-cells
	 */
	void explode(const double dt, const unsigned long step, const int x, const int y, const int z, const double radius, const double force);

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
	 */
	void advance_particles(const double dt);

private:

	// The configuration
	const SimConfig cfg_;

#ifdef WRITE_REFERENCE
	// NcWriter object to write reference data
	NcWriter* ncWriter_;
#endif
	
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
	double* d_;

	// Pressure vector
	std::vector<double> p;

	// Conjugate Gradient Solver
	ICConjugateGradientSolver cg_solver;
	
	/** Compute the weight using the SPH Kernels multiplied by the norm 
	 * 	||x_p - x_uij||
	 * Params:
	 * - r2 is the squared norm of x_p - x_uij
	 * - h2 is the squared value of h
	 * - h is the particle-grid interaction threshold in meters
	 */
	double compute_weight( const double r2,
						   const double h2,
						   const double h );
	
	/** Accumulate velocities and weights
	 * Params:
	 * - vel_grid is a pointer to the velocity field on the grid
	 * - weights is a pointer to an array of weights used to scale velocities 
	 *   during accumulation
	 * - vel_particle is the velocity of the current particle
	 * - x_, y_, z_particle are the coordinates of the current particle
	 * - face_coord_x, _y, _z are the coordinates on the current cell face where 
	 *   the velocities are stored
	 * - h is the particle-grid interaction threshold in meters
	 * - idx is the global index of the current gridcell
	 */
	void accumulate_vel( double* const vel_grid,
						 double* const weights,
						 const double vel_particle,
						 const double x_particle,
						 const double y_particle,
						 const double z_particle,
						 const double face_coord_x,
						 const double face_coord_y,
						 const double face_coord_z,
						 const double h,
						 const Mac3d::globalCellIdx_t idx );
	
	/** Normalize accumulated velocities
	 * Params:
	 * - visited_u is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from particle_to_grid
	 * - visited_v is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from particle_to_grid
	 * - visited_w is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from particle_to_grid
	 * - nx, ny, nz are the grid dimensions
	 */
	void normalize_accumulated_vels( bool* const visited_u, 
									 bool* const visited_v, 
									 bool* const visited_w, 
									 Mac3d::cellIdx_t nx, 
									 Mac3d::cellIdx_t ny, 
									 Mac3d::cellIdx_t nz );
	
	/** Extrapolate velocities into air cells
	 * Params:
	 * - vel is the velocity field on the grid
	 * - visited_vel is a lists of flags for visited grid-velocities: 
	 * 	 1 -> visited from particle_to_grid
	 * - n, m, l are the dimensions of the velocity fields (beware that 
	 *   velocities are stored on faces, thus a +1 is needed on the primary 
	 *   dimension, e.g. size(u) = (nx+1)*ny*nz -> n=nx+1, m=ny, l=nz)
	 */
	void extrapolate_vel( double* const vel,
						  const bool* const visited_vel,
						  const Mac3d::cellIdx_t n,
						  const Mac3d::cellIdx_t m,
						  const Mac3d::cellIdx_t l );

	void compute_pressure_matrix();
	void compute_pressure_rhs(const double dt);
	void apply_pressure_gradients(const double dt);
};

#endif
