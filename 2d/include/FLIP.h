#ifndef FLIP_H
#define FLIP_H

#include "Particle.h"
#include "Mac2d.h"
#include <Eigen/Sparse>

using SparseMat_t = Eigen::SparseMatrix<double>;
using Triplet_t = Eigen::Triplet<double>;

class FLIP {
	public:
	
	// Default Constructor
	// TODO
	FLIP() {}
	
	// Copy Constructor
	//~ TODO
	
	// Init-Constructor
	FLIP( Particle* particles, Mac2d MACGrid, SparseMat_t& A );
	
	// Forward Euler
	// - Use velocity = particle.get_Velocties()
	// - For gravity use ext_forces = (0, 0, -9.81)
	void fwd_Euler( Eigen::Vector3d& velocity,
					const Eigen::Vector3d& ext_forces,
					const double dt);
	
	// Transfer particle velocities to grid
	// and save a copy of the grid velocities
	//~ TODO: void transfer_Velocities();
	
	// Construct Pressure-Matrix
	// - Reserve space?
	// - Use triplets
	//~ TODO: void construct_A();
	
	// Compute rhs of Ap = d
	//~ TODO: void compute_Rhs();
	
	// Solve Ap = d with MICCG(0)
	//~ TODO: void solve_MICCG();
	
	// Update grid velocities with pressure gradients
	//~ TODO: void update_GridVelocities();
	
	// Update particle velocities by mixing FLIP and PIC
	//~ TODO: void update_ParticleVelocities();
	
	// Perform one FLIP step
	//~ TODO: void step_FLIP();
	
	private:
	
	// Array of all particles
	Particle* particles_;
	
	// MAC Grid
	Mac2d MACGrid_;
	
	// Pressure-Matrix
	SparseMat_t A;
};

#endif
