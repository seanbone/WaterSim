#include "FLIP.h"

FLIP::FLIP( Particle* particles, Mac2d MACGrid, SparseMat_t& A ){
	
}

void FLIP::fwd_Euler( Eigen::Vector3d& velocity,
					  const Eigen::Vector3d& ext_forces,
					  const double dt)
{
	//~ TODO
}

//~ void FLIP::transfer_Velocities(){}

//~ void FLIP::construct_A(){}

//~ void FLIP::compute_Rhs(){}

//~ void FLIP::solve_MICCG(){}

//~ void FLIP::update_GridVelocities(){}

//~ void FLIP::update_ParticleVelocities(){}

//~ void FLIP::step_FLIP(){}
