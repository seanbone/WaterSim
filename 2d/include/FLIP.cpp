#include "FLIP.h"

FLIP::FLIP(Particle* particles, Mac2d& MACGrid) {
	
}

/**
 * Advance FLIP simulation by one frame
 */
void FLIP::step_FLIP(double dt, double time, unsigned long step) {
	/** One FLIP step:
	 * 1. Forward Euler (/RK2) to update particle positions
	 * 2. Transfer particle velocities to grid $\rightarrow u^*$ & save copy
	 * 3. Apply external forces (gravity etc) with Forward Euler: $u^* += \vec g\Delta t$
	 * 4. Construct pressure matrix $A$
	 * 5. Compute right-hand side (divergence) $d$
	 * 6. Solve $Ap = d$ with MICCG(0)
	 * 7. Update grid velocities with pressure gradients $\rightarrow u^{n+1}$
	 * 8. Update particle velocities by mixing FLIP and PIC
	 */

	// 1.
	advance_particles();

	// 2.
	particle_to_grid();

	// 3.
	apply_forces();

	// 4. to 6.
	update_pressures();

	// 7.
	apply_pressure_gradients();

	// 8.
	grid_to_particle();

}

void FLIP::advance_particles() {
	// TODO: forward euler to update particle positions
}

void FLIP::particle_to_grid() {
	// TODO: particle to grid transfer of velocities
}

void FLIP::apply_forces() {
	// TODO: apply external forces to velocities on grid via forward euler
}

void FLIP::update_pressures() {
	// TODO: construct and solve LSE for new pressures on grid
}

void FLIP::apply_pressure_gradients() {
	// TODO: update velocities on grid 
}

void FLIP::grid_to_particle() {
	// TODO: FLIP grid to particle transfer
}

//void FLIP::fwd_Euler( Eigen::Vector3d& velocity,
//					  const Eigen::Vector3d& ext_forces,
//					  const double dt)
//{
//	//~ TODO
//}

//~ void FLIP::transfer_Velocities(){}

//~ void FLIP::construct_A(){}

//~ void FLIP::compute_Rhs(){}

//~ void FLIP::solve_MICCG(){}

//~ void FLIP::update_GridVelocities(){}

//~ void FLIP::update_ParticleVelocities(){}
