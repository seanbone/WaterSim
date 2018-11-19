#include "FLIP.h"

FLIP::FLIP(Particle* particles, const unsigned num_particles, Mac2d* MACGrid) 
	: particles_(particles), num_particles_(num_particles), MACGrid_(MACGrid) {
	
}

/**
 * Advance FLIP simulation by one frame
 */
void FLIP::step_FLIP(const double dt, const double time, const unsigned long step) {
	/** One FLIP step:
	 * 1. Compute velocity field (particle-to-grid transfer)
	 *    - Particle-to-grid transfer
	 *    - Classify cells (fluid/air)
	 *    - Extrapolate velocity field into air region
	 * 2. Apply external forces (fwd euler on field)
	 * 3. Compute & apply pressure gradients
	 * 4. Update particle velocities
	 * 8. Update particle positions
	 */

	// TODO: subsample time interval to satisfy CFL condition
	double cur_time = time;

	// 1.
	compute_velocity_field();

	// 2.
	apply_forces();

	// 3.
	do_pressures();

	// 4.
	grid_to_particle();

	// 5.
	advance_particles();
}

void FLIP::compute_velocity_field() {
	// TODO: 1. Compute the velocity field (velocities on grid)
	//  1a. particle-to-grid transfer
	//  1b. classify nonsolid cells as fluid or air
	//  1c. extrapolate velocity field to air cells
	//    -> see SIGGRAPH §6.3
}

void FLIP::apply_forces() {
	// TODO: compute&apply external forces (gravity, vorticity confinement, ...) 
	// Apply them to the velocity field via forward euler
	// Only worry about gravity for now
}

void FLIP::do_pressures() {
	// TODO: 3. Compute & apply pressure gradients to field
	//   3a. Compute A matrix
	//   3b. Compute rhs d
	//   3c. Solve for p: Ap = d (MICCG(0))
	//   3d. Apply pressure gradients to velocity field
	//     -> see SIGGRAPH §4
	// Note: boundary conditions are handles here by setting the pressures
	//    such that no particle exits the system.
}

void FLIP::grid_to_particle() {
	// TODO: FLIP grid to particle transfer
	//  -> See slides Fluids II, FLIP_explained.pdf
}

void FLIP::advance_particles() {
	// TODO: update particle positions 
	//  - Use RK2 interpolator
}

