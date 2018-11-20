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
	double alpha = 0.00;
	
	//Store the initial positions of the particles
	std::vector<Eigen::Vector3d> initial_positions(num_particles_);
	std::vector<Eigen::Vector3d> initial_velocities(num_particles_);
	for(int i = 0; i < num_particles_; ++i){
		initial_positions[i] = (particles_+i)->get_position();
		initial_velocities[i] = (particles_+i)->get_velocity();
	}
	
	//Applying the forces we can compute u* and already make the interpolation interp(u*, x_p)
	apply_forces();
	std::vector<Eigen::Vector3d> interp_u_star(num_particles_);
	for(int i = 0; i < num_particles_; ++i){
		Eigen::Vector3d temp1 = initial_positions[i];
		double x = temp1[0];
		double y = temp1[1];
		Mac2d::Pair_t indices = MACGrid_->index_from_coord(x,y);
		int x1, x2, y1, y2;
		Eigen::Vector3d& temp2 = interp_u_star[i];
		
		//Update the u-velocity
		x1 = indices.first;
		x2 = x1 + 1;
		if(y > indices.second + 0.5){
			y1 = indices.second;
			y2 = y1 + 1;
		}
		else{
			y2 = indices.second;
			y1 = y2 - 1;
		}
		temp2[0] = 1/((x2 - x)*(y2-y))*((MACGrid_->get_u(x1,y1))*(x2 - x)*(y2-y) 
							+ (MACGrid_->get_u(x2,y1))*(x - x)*(y2-y) 
							+ (MACGrid_->get_u(x1,y2))*(x2 - x)*(y-y1) 
							+ (MACGrid_->get_u(x2,y2))*(x - x1)*(y-y1));
							
		//Update the v-velocity
		y1 = indices.second;
		y2 = y1 + 1;
		if(x > indices.first + 0.5){
			x1 = indices.first;
			x2 = x1 + 1;
		}
		else{
			x2 = indices.first;
			x1 = x2 - 1;
		}
		temp2[1] = 1/((x2 - x)*(y2-y))*((MACGrid_->get_u(x1,y1))*(x2 - x)*(y2-y) 
							+ (MACGrid_->get_u(x2,y1))*(x - x)*(y2-y) 
							+ (MACGrid_->get_u(x1,y2))*(x2 - x)*(y-y1) 
							+ (MACGrid_->get_u(x2,y2))*(x - x1)*(y-y1));	
	}
	
	//Applying the pressures we can compute u^n+1 and make the interpolation interp(u^n+1, x_p)
	do_pressures();
	std::vector<Eigen::Vector3d> interp_u_n1(num_particles_);
	for(int i = 0; i < num_particles_; ++i){
		Eigen::Vector3d temp1 = initial_positions[i];
		double x = temp1[0];
		double y = temp1[1];
		Mac2d::Pair_t indices = MACGrid_->index_from_coord(x,y);
		int x1, x2, y1, y2;
		Eigen::Vector3d& temp2 = interp_u_n1[i];
		
		//Update the u-velocity
		x1 = indices.first;
		x2 = x1 + 1;
		if(y > indices.second + 0.5){
			y1 = indices.second;
			y2 = y1 + 1;
		}
		else{
			y2 = indices.second;
			y1 = y2 - 1;
		}
		temp2[0] = 1/((x2 - x)*(y2-y))*((MACGrid_->get_u(x1,y1))*(x2 - x)*(y2-y) 
							+ (MACGrid_->get_u(x2,y1))*(x - x)*(y2-y) 
							+ (MACGrid_->get_u(x1,y2))*(x2 - x)*(y-y1) 
							+ (MACGrid_->get_u(x2,y2))*(x - x1)*(y-y1));
							
		//Update the v-velocity
		y1 = indices.second;
		y2 = y1 + 1;
		if(x > indices.first + 0.5){
			x1 = indices.first;
			x2 = x1 + 1;
		}
		else{
			x2 = indices.first;
			x1 = x2 - 1;
		}
		temp2[1] = 1/((x2 - x)*(y2-y))*((MACGrid_->get_u(x1,y1))*(x2 - x)*(y2-y) 
							+ (MACGrid_->get_u(x2,y1))*(x - x)*(y2-y) 
							+ (MACGrid_->get_u(x1,y2))*(x2 - x)*(y-y1) 
							+ (MACGrid_->get_u(x2,y2))*(x - x1)*(y-y1));	
	}
	
	//Compute the updated velocity for every particle
	for(int i = 0; i < num_particles_; ++i){
		Eigen::Vector3d temp3 = initial_velocities[i] + interp_u_n1[i] + interp_u_star[1]*(alpha - 1);
		(particles_ + i)->set_velocity(temp3);
	}
}

void FLIP::advance_particles() {
	// TODO: update particle positions 
	//  - Use RK2 interpolator	
}

