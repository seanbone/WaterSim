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
	
/**********************************
 * Particle-to-Grid Transfer
***********************************/

	// Set all grid velocities to zero
	MACGrid_->set_velocities_to_zero();
	MACGrid_->set_weights_to_zero();
	
	// Positions and velocities of a single particle
	Eigen::Vector3d pos;
	Eigen::Vector3d vel;
	
	// Coordinates of a cell
	Eigen::Vector3d cell_coord;
	
	// Sizes of the edges of a cell (in meters)
	double cell_sizex = MACGrid_->get_cell_sizex();
	double cell_sizey = MACGrid_->get_cell_sizey();
	
	// Threshold h and h scaled so that it is equal to the distance expressed in number of cells
	double h = cell_sizex;
	int h_scaledx = h/cell_sizex;
	int h_scaledy = h/cell_sizey;
	
	// Iterate over all particles and add weighted particles velocities
	// to grid points within a threshold h (in this case equal to the 
	// length of an edge of a cell)
	for( unsigned int n = 0; n < num_particles_; ++n ){
		pos = (*(particles_ + n)).get_position();
		vel = (*(particles_ + n)).get_velocity();
		
		Pair_t tmp = MACGrid_->index_from_coord(pos(0), pos(1));
		cell_coord << tmp.first, tmp.second, 0;
		
		for( int j = cell_coord(1) - h_scaledy; j < cell_coord(1) + h_scaledy; ++j ){
			for( int i = cell_coord(0) - h_scaledx; i < cell_coord(0) + h_scaledx; ++i ){
				// Coordinates of the points on the grid edges
				Eigen::Vector3d grid_coord;
				
				// Left edge
				grid_coord << cell_coord(0)*cell_sizex, (cell_coord(1) + 0.5)*cell_sizey, 0;
				accumulate_u(pos, vel, grid_coord, h, i, j);
				
				// Lower edge
				grid_coord(1) -= 0.5 * cell_sizey;
				grid_coord(0) += 0.5 * cell_sizex;
				accumulate_v(pos, vel, grid_coord, h, i, j);
			}
		}
	}
	
	
}

bool FLIP::check_threshold( const Eigen::Vector3d& particle_coord, 
					  const Eigen::Vector3d& grid_coord, 
					  const double h )
{
	if ( (particle_coord - grid_coord).norm() < h )
		return true;
	
	return false;
}

double FLIP::compute_weight( const Eigen::Vector3d& particle_coord, 
							 const Eigen::Vector3d& grid_coord, 
							 const double h )
{
	double r = (particle_coord - grid_coord).norm();
	return ( (315/(64 * M_PI * pow(h, 9))) * (pow(h, 2) - pow(r, 2)) ) * r;
}

// Accumulate velocities and weights for u					
void FLIP::accumulate_u( const Eigen::Vector3d& pos,
						 const Eigen::Vector3d& vel,
						 const Eigen::Vector3d& grid_coord,
						 const double h,
						 const int i,
						 const int j )
{
	if ( check_threshold(pos, grid_coord, h) ){
		double u_pred = MACGrid_->get_u(i, j);
		double W_u = compute_weight(pos, grid_coord, h);
		double u_curr = u_pred + (W_u * vel(0));
		
		// Accumulate velocities
		MACGrid_->set_u(i, j, u_curr);
		
		// Accumulate weights
		double W_u_pred = MACGrid_->get_weights_u(i, j);
		double W_u_curr = W_u_pred + W_u;
		MACGrid_->set_weights_u(i, j, W_u_curr);
	}
}

// Accumulate velocities and weights for v
void FLIP::accumulate_v( const Eigen::Vector3d& pos,
						 const Eigen::Vector3d& vel,
						 const Eigen::Vector3d& grid_coord,
						 const double h,
						 const int i,
						 const int j )
{
	if ( check_threshold(pos, grid_coord, h) ){
		double v_pred = MACGrid_->get_v(i, j);
		double W_v = compute_weight(pos, grid_coord, h);
		double v_curr = v_pred + (W_v * vel(1));
		
		// Accumulate velocities
		MACGrid_->set_v(i, j, v_curr);
		
		// Accumulate weights
		double W_v_pred = MACGrid_->get_weights_v(i, j);
		double W_v_curr = W_v_pred + W_v;
		MACGrid_->set_weights_v(i, j, W_v_curr);
	}
}

void FLIP::normalize_accumulated_u(){
	for( int j = 0; j < M_+1; ++j ){
		for( int i = 0; i < N_+1; ++i ){
			double W_u = MACGrid_->get_weights_u(i, j);
			if ( W_u != 0 ){
				double u_pred = MACGrid_->get_u(i, j);
				double u_curr = u_pred/W_u;
				MACGrid_->set_u(i, j, u_curr);
			}
		}	
	}	
}

void FLIP::normalize_accumulated_v(){
	for( int j = 0; j < M_+1; ++j ){
		for( int i = 0; i < N_+1; ++i ){
			double W_v = MACGrid_->get_weights_v(i, j);
			if ( W_v != 0 ){
				double v_pred = MACGrid_->get_v(i, j);
				double v_curr = v_pred/W_v;
				MACGrid_->set_v(i, j, v_curr);
			}
		}	
	}	
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

