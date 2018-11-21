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
	 * 1a. Copy velocity field to intermediate velocity field u^*
	 * 2. Apply external forces (fwd euler on field)
	 * 3. Compute & apply pressure gradients
	 * 4. Update particle velocities
	 * 5. Update particle positions
	 */

	// TODO: subsample time interval to satisfy CFL condition
	double cur_time = time;

	//std::cout << MACGrid_->get_v(10, 40) << "\n\n";
	// 1.
	compute_velocity_field();

	// 1a.
	MACGrid_->set_uv_star();

	//std::cout << MACGrid_->get_v(10, 40) << "\n\n";
	//std::cout << particles_->get_velocity() << std::endl;
	//std::cout << ".....\n";
	// 2.
	apply_forces(dt);

	// 3.
	//~ do_pressures(dt);

	//~ std::cout << "Pressure for cell (5, 1): " << MACGrid_->get_pressure(5, 1) << std::endl;

	// 4.
	grid_to_particle();

	// 5.
	advance_particles(dt, step);
}

/*** COMPUTE VELOCITY FIELD ***/
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
	
	// Reset all fluid flags
	MACGrid_->reset_fluid();
	
	// Iterate over all particles and add weighted particles velocities
	// to grid points within a threshold h (in this case equal to the 
	// length of an edge of a cell)
	for( unsigned int n = 0; n < num_particles_; ++n ){
		pos = (*(particles_ + n)).get_position();
		vel = (*(particles_ + n)).get_velocity();
		
		Mac2d::Pair_t tmp = MACGrid_->index_from_coord(pos(0), pos(1));
		cell_coord << tmp.first, tmp.second, 0;
		
		// Set the cell of the current particle to a fluid-cell
		if ( !(MACGrid_->is_fluid(cell_coord(0), cell_coord(1))) and !(MACGrid_->is_solid(cell_coord(0), cell_coord(1))) ){
			MACGrid_->set_fluid(cell_coord(0), cell_coord(1));
		}
		
		// Coordinates of the points on the grid edges
		Eigen::Vector3d grid_coord;
		grid_coord << 0, 0, 0;
		for( int j = cell_coord(1) - h_scaledy; j < cell_coord(1) + h_scaledy + 1; ++j ){
			for( int i = cell_coord(0) - h_scaledx; i < cell_coord(0) + h_scaledx + 1; ++i ){
				if ( ( i >= 0 and j >= 0 ) ){
					if ( ( i <= MACGrid_->get_num_cells_x() and j < MACGrid_->get_num_cells_y() ) ){
					
						// Left edge
						grid_coord(0) = i * cell_sizex;
						grid_coord(1) = (j + 0.5) * cell_sizey;
						accumulate_u(pos, vel, grid_coord, h, i, j);
					}
					
					if ( ( i < MACGrid_->get_num_cells_x() and j <= MACGrid_->get_num_cells_y() ) ){
						// Lower edge
						grid_coord(0) += 0.5 * cell_sizex;
						grid_coord(1) -= 0.5 * cell_sizey;
						accumulate_v(pos, vel, grid_coord, h, i, j);
					}
				}
			}
		}
	}
	
	// Normalize grid-velocities
	normalize_accumulated_u();
	normalize_accumulated_v();
}

bool FLIP::check_threshold( const Eigen::Vector3d& particle_coord, 
					  const Eigen::Vector3d& grid_coord, 
					  const double h )
{
	if ( (particle_coord - grid_coord).norm() < h ) {
		return true;
	}
	
	return false;
}

double FLIP::compute_weight( const Eigen::Vector3d& particle_coord, 
							 const Eigen::Vector3d& grid_coord, 
							 const double h )
{
	double r = (particle_coord - grid_coord).norm();
	return ( (315/(64 * M_PI * pow(h, 9))) * (pow(h, 2) - pow(r, 2)) );
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
		double u_prev = MACGrid_->get_u(i, j);
		double W_u = compute_weight(pos, grid_coord, h);
		double u_curr = u_prev + (W_u * vel(0));
		
		// Accumulate velocities
		MACGrid_->set_u(i, j, u_curr);
		
		// Accumulate weights
		double W_u_prev = MACGrid_->get_weights_u(i, j);
		double W_u_curr = W_u_prev + W_u;
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
		double v_prev = MACGrid_->get_v(i, j);
		double W_v = compute_weight(pos, grid_coord, h);
		double v_curr = v_prev + (W_v * vel(1));
		
		// Accumulate velocities
		MACGrid_->set_v(i, j, v_curr);
		
		// Accumulate weights
		double W_v_prev = MACGrid_->get_weights_v(i, j);
		double W_v_curr = W_v_prev + W_v;
		MACGrid_->set_weights_v(i, j, W_v_curr);
	}
}

void FLIP::normalize_accumulated_u(){
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	for( unsigned j = 0; j < M; ++j ){
		for( unsigned i = 0; i < N+1; ++i ){
			double W_u = MACGrid_->get_weights_u(i, j);
			if ( W_u != 0 ){
				double u_prev = MACGrid_->get_u(i, j);
				double u_curr = u_prev/W_u;
				MACGrid_->set_u(i, j, u_curr);
			}
		}	
	}	
}

void FLIP::normalize_accumulated_v(){
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	for( unsigned j = 0; j < M+1; ++j ){
		for( unsigned i = 0; i < N; ++i ){
			double W_v = MACGrid_->get_weights_v(i, j);
			if ( W_v != 0 ){
				double v_prev = MACGrid_->get_v(i, j);
				double v_curr = v_prev/W_v;
				MACGrid_->set_v(i, j, v_curr);
			}
		}	
	}	
}

/*** APPLY EXTERNAL FORCES ***/
void FLIP::apply_forces(const double dt) {
	// Compute&apply external forces (gravity, vorticity confinement, ...) 
	// Apply them to the velocity field via forward euler
	// Only worry about gravity for now
	auto& g = MACGrid_;
	const unsigned N = g->get_num_cells_x();
	const unsigned M = g->get_num_cells_y();
	
	// Iterate over cells & update: dv = dt*g
	for (unsigned j = 0; j <= M; j++) {
		for (unsigned i = 0; i < N; i++) {
			g->set_v(i, j, g->get_v(i,j) + dt*gravity_mag_);
		}
	}
}

/*** PRESSURE SOLVING ***/
void FLIP::do_pressures(const double dt) {
	// 3. Compute & apply pressure gradients to field
	
	//   3a. Compute A matrix
	compute_pressure_matrix();

	//   3b. Compute rhs d
	compute_pressure_rhs(dt);

	std::cout << "nnz: " << A_.nonZeros() << std::endl;

	//   3c. Solve for p: Ap = d (MICCG(0))
	// *TODO: use MICCG(0) solver
	using namespace Eigen;
	using solver_t = ConjugateGradient<SparseMatrix<double>, Lower|Upper>;
	solver_t solver;
	solver.compute(A_);
	VectorXd p = solver.solve(d_);

	// Copy pressures to MAC grid
	//MACGrid_->set_pressure(p);
	
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	unsigned cellidx = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (MACGrid_->is_fluid(i, j)) {
				MACGrid_->set_pressure(i, j, p(cellidx));
			} else {
				MACGrid_->set_pressure(i, j, 0);
			} // is_fluid(i,j)
		}
	}

	//  3d. Apply pressure gradients to velocity field
	//     -> see SIGGRAPH §4
	//~ std::cout << "v before" << std::endl;
	//~ std::cout << particles_->get_velocity() << std::endl;
	apply_pressure_gradients(dt);
	//~ std::cout << "v after" << std::endl;
	//~ std::cout << particles_->get_velocity() << std::endl;

	// Note: boundary conditions are handles here by setting the pressures
	//    such that no particle exits the system.
}

void FLIP::compute_pressure_matrix() {
	// Compute matrix for pressure solve and store in A_
	//  See eq. (4.19) and (4.24) in SIGGRAPH notes
	
	std::vector< Mac2d::Triplet_t > triplets;
	
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();


	unsigned num_fluid_cells = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (MACGrid_->is_fluid(i, j))
				num_fluid_cells++;
		}
	}

	unsigned cellidx = 0;
	// *TODO: reduce matrix size from (N*M)x(N*M) to (#fluidcells)x(#fluidcells)
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (MACGrid_->is_fluid(i, j)) {
				// Copy diagonal entry
				auto diag_e = MACGrid_->get_a_diag()[i + j*nx];
				triplets.push_back(Mac2d::Triplet_t(cellidx, cellidx, diag_e.value()));
				
				// Compute off-diagonal entries
				// x-adjacent cells
				if (cellidx+1 < num_fluid_cells && MACGrid_->is_fluid(i+1, j)) {
						triplets.push_back(Mac2d::Triplet_t(cellidx, cellidx+1, -1));
						// Use symmetry to avoid computing (i-1,j) separately
						triplets.push_back(Mac2d::Triplet_t(cellidx+1, cellidx, -1));
				}
				// y-adjacent cells
				if (cellidx+ny < num_fluid_cells && MACGrid_->is_fluid(i, j+1)) {
						triplets.push_back(Mac2d::Triplet_t(cellidx, cellidx + ny, -1));
						// Use symmetry to avoid computing (i,j-1) separately
						triplets.push_back(Mac2d::Triplet_t(cellidx + ny, cellidx, -1));
				}

				cellidx++;
			} // if is_fluid(i,j)
		}
	} // outer for

	A_.setZero();
	A_.resize(num_fluid_cells, num_fluid_cells);
	A_.setFromTriplets(triplets.begin(), triplets.end());
	//~ std::cout << "nnz: " << A_.nonZeros() << std::endl;
}

void FLIP::compute_pressure_rhs(const double dt) {
	// Compute right-hand side of the pressure equations and store in d_
	//  See eq. (4.19) and (4.24) in SIGGRAPH notes
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();

	// Alias for MAC grid
	auto& g = MACGrid_;

	unsigned num_fluid_cells = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (MACGrid_->is_fluid(i, j))
				num_fluid_cells++;
		}
	}

	d_.setZero();
	d_.resize(num_fluid_cells);
	
	// *TODO: reduce vector size from (N*M) to (#fluidcells)
	unsigned cellidx = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (g->is_fluid(i,j)) {
				// get_u(i,j) = u_{ (i-1/2, j) }
				double d_ij = -(g->get_u(i+1,j) - g->get_u(i,j));
				d_ij -= g->get_v(i,j+1) - g->get_v(i,j);
	
				// Note: u_{solid} = 0
				// (i+1, j)
				if ((i < (nx-1) && g->is_solid(i+1,j)) || i == nx-1) {
					d_ij += g->get_u(i+1,j);
				}
				// (i-1, j)
				if ((i > 0 && g->is_solid(i-1,j)) || i == 0) {
					d_ij += g->get_u(i-1,j);
				}
	
				// (i, j+1)
				if ((j < (ny-1) && g->is_solid(i,j+1)) || j == ny-1) {
					d_ij += g->get_v(i,j+1);
				}
				// (i, j-1)
				if ((j > 0 && g->is_solid(i,j-1)) || j == 0) {
					d_ij += g->get_v(i,j-1);
				}
	
				// *TODO: cell x and y size should always be the same -> enforce via sim params?
				d_(cellidx) = fluid_density_ * g->get_cell_sizex() * d_ij / dt;

				cellidx++;
			} // if is_fluid(i,j)
		}
	}
}

void FLIP::apply_pressure_gradients(const double dt) {
	// Apply pressure gradients to velocity field
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	// MACGrid alias
	auto& g = MACGrid_;
	// *TODO: correct cell dimensions
	double dx = g->get_cell_sizex();
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (i != 0) {
				// get_u(i,j) = u_{ (i-1/2, j) }
				// See SIGGRAPH eq. (4.4)
				double du = (g->get_pressure(i,j) - g->get_pressure(i-1,j));
				du *= (dt/(dx*fluid_density_));
				g->set_u(i,j, g->get_u(i,j) - du);
			}
			if (j != 0) {
				// get_v(i,j) = v_{ (i, j-1/2) }
				// See SIGGRAPH eq. (4.5)
				double dv = (g->get_pressure(i,j) - g->get_pressure(i,j-1));
				dv *= (dt/(dx*fluid_density_));
				g->set_v(i,j, g->get_v(i,j) - dv);
			}
		}
	}
}



/*** UPDATE PARTICLE VELOCITIES & MOVE PARTICLES ***/
void FLIP::grid_to_particle(){
	// FLIP grid to particle transfer
	//  -> See slides Fluids II, FLIP_explained.pdf
	double alpha = 0.0;
	
	for(int i = 0; i < num_particles_; ++i){
		//Store the initial positions and velocities of the particles
		Eigen::Vector3d initial_position = (particles_+i)->get_position();
		Eigen::Vector3d initial_velocity = (particles_+i)->get_velocity();
		
		//Initialization of the variables
		Eigen::Vector3d interp_u_star;
		interp_u_star.setZero();
		Eigen::Vector3d interp_u_n1;
		interp_u_n1.setZero();
		Eigen::Vector3d u_update;
		u_update.setZero();

		double x = initial_position[0];
		double y = initial_position[1];
		double sx = MACGrid_->get_cell_sizex();
		double sy = MACGrid_->get_cell_sizex();
		Mac2d::Pair_t indices = MACGrid_->index_from_coord(x,y);
		double x1, x2, y1, y2;
		int ix1, ix2, iy1, iy2;
		double u11, u12, u21, u22;
		double v11, v12, v21, v22;
		
		//With u* and v* we can make the interpolation interp(u*, x_p),
		//with the new u and v we can make the interpolation interp(u_n1, x_p)
		
		//Update the u-velocity (bilinear interpolation)
		ix1 = indices.first;
		ix2 = ix1 + 1;
		if(y > (indices.second + 0.5) * sy){
			iy1 = indices.second;
			iy2 = iy1 + 1;
		}
		else{
			iy2 = indices.second;
			iy1 = iy2 - 1;
		}
		x1 = ix1 * sx;
		x2 = ix2 * sx;
		y1 = iy1 * sy;
		y2 = iy2 * sy;

		u11 = MACGrid_->get_u_star(ix1,iy1);
		u12 = MACGrid_->get_u_star(ix1,iy2);
		u21 = MACGrid_->get_u_star(ix2,iy1);
		u22 = MACGrid_->get_u_star(ix2,iy2);
		interp_u_star[0] = 1/((x2 - x1)*(y2-y1))*(u11*(x2 - x)*(y2-y) 
							+ u21*(x - x1)*(y2-y) 
							+ u12*(x2 - x)*(y-y1) 
							+ u22*(x - x1)*(y-y1));
		
		u11 = MACGrid_->get_u(ix1,iy1);
		u12 = MACGrid_->get_u(ix1,iy2);
		u21 = MACGrid_->get_u(ix2,iy1);
		u22 = MACGrid_->get_u(ix2,iy2);
		interp_u_n1[0] = 1/((x2 - x1)*(y2-y1))*(u11*(x2 - x)*(y2-y) 
							+ u21*(x - x1)*(y2-y) 
							+ u12*(x2 - x)*(y-y1) 
							+ u22*(x - x1)*(y-y1));
		
		//Update the v-velocity (bilinear interpolation)
		iy1 = indices.second;
		iy2 = iy1 + 1;
		if(x > (indices.first + 0.5) * sx){
			ix1 = indices.first;
			ix2 = ix1 + 1;
		}
		else{
			ix2 = indices.first;
			ix1 = ix2 - 1;
		}
		x1 = ix1 * sx;
		x2 = ix2 * sx;
		y1 = iy1 * sy;
		y2 = iy2 * sy;
		
		v11 = MACGrid_->get_v_star(ix1,iy1);
		v12 = MACGrid_->get_v_star(ix1,iy2);
		v21 = MACGrid_->get_v_star(ix2,iy1);
		v22 = MACGrid_->get_v_star(ix2,iy2);
		interp_u_star[1] = 1/((x2 - x1)*(y2-y1))*(v11*(x2 - x)*(y2-y) 
							+ v21*(x - x1)*(y2-y) 
							+ v12*(x2 - x)*(y-y1) 
							+ v22*(x - x1)*(y-y1));

		v11 = MACGrid_->get_v(ix1,iy1);
		v12 = MACGrid_->get_v(ix1,iy2);
		v21 = MACGrid_->get_v(ix2,iy1);
		v22 = MACGrid_->get_v(ix2,iy2);
		interp_u_n1[1] = 1/((x2 - x1)*(y2-y1))*(v11*(x2 - x)*(y2-y) 
							+ v21*(x - x1)*(y2-y) 
							+ v12*(x2 - x)*(y-y1) 
							+ v22*(x - x1)*(y-y1));
		
		//Update the final velocity of the particles
		u_update = initial_velocity*(1 - alpha) + interp_u_n1 + interp_u_star*(alpha - 1);
		(particles_ + i)->set_velocity(u_update);
	}
}


void FLIP::advance_particles(const double dt, const unsigned step) {
	//Se una particles esce dal sistema o entra in un solido, rispingerla dentro.
	// TODO: update particle positions 
	//  - Use RK2 interpolator	
	for(unsigned n = 0; n < num_particles_; ++n){
		Eigen::Vector3d pos_prev = (particles_ + n)->get_prev_position();
		Eigen::Vector3d pos_curr = (particles_ + n)->get_position();
		Eigen::Vector3d vel = (particles_ + n)->get_velocity();
		
		Eigen::Vector3d pos_next;
		if (step == 0) {
			// Euler
			pos_next = pos_curr + dt*vel;
		} else {
			// Leapfrog
			pos_next = pos_prev + 2*dt*vel;
		}
		if (n == 0) {
			std::cout << "pos_prev:\n" << pos_prev;
			std::cout << "\npos_curr:\n" << pos_curr;
			std::cout << "\npos_next:\n" << pos_next;
			std::cout << "\nvel:\n" << vel;
			std::cout << "\n------------\n";
		}
		
		double x = pos_next(0);
		double y = pos_next(1);
		//Check if the particle exits the grid
		double size_x = MACGrid_->get_cell_sizex() * MACGrid_->get_num_cells_x();
		double size_y = MACGrid_->get_cell_sizey() * MACGrid_->get_num_cells_y();
		if(x < 0) {
			pos_next(0) = 0;
			vel(0) = 0;
		}
		if(x > size_x) {
			pos_next(0) = size_x;
			vel(0) = 0;
		}
		if(y < 0) {
			pos_next(1) = 0;
			vel(1) = 0;
		}
		if(y > size_y) {
			pos_next(1) = size_y;
			vel(1) = 0;
		}
		//(particles_ + n)->set_velocity(vel);
		

		//Check if the particle enters in a solid
		Mac2d::Pair_t prev_indices = MACGrid_->index_from_coord(pos_curr(0),pos_curr(1));
		Mac2d::Pair_t new_indices = MACGrid_->index_from_coord(pos_next(0), pos_next(1));
		double sx = MACGrid_->get_cell_sizex();
		double sy = MACGrid_->get_cell_sizey();
		//TODO: correctly shift particles & velocities
		if (MACGrid_->is_solid(new_indices.first, new_indices.second)) {
			if (prev_indices.first  > new_indices.first)
				pos_next(0) = (prev_indices.first - 0.25) * sx;
			else if (prev_indices.first < new_indices.first)
				pos_next(0) = (prev_indices.first + 0.25) * sx;

			if (prev_indices.second > new_indices.second)
				pos_next(1) = (prev_indices.second - 0.25) * sy;
			else if (prev_indices.second < new_indices.second)
				pos_next(1) = (prev_indices.second + 0.25) * sy;
		} else {
			//(particles_ + n)->set_velocity(vel(0), 0, vel(2));
		}
		(particles_ + n)->set_position(pos_next);
	}
}

