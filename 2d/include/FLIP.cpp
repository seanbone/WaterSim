#include "FLIP.h"

FLIP::FLIP(Particle* particles, const unsigned num_particles, Mac2d* MACGrid,
		   const double density, const double gravity, const double alpha) 
	: particles_(particles), num_particles_(num_particles), MACGrid_(MACGrid),
	  fluid_density_(density), gravity_mag_(gravity), alpha_(alpha) {
	
}

/**
 * Advance FLIP simulation by one frame
 */
void FLIP::step_FLIP(const double dt, const unsigned long step) {
	/** One FLIP step:
	 * 1. Compute velocity field (particle-to-grid transfer)
	 *    - Particle-to-grid transfer
	 *    - Classify cells (fluid/air)
	 *    - Extrapolate velocity field into air region
	 * 1a. Copy velocity field to intermediate velocity field u^*
	 * 2. Apply external forces (fwd euler on field)
	 * 3. Enforce boundary conditions for grid boundaries
	 * 4. Compute & apply pressure gradients
	 * 5. Update particle velocities
	 * 6. Update particle positions
	 */
	
	// 1.
	compute_velocity_field();

	// 1a.
	MACGrid_->set_uv_star();

	// 2.
	apply_forces(dt);
	
	// Uncomment to add "meteor splash"
	//if ( step >= 5 and step <= 20 ){
	//	explode(dt, step, 300);
	//}
	
	// 3.
	apply_boundary_conditions();

	// 4.
	do_pressures(dt);

	// 5.
	grid_to_particle();
	
	// 6. subsample time interval to satisfy CFL condition
	double dt_new = compute_timestep(dt);
	double num_substeps = std::ceil(dt/dt_new);
	for( int s = 0; s < num_substeps ; ++s ){
		
		// 7.
		advance_particles(dt/num_substeps, step);
	}
}

double FLIP::compute_timestep( const double dt ){
	double dt_new;
	
	Eigen::Vector3d vel;
	double u_max = 0;
	double v_max = 0;
	for( unsigned int n = 0; n < num_particles_; ++n ){
		vel = (particles_ + n)->get_velocity();
		if ( std::abs(vel(0)) > std::abs(u_max) ){
			u_max = vel(0);
		}
		if ( std::abs(vel(1)) > std::abs(v_max) ){
			v_max = vel(1);
		}
	}
	
	if ( u_max == 0 ){
		dt_new = dt;
	} else {
		dt_new = std::abs(MACGrid_->get_cell_sizex()/u_max);
		if ( dt_new > dt){
			dt_new = dt;
		}
	}
	
	if ( v_max == 0 ){
		dt_new = dt;
	} else {
		double tmp = std::abs(MACGrid_->get_cell_sizey()/v_max);
		if ( tmp < dt_new){
			dt_new = tmp;
		}
		if ( dt_new > dt ){
			dt_new = dt;
		}
	}
	
	return dt_new;
}

/*** COMPUTE VELOCITY FIELD ***/
void FLIP::compute_velocity_field() {
	// 1. Compute the velocity field (velocities on grid)
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
	double h = 2*cell_sizex;
	int h_scaledx = std::ceil(h/cell_sizex);
	int h_scaledy = std::ceil(h/cell_sizey);
	
	// Lists of flags for visited grid-velocities: 1 -> visited
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned M = MACGrid_->get_num_cells_y();
	bool* visited_u = new bool[M*(N+1)];
	bool* visited_v = new bool[N*(M+1)];
	std::fill(visited_u, visited_u + M*(N+1), 0);
	std::fill(visited_v, visited_v + N*(M+1), 0);
	
	// Reset all fluid flags
	MACGrid_->reset_fluid();
	
	// Iterate over all particles and add weighted particles velocities
	// to grid points within a threshold h (in this case equal to the 
	// length of an edge of a cell)
	for( unsigned int n = 0; n < num_particles_; ++n ){
		pos = (particles_ + n)->get_position();
		vel = (particles_ + n)->get_velocity();

		Mac2d::Pair_t tmp = MACGrid_->index_from_coord(pos(0), pos(1));
		cell_coord << tmp.first, tmp.second, 0;
		
		// Set the cell of the current particle to a fluid-cell
		if ( !(MACGrid_->is_fluid(cell_coord(0), cell_coord(1))) and !(MACGrid_->is_solid(cell_coord(0), cell_coord(1))) ){
			MACGrid_->set_fluid(cell_coord(0), cell_coord(1));
		}
		
		// Coordinates of the points on the grid edges
		Eigen::Vector3d grid_coord;
		grid_coord << 0, 0, 0;
		int nx = MACGrid_->get_num_cells_x();
		int ny = MACGrid_->get_num_cells_y();
		for( int j = cell_coord(1) - h_scaledy; j <= cell_coord(1) + h_scaledy + 1; ++j ){
			for( int i = cell_coord(0) - h_scaledx; i <= cell_coord(0) + h_scaledx + 1; ++i ){
				if ( ( i >= 0 and j >= 0 ) ){
					if ( ( i <= nx and j < ny ) ){
					
						// Left edge
						grid_coord(0) = (i - 0.5) * cell_sizex;
						grid_coord(1) = j * cell_sizey;
						accumulate_u(pos, vel, grid_coord, h, i, j);
					}
					
					if ( ( i < nx and j <= ny ) ){
						
						// Lower edge
						grid_coord(0) = i * cell_sizex;
						grid_coord(1) = (j - 0.5) * cell_sizey;
						accumulate_v(pos, vel, grid_coord, h, i, j);
					}
				}
			}
		}
	}
	
	// Normalize grid-velocities
	normalize_accumulated_u( visited_u );
	normalize_accumulated_v( visited_v );
	
	// Extrapolate velocities
	extrapolate_u( visited_u );
	extrapolate_v( visited_v );
	
	delete[] visited_u;
	delete[] visited_v;
}

bool FLIP::check_threshold( const Eigen::Vector3d& particle_coord, 
					  const Eigen::Vector3d& grid_coord, 
					  const double h )
{
	if ( (particle_coord - grid_coord).norm() <= h ) {
		return true;
	}
	
	return false;
}

double FLIP::compute_weight( const Eigen::Vector3d& particle_coord, 
							 const Eigen::Vector3d& grid_coord, 
							 const double h )
{
	double r = (particle_coord - grid_coord).norm();
	double diff = std::pow(h, 2) - std::pow(r, 2);
	return ( (315/(64 * M_PI * std::pow(h, 9))) * std::pow(diff, 3) );
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

void FLIP::normalize_accumulated_u( bool* const visited_u ){
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	for( unsigned j = 0; j < M; ++j ){
		for( unsigned i = 0; i < N+1; ++i ){
			double W_u = MACGrid_->get_weights_u(i, j);
			if ( W_u != 0 ){
				double u_prev = MACGrid_->get_u(i, j);
				double u_curr = u_prev/W_u;
				MACGrid_->set_u(i, j, u_curr);
				*(visited_u + (N+1)*j + i) = 1;
			}
		}	
	}	
}

void FLIP::normalize_accumulated_v( bool* const visited_v ){
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	for( unsigned j = 0; j < M+1; ++j ){
		for( unsigned i = 0; i < N; ++i ){
			double W_v = MACGrid_->get_weights_v(i, j);
			if ( W_v != 0 ){
				double v_prev = MACGrid_->get_v(i, j);
				double v_curr = v_prev/W_v;
				MACGrid_->set_v(i, j, v_curr);
				*(visited_v + N*j + i) = 1;
			}
		}	
	}	
}

void FLIP::extrapolate_u( const bool* const visited_u ){
	// Do the cases for upper and right bound
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned* counter = new unsigned[M*(N+1)];
	std::fill(counter, counter + M*(N+1), 0);
	for( unsigned j = 0; j < M; ++j ){
		for( unsigned i = 0; i < N+1; ++i ){
			if ( *(visited_u + (N+1)*j + i) ){
				if ( i != 0 and !(*(visited_u + (N+1)*j + (i-1))) ){
					double tmp = MACGrid_->get_u(i-1, j) * *(counter + (N+1)*j + (i-1));
					*(counter + (N+1)*j + (i-1)) += 1;
					MACGrid_->set_u(i-1, j, (tmp + MACGrid_->get_u(i, j))/(*(counter + (N+1)*j + (i-1))));
				}
				if ( j != 0 and !(*(visited_u + (N+1)*(j-1) + i)) ){
					double tmp = MACGrid_->get_u(i, j-1) * *(counter + (N+1)*(j-1) + i);
					*(counter + (N+1)*(j-1) + i) += 1;
					MACGrid_->set_u(i, j-1, (tmp + MACGrid_->get_u(i, j))/(*(counter + (N+1)*(j-1) + i)));
				}
				if ( i != N and !(*(visited_u + (N+1)*j + (i+1))) ){
					double tmp = MACGrid_->get_u(i+1, j) * *(counter + (N+1)*j + (i+1));
					*(counter + (N+1)*j + (i+1)) += 1;
					MACGrid_->set_u(i+1, j, (tmp + MACGrid_->get_u(i, j))/(*(counter + (N+1)*j + (i+1))));
				}
				if ( j != M-1 and !(*(visited_u + (N+1)*(j+1) + i)) ){
					double tmp = MACGrid_->get_u(i, j+1) * *(counter + (N+1)*(j+1) + i);
					*(counter + (N+1)*(j+1) + i) += 1;
					MACGrid_->set_u(i, j+1, (tmp + MACGrid_->get_u(i, j))/(*(counter + (N+1)*(j+1) + i)));
				}
			}
		}
	}
	
	delete[] counter;
}

void FLIP::extrapolate_v( const bool* const visited_v ){
	// Do the cases for upper and right bound
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned* counter = new unsigned[N*(M+1)];
	std::fill(counter, counter + N*(M+1), 0);
	for( unsigned j = 0; j < M+1; ++j ){
		for( unsigned i = 0; i < N; ++i ){
			if ( *(visited_v + N*j + i) ){
				if ( i != 0 and !(*(visited_v + N*j + (i-1))) ){
					double tmp = MACGrid_->get_v(i-1, j) * *(counter + N*j + (i-1));
					*(counter + N*j + (i-1)) += 1;
					MACGrid_->set_v(i-1, j, (tmp + MACGrid_->get_v(i, j))/(*(counter + N*j + (i-1))));
				}
				if ( j != 0 and !(*(visited_v + N*(j-1) + i)) ){
					double tmp = MACGrid_->get_v(i, j-1) * *(counter + N*(j-1) + i);
					*(counter + N*(j-1) + i) += 1;
					MACGrid_->set_v(i, j-1, (tmp + MACGrid_->get_v(i, j))/(*(counter + N*(j-1) + i)));
				}
				if ( i != N-1 and !(*(visited_v + N*j + (i+1))) ){
					double tmp = MACGrid_->get_v(i+1, j) * *(counter + N*j + (i+1));
					*(counter + N*j + (i+1)) += 1;
					MACGrid_->set_v(i+1, j, (tmp + MACGrid_->get_v(i, j))/(*(counter + N*j + (i+1))));
				}
				if ( j != M and !(*(visited_v + N*(j+1) + i)) ){
					double tmp = MACGrid_->get_v(i, j+1) * *(counter + N*(j+1) + i);
					*(counter + N*(j+1) + i) += 1;
					MACGrid_->set_v(i, j+1, (tmp + MACGrid_->get_v(i, j))/(*(counter + N*(j+1) + i)));
				}
			}
		}
	}
	
	delete[] counter;
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
			g->set_v(i, j, g->get_v(i,j) - dt*gravity_mag_);
		}
	}
}


/*** BOUNDARY CONDITIONS ***/
void FLIP::apply_boundary_conditions() {
	// Enforce boundary conditions for grid & solid boundaries
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	// Solid walls
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			bool ij_solid = MACGrid_->is_solid(i,j);
			if (ij_solid || MACGrid_->is_solid(i+1,j))
				MACGrid_->set_u(i+1, j, 0);
			if (ij_solid || MACGrid_->is_solid(i,j+1))
				MACGrid_->set_v(i, j+1, 0);
		}
	}

	// Outer (system) boundaries
	for (unsigned i = 0; i < nx; i++) {
		MACGrid_->set_v(i, 0, 0);
		MACGrid_->set_v(i, ny, 0);
	}
	for (unsigned j = 0; j < ny; j++) {
		MACGrid_->set_u(0, j, 0);
		MACGrid_->set_u(nx, j, 0);
	}
}


/*** PRESSURE SOLVING ***/
void FLIP::do_pressures(const double dt) {
	// Compute & apply pressure gradients to field
	
	// Compute A matrix
	compute_pressure_matrix();

	// Compute rhs d
	compute_pressure_rhs(dt);

	// Solve for p: Ap = d (MICCG(0))
	using namespace Eigen;
	using solver_t = ConjugateGradient<SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double> >;

	//MatrixXd A = A_;
	solver_t solver;
	solver.setMaxIterations(100);
	solver.compute(A_);
	VectorXd p = solver.solve(d_);

	// Copy pressures to MAC grid
	MACGrid_->set_pressure(p);

	// Apply pressure gradients to velocity field
	//     -> see SIGGRAPH §4
	apply_pressure_gradients(dt);
}

void FLIP::compute_pressure_matrix() {
	// Compute matrix for pressure solve and store in A_
	//  See eq. (4.19) and (4.24) in SIGGRAPH notes
	
	std::vector< Mac2d::Triplet_t > triplets;
	
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();

	unsigned cellidx = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++, cellidx++) {
			// Copy diagonal entry
			auto& diag_e = MACGrid_->get_a_diag()[i + j*nx];
			triplets.push_back(diag_e);

			// Compute off-diagonal entries
			if (MACGrid_->is_fluid(i, j)) {
				// x-adjacent cells
				if (i+1 < nx && MACGrid_->is_fluid(i+1, j)) {
						triplets.push_back(Mac2d::Triplet_t(cellidx, cellidx+1, -1));
						// Use symmetry to avoid computing (i-1,j) separately
						triplets.push_back(Mac2d::Triplet_t(cellidx+1, cellidx, -1));
				}
				// y-adjacent cells
				if (j+1 < ny && MACGrid_->is_fluid(i, j+1)) {
						triplets.push_back(Mac2d::Triplet_t(cellidx, cellidx + nx, -1));
						// Use symmetry to avoid computing (i,j-1) separately
						triplets.push_back(Mac2d::Triplet_t(cellidx + nx, cellidx, -1));
				}
			} // if is_fluid(i,j)
		}
	} // outer for

	A_.resize(nx*ny, nx*ny);
	A_.setZero();
	A_.setFromTriplets(triplets.begin(), triplets.end());
}

void FLIP::compute_pressure_rhs(const double dt) {
	// Compute right-hand side of the pressure equations and store in d_
	//  See eq. (4.19) and (4.24) in SIGGRAPH notes
	//   Note: u_{solid} = 0
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();

	// Alias for MAC grid
	auto& g = MACGrid_;
	
	d_.resize(nx*ny);
	d_.setZero();
	
	unsigned cellidx = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++, cellidx++) {
			if (g->is_fluid(i,j)) {
				// get_u(i,j) = u_{ (i-1/2, j) }
				double d_ij = -(g->get_u(i+1,j) - g->get_u(i,j));
				d_ij -= g->get_v(i,j+1) - g->get_v(i,j);
				
				// Note: u_{solid} = 0
	
				// Check each adjacent cell. If solid, alter term as in (4.24)
				// Consider cells outside of the boundary as solid
				// (i+1, j)
				if ((i < (nx-1) && g->is_solid(i+1,j)) || i == nx-1) {
					d_ij += g->get_u(i+1,j);
				}
				
				// (i-1, j)
				if ((i > 0 && g->is_solid(i-1,j)) || i == 0) {
					d_ij += g->get_u(i,j);
				}
	
				// (i, j+1)
				if ((j < (ny-1) && g->is_solid(i,j+1)) || j == ny-1) {
					d_ij += g->get_v(i,j+1);
				}
				
				// (i, j-1)
				if ((j > 0 && g->is_solid(i,j-1)) || j == 0) {
					d_ij += g->get_v(i,j);
				}
				
				d_(cellidx) = fluid_density_ * g->get_cell_sizex() * d_ij / dt;
			} else { // if is_fluid(i,j)
				d_(cellidx) = 0;
			}
		}
	}
}

void FLIP::apply_pressure_gradients(const double dt) {
	// Apply pressure gradients to velocity field
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	// MACGrid alias
	auto& g = MACGrid_;
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
	// FLIP: alpha = 0.
	// PIC: alpha = 1.
	double alpha = alpha_;
	int nx = MACGrid_->get_num_cells_x();
	int ny = MACGrid_->get_num_cells_y();
	
	for(unsigned i = 0; i < num_particles_; ++i){
		//Store the initial positions and velocities of the particles
		Eigen::Vector3d initial_position = (particles_+i)->get_position();
		Eigen::Vector3d initial_velocity = (particles_+i)->get_velocity();
		auto initial_idx = MACGrid_->index_from_coord(initial_position(0), 
													  initial_position(1));
		
		//Initialization of the variables
		Eigen::Vector3d interp_u_star;
		interp_u_star.setZero();
		Eigen::Vector3d interp_u_n1;
		interp_u_n1.setZero();
		Eigen::Vector3d u_update;
		u_update.setZero();

		double x = initial_position[0];
		double y = initial_position[1];
		
		//With u* and v* we can make the interpolation interp(u*, x_p),
		//with the new u and v we can make the interpolation interp(u_n1, x_p)
		
		//Update the u-velocity (bilinear interpolation)
		interp_u_star[0] = MACGrid_->get_interp_u_star(x,y);
		interp_u_n1[0] = MACGrid_->get_interp_u(x,y);
		
		//Update the v-velocity (bilinear interpolation)
		interp_u_star[1] = MACGrid_->get_interp_v_star(x,y);
		interp_u_n1[1] = MACGrid_->get_interp_v(x,y);
		

		//Update the final velocity of the particles
		
		// Use pure PIC on boundary, blend PIC+FLIP elsewhere
		if (initial_idx.first == 0 or initial_idx.first == nx-1
		 or initial_idx.second == 0 or initial_idx.second == ny-1){
			//~ u_update = initial_velocity*(1 - std::min(1., 2*alpha)) + interp_u_n1 - interp_u_star*(1 - std::min(1., 2*alpha));
			u_update = initial_velocity*(1 - alpha) + interp_u_n1 - interp_u_star*(1 - alpha);
		} else {
			u_update = initial_velocity*(1 - alpha) + interp_u_n1 - interp_u_star*(1 - alpha);
		}
		
		(particles_ + i)->set_velocity(u_update);
	}
}


void FLIP::advance_particles(const double dt, const unsigned step) {
	for(unsigned n = 0; n < num_particles_; ++n){
		Eigen::Vector3d pos_curr = (particles_ + n)->get_position();
		Eigen::Vector3d vel = (particles_ + n)->get_velocity();
		
		Eigen::Vector3d pos_next;
		
		// Euler estimate
		Eigen::Vector3d pos_half = pos_curr + 0.5*dt*vel;
		
		double size_x = MACGrid_->get_grid_size()(0);
		double size_y = MACGrid_->get_grid_size()(1);
		double cell_sizex = MACGrid_->get_cell_sizex();
		double cell_sizey = MACGrid_->get_cell_sizey();
		double x_half = pos_half(0);
		double y_half = pos_half(1);
		
		// Check if pos_half is out of the grid
		if ((x_half <= -0.5*cell_sizex) or (x_half >= size_x - 0.5*cell_sizex) or (y_half <= -0.5*cell_sizey) or (y_half >= size_y - 0.5*cell_sizex)){
			continue;
		}
		
		// RK2
		pos_next(0) = pos_curr(0) + dt*MACGrid_->get_interp_u(x_half, y_half);
		pos_next(1) = pos_curr(1) + dt*MACGrid_->get_interp_v(x_half, y_half);
		
		double x = pos_next(0);
		double y = pos_next(1);
		
		// Check if the particle exits the grid
		if (x <= -0.5*cell_sizex) {
			pos_next(0) = 0.;
		}
		if (x >= size_x - 0.5*cell_sizex) {
			pos_next(0) = size_x - cell_sizex;
		}
		if (y <= -0.5*cell_sizey) {
			pos_next(1) = 0.;
		}
		if (y >= size_y - 0.5*cell_sizey) {
			pos_next(1) = size_y - cell_sizey;
		}

		(particles_ + n)->set_position(pos_next);
	}
}

void FLIP::explode(const double dt, const unsigned long step, const double value){
	// Apply external forces to simulate meteorite crash
	const double alpha = 2;
	const double beta = 5;
	const double force = value;
	const double spin = -0.25*value;
	const double sty = step - 5;
	const double stx = 3*sty;
	
	// Coordinates of the bottom left point of the square containing the meteorite
	const int x_bl = 46;
	const int y_bl = 15;
	
	// Bottom
	MACGrid_->set_v(x_bl+1-stx, y_bl-sty, MACGrid_->get_v(x_bl+1-stx,y_bl-sty) - dt*force*alpha);
	MACGrid_->set_v(x_bl+2-stx, y_bl-sty, MACGrid_->get_v(x_bl+2-stx,y_bl-sty) - dt*force*alpha);
	
	// Bottom Half
	MACGrid_->set_v(x_bl-stx, y_bl+1-sty, MACGrid_->get_v(x_bl-stx,y_bl+1-sty) - dt*force*alpha);
	MACGrid_->set_v(x_bl+1-stx, y_bl+1-sty, MACGrid_->get_v(x_bl+1-stx,y_bl+1-sty) - dt*force*beta);
	MACGrid_->set_v(x_bl+2-stx, y_bl+1-sty, MACGrid_->get_v(x_bl+2-stx,y_bl+1-sty) - dt*force*beta);
	MACGrid_->set_v(x_bl+3-stx, y_bl+1-sty, MACGrid_->get_v(x_bl+3-stx,y_bl+1-sty) - dt*force*alpha);
	
	// Half
	MACGrid_->set_v(x_bl-stx, y_bl+2-sty, MACGrid_->get_v(x_bl-stx,y_bl+2-sty) - dt*spin);
	MACGrid_->set_v(x_bl+1-stx, y_bl+2-sty, MACGrid_->get_v(x_bl+1-stx,y_bl+2-sty) - dt*spin*beta);
	MACGrid_->set_v(x_bl+2-stx, y_bl+2-sty, MACGrid_->get_v(x_bl+2-stx,y_bl+2-sty) + dt*spin*beta);
	MACGrid_->set_v(x_bl+3-stx, y_bl+2-sty, MACGrid_->get_v(x_bl+3-stx,y_bl+2-sty) + dt*spin);
	
	// Top Half
	MACGrid_->set_v(x_bl-stx, y_bl+3-sty, MACGrid_->get_v(x_bl-stx,y_bl+3-sty) + dt*force);
	MACGrid_->set_v(x_bl+1-stx, y_bl+3-sty, MACGrid_->get_v(x_bl+1-stx,y_bl+3-sty) + dt*force*beta);
	MACGrid_->set_v(x_bl+2-stx, y_bl+3-sty, MACGrid_->get_v(x_bl+2-stx,y_bl+3-sty) + dt*force*beta);
	MACGrid_->set_v(x_bl+3-stx, y_bl+3-sty, MACGrid_->get_v(x_bl+3-stx,y_bl+3-sty) + dt*force);
	
	// Top
	MACGrid_->set_v(x_bl+1-stx, y_bl+4-sty, MACGrid_->get_v(x_bl+1-stx,y_bl+4-sty) + dt*force);
	MACGrid_->set_v(x_bl+2-stx, y_bl+4-sty, MACGrid_->get_v(x_bl+2-stx,y_bl+4-sty) + dt*force);
	
	// Left
	MACGrid_->set_u(x_bl-stx, y_bl+1-sty, MACGrid_->get_u(x_bl-stx,y_bl+1-sty) - dt*force*alpha);
	MACGrid_->set_u(x_bl-stx, y_bl+2-sty, MACGrid_->get_u(x_bl-stx,y_bl+2-sty) - dt*force*alpha);
	
	// Left Half
	MACGrid_->set_u(x_bl+1-stx, y_bl-sty, MACGrid_->get_u(x_bl+1-stx,y_bl-sty) - dt*force*alpha);
	MACGrid_->set_u(x_bl+1-stx, y_bl+1-sty, MACGrid_->get_u(x_bl+1-stx,y_bl+1-sty) - dt*force*beta);
	MACGrid_->set_u(x_bl+1-stx, y_bl+2-sty, MACGrid_->get_u(x_bl+1-stx,y_bl+2-sty) - dt*force*beta);
	MACGrid_->set_u(x_bl+1-stx, y_bl+3-sty, MACGrid_->get_u(x_bl+1-stx,y_bl+3-sty) - dt*force*alpha);
	
	// Half
	MACGrid_->set_u(x_bl+2-stx, y_bl-sty, MACGrid_->get_u(x_bl+2-stx,y_bl-sty) + dt*spin);
	MACGrid_->set_u(x_bl+2-stx, y_bl+1-sty, MACGrid_->get_u(x_bl+2-stx,y_bl+1-sty) + dt*spin*beta);
	MACGrid_->set_u(x_bl+2-stx, y_bl+2-sty, MACGrid_->get_u(x_bl+2-stx,y_bl+2-sty) - dt*spin*beta);
	MACGrid_->set_u(x_bl+2-stx, y_bl+3-sty, MACGrid_->get_u(x_bl+2-stx,y_bl+3-sty) - dt*spin);
	
	// Right Half
	MACGrid_->set_u(x_bl+3-stx, y_bl-sty, MACGrid_->get_u(x_bl+3-stx,y_bl-sty) + dt*force);
	MACGrid_->set_u(x_bl+3-stx, y_bl+1-sty, MACGrid_->get_u(x_bl+3-stx,y_bl+1-sty) + dt*force*beta);
	MACGrid_->set_u(x_bl+3-stx, y_bl+2-sty, MACGrid_->get_u(x_bl+3-stx,y_bl+2-sty) + dt*force*beta);
	MACGrid_->set_u(x_bl+3-stx, y_bl+3-sty, MACGrid_->get_u(x_bl+3-stx,y_bl+3-sty) + dt*force);
	
	// Right
	MACGrid_->set_u(x_bl+4-stx, y_bl+1-sty, MACGrid_->get_u(x_bl+4-stx,y_bl+1-sty) + dt*force);
	MACGrid_->set_u(x_bl+4-stx, y_bl+2-sty, MACGrid_->get_u(x_bl+4-stx,y_bl+2-sty) + dt*force);
}
