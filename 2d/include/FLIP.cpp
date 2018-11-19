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
	do_pressures(dt);

	// 4.
	grid_to_particle();

	// 5.
	advance_particles();
}

/*** COMPUTE VELOCITY FIELD ***/
void FLIP::compute_velocity_field() {
	// TODO: 1. Compute the velocity field (velocities on grid)
	//  1a. particle-to-grid transfer
	//  1b. classify nonsolid cells as fluid or air
	//  1c. extrapolate velocity field to air cells
	//    -> see SIGGRAPH §6.3
}

/*** APPLY EXTERNAL FORCES ***/
void FLIP::apply_forces() {
	// TODO: compute&apply external forces (gravity, vorticity confinement, ...) 
	// Apply them to the velocity field via forward euler
	// Only worry about gravity for now
}

/*** PRESSURE SOLVING ***/
void FLIP::do_pressures(const double dt) {
	// 3. Compute & apply pressure gradients to field
	
	//   3a. Compute A matrix
	compute_pressure_matrix();

	//   3b. Compute rhs d
	compute_pressure_rhs(dt);

	//   3c. Solve for p: Ap = d (MICCG(0))
	// *TODO: use MICCG(0) solver
	using namespace Eigen;
	using solver_t = ConjugateGradient<SparseMatrix<double>, Lower|Upper>;
	solver_t solver;
	solver.compute(A_);
	VectorXd p = solver.solve(d_);

	// Copy pressures to MAC grid
	MACGrid_->set_pressure(p);

	//  3d. Apply pressure gradients to velocity field
	//     -> see SIGGRAPH §4
	apply_pressure_gradients(dt);

	// Note: boundary conditions are handles here by setting the pressures
	//    such that no particle exits the system.
}

void FLIP::compute_pressure_matrix() {
	// Compute matrix for pressure solve and store in A_
	//  See eq. (4.19) and (4.24) in SIGGRAPH notes
	
	std::vector< Mac2d::Triplet_t > triplets;
	
	unsigned nx = MACGrid_->get_cell_sizex();
	unsigned ny = MACGrid_->get_cell_sizey();

	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			if (MACGrid_->is_fluid(i, j)) {
				// Copy diagonal entry
				triplets.push_back(MACGrid_->get_a_diag()[i + j*nx]);
				
				// Compute off-diagonal entries
				// x-adjacent cells
				if (i < nx-1 && MACGrid_->is_fluid(i+1, j)) {
						triplets.push_back(Mac2d::Triplet_t(i+j*ny, i+1+j*ny, -1));
						// Use symmetry to avoid computing (i-1,j) separately
						triplets.push_back(Mac2d::Triplet_t(i+1+j*ny, i+j*ny, -1));
				}
				// y-adjacent cells
				if (j < ny-1 && MACGrid_->is_fluid(i, j+1)) {
						triplets.push_back(Mac2d::Triplet_t(i+j*ny, i+(j+1)*ny, -1));
						// Use symmetry to avoid computing (i,j-1) separately
						triplets.push_back(Mac2d::Triplet_t(i+(j+1)*ny, i+j*ny, -1));
				}
			} // if is_fluid(i,j)
		}
	} // outer for

	A_.setZero();
	A_.setFromTriplets(triplets.begin(), triplets.end());
}

void FLIP::compute_pressure_rhs(const double dt) {
	// Compute right-hand side of the pressure equations and store in d_
	//  See eq. (4.19) and (4.24) in SIGGRAPH notes
	d_.setZero();

	unsigned nx = MACGrid_->get_cell_sizex();
	unsigned ny = MACGrid_->get_cell_sizey();
	// Alias for MAC grid
	auto& g = MACGrid_;

	unsigned idx = 0;
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
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
			d_(idx) = fluid_density_ * dt * d_ij / g->get_cell_sizex();
			idx++;
		}
	}
}

void FLIP::apply_pressure_gradients(const double dt) {
	// TODO: apply pressure gradients to velocity field
	unsigned nx = MACGrid_->get_cell_sizex();
	unsigned ny = MACGrid_->get_cell_sizey();
	for (unsigned j = 1; j < ny; j++) {
		for (unsigned i = 1; i < nx; i++) {
			
		}
	}
}



/*** UPDATE PARTICLE VELOCITIES & MOVE PARTICLES ***/
void FLIP::grid_to_particle() {
	// TODO: FLIP grid to particle transfer
	//  -> See slides Fluids II, FLIP_explained.pdf
}

void FLIP::advance_particles() {
	// TODO: update partimacle positions 
	//  - Use RK2 interpolator
}

