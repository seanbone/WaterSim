#include "FLIP.h"
#include "tsc_x86.hpp"

FLIP::FLIP(Particle* particles, unsigned num_particles, Mac3d* MACGrid,
		   const SimConfig& cfg)
	: particles_(particles), num_particles_(num_particles), MACGrid_(MACGrid), cfg_(cfg),
	  fluid_density_(cfg.getDensity()), gravity_mag_(cfg.getGravity()), alpha_(cfg.getAlpha()) {
	
}


/*** PERFORM ONE STEP ***/
void FLIP::step_FLIP(double dt, unsigned long step) {
	/** One FLIP step:
	 * 1. Compute velocity field (particle-to-grid transfer)
	 *    - Particle-to-grid transfer
	 *    - Classify cells (fluid/air)
	 *    - Extrapolate velocity field into air region
	 * 1a. Copy velocity field to intermediate velocity field u^*
	 * 2. Apply external forces (fwd euler on field)
	 * 3. Enforce boundary conditions for grid & solid boundaries
	 * 4. Compute & apply pressure gradients
	 * 5. Update particle velocities from grid-velocities
	 * 6. Subsample time interval to satisfy CFL condition
	 * 7. Update particle positions
	 */

	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");

	// 1.
	tsctimer.start_timing("velocity_field");
	compute_velocity_field();

	// 1a.
	MACGrid_->set_uvw_star();
	tsctimer.stop_timing("velocity_field", true, "");

	// 2.
	tsctimer.start_timing("apply_forces");
	apply_forces(dt);
	
	if (cfg_.getApplyMeteorForce() && step <= 200 ){
		explode(dt, step, 15, 0, 15, 2, 800);
	}

	tsctimer.stop_timing("apply_forces", true, "");

	// 3.
	tsctimer.start_timing("apply_boundary_conditions");
	apply_boundary_conditions();
	tsctimer.stop_timing("apply_boundary_conditions", true, "");
	
	// 4.
	tsctimer.start_timing("do_pressures");
	do_pressures(dt);
	tsctimer.stop_timing("do_pressures", true, "");
	
	// 5.
	tsctimer.start_timing("grid_to_particle");
	grid_to_particle();
	tsctimer.stop_timing("grid_to_particle", true, "");

	// 6.
	tsctimer.start_timing("advance_particles");
	double dt_new = compute_timestep(dt);
	double num_substeps = std::ceil(dt/dt_new);
	for( int s = 0; s < num_substeps ; ++s ){
		
		// 7.
		advance_particles(dt/num_substeps, step);
	}
	tsctimer.stop_timing("advance_particles", true, "");
}


/*** COMPUTE VELOCITY FIELD ***/
void FLIP::compute_velocity_field() {

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
	double cell_sizez = MACGrid_->get_cell_sizez();
	
	// Threshold h and h_scaled so that it is equal to the distance 
	// expressed in number of cells
	double h = 2*cell_sizex;
	int h_scaledx = std::ceil(h/cell_sizex);
	int h_scaledy = std::ceil(h/cell_sizey);
	int h_scaledz = std::ceil(h/cell_sizez);
	
	// Lists of flags for visited grid-velocities: 1 -> visited
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned L = MACGrid_->get_num_cells_z();
	bool* visited_u = new bool[M*L*(N+1)];
	bool* visited_v = new bool[N*L*(M+1)];
	bool* visited_w = new bool[N*M*(L+1)];
	std::fill(visited_u, visited_u + M*L*(N+1), 0);
	std::fill(visited_v, visited_v + N*L*(M+1), 0);
	std::fill(visited_w, visited_w + N*M*(L+1), 0);
	
	// Reset all fluid flags
	MACGrid_->reset_fluid();
	
	// Iterate over all particles and add weighted particle velocities
	// to grid points within a threshold h (in this case equal to double
	// the length of an edge of a cell)
	for( unsigned int n = 0; n < num_particles_; ++n ){
		
		// Get the position and velocity vector of the current particle
		pos = (particles_ + n)->get_position();
		vel = (particles_ + n)->get_velocity();

		// Get the indices corresponding to the cell containing the 
		// current particle
		Eigen::Vector3d tmp = MACGrid_->index_from_coord(pos(0), pos(1), pos(2));
		cell_coord << tmp(0), tmp(1), tmp(2);
		
		// Set the cell of the current particle to a fluid-cell
		if ( !(MACGrid_->is_fluid(cell_coord(0), cell_coord(1), cell_coord(2))) and !(MACGrid_->is_solid(cell_coord(0), cell_coord(1), cell_coord(2))) ){
			MACGrid_->set_fluid(cell_coord(0), cell_coord(1), cell_coord(2));
		}
		
		// Coordinates of the points on the grid faces (in meters)
		Eigen::Vector3d grid_coord;
		grid_coord << 0, 0, 0;
		
		// Get total number of cells on each axis
		int nx = MACGrid_->get_num_cells_x();
		int ny = MACGrid_->get_num_cells_y();
		int nz = MACGrid_->get_num_cells_z();
		
		// For each particle iterate only over the grid-velocities in a
		// h_scaled neighborhood
		for( int k = cell_coord(2) - h_scaledz; k <= cell_coord(2) + h_scaledz + 1; ++k ){	
			for( int j = cell_coord(1) - h_scaledy; j <= cell_coord(1) + h_scaledy + 1; ++j ){
				for( int i = cell_coord(0) - h_scaledx; i <= cell_coord(0) + h_scaledx + 1; ++i ){
					if ( ( i >= 0 and j >= 0 and k >= 0 ) ){
						if ( ( i <= nx and j < ny and k < nz ) ){
						
							// Left Face
							grid_coord(0) = (i - 0.5) * cell_sizex;
							grid_coord(1) = j * cell_sizey;
							grid_coord(2) = k * cell_sizez;
							accumulate_u(pos, vel, grid_coord, h, i, j, k);
						}
						
						if ( ( i < nx and j <= ny and k < nz ) ){
							
							// Lower Face
							grid_coord(0) = i * cell_sizex;
							grid_coord(1) = (j - 0.5) * cell_sizey;
							grid_coord(2) = k * cell_sizez;
							accumulate_v(pos, vel, grid_coord, h, i, j, k);
						}						
						
						if ( ( i < nx and j < ny and k <= nz ) ){
							
							// Farthest Face (the closer to the origin)
							grid_coord(0) = i * cell_sizex;
							grid_coord(1) = j * cell_sizey;
							grid_coord(2) = (k - 0.5) * cell_sizez;
							accumulate_w(pos, vel, grid_coord, h, i, j, k);
						}
					}
				}
			}
		}
	}
	
	// Normalize grid-velocities
	normalize_accumulated_u( visited_u );
	normalize_accumulated_v( visited_v );
	normalize_accumulated_w( visited_w );
	
	// Extrapolate velocities
	extrapolate_u( visited_u );
	extrapolate_v( visited_v );
	extrapolate_w( visited_w );
	
	// Clear the Heap
	delete[] visited_u;
	delete[] visited_v;
	delete[] visited_w;
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
	// Distance between the particle and the location on which the 
	// grid-velocity is saved (the center of a face)
	double r = (particle_coord - grid_coord).norm();
	
	// Compute h^9 (std::pow() is inefficient)
	double h2 = h*h;
	double h4 = h2 * h2;
	double h9 = h4 * h4 * h;
	
	double diff = h2 - r*r;
	double diff3 = diff*diff*diff;
	
	double coeff = 315/(64 * M_PI * h9);
	
	return coeff * diff3;
}

				
void FLIP::accumulate_u( const Eigen::Vector3d& pos,
						 const Eigen::Vector3d& vel,
						 const Eigen::Vector3d& grid_coord,
						 const double h,
						 const int i,
						 const int j,
						 const int k )
{	
	// If the current particle is within the given threshold h, update 
	// the horizontal grid-velocity at grid_coord with the weighted 
	// particle-velocity
	if ( check_threshold(pos, grid_coord, h) ){
		double u_prev = MACGrid_->get_u(i, j, k);
		double W_u = compute_weight(pos, grid_coord, h);
		double u_curr = u_prev + (W_u * vel(0));
		
		// Accumulate velocities
		MACGrid_->set_u(i, j, k, u_curr);
		
		// Accumulate weights
		double W_u_prev = MACGrid_->get_weights_u(i, j, k);
		double W_u_curr = W_u_prev + W_u;
		MACGrid_->set_weights_u(i, j, k, W_u_curr);
	}
}


void FLIP::accumulate_v( const Eigen::Vector3d& pos,
						 const Eigen::Vector3d& vel,
						 const Eigen::Vector3d& grid_coord,
						 const double h,
						 const int i,
						 const int j,
						 const int k )
{
	// If the current particle is within the given threshold h, update 
	// the vertical grid-velocity at grid_coord with the weighted 
	// particle-velocity
	if ( check_threshold(pos, grid_coord, h) ){
		double v_prev = MACGrid_->get_v(i, j, k);
		double W_v = compute_weight(pos, grid_coord, h);
		double v_curr = v_prev + (W_v * vel(1));
		
		// Accumulate velocities
		MACGrid_->set_v(i, j, k, v_curr);
		
		// Accumulate weights
		double W_v_prev = MACGrid_->get_weights_v(i, j, k);
		double W_v_curr = W_v_prev + W_v;
		MACGrid_->set_weights_v(i, j, k, W_v_curr);
	}
}


void FLIP::accumulate_w( const Eigen::Vector3d& pos,
						 const Eigen::Vector3d& vel,
						 const Eigen::Vector3d& grid_coord,
						 const double h,
						 const int i,
						 const int j,
						 const int k )
{
	// If the current particle is within the given threshold h, update 
	// the outgoing grid-velocity at grid_coord with the weighted 
	// particle-velocity
	if ( check_threshold(pos, grid_coord, h) ){
		double w_prev = MACGrid_->get_w(i, j, k);
		double W_w = compute_weight(pos, grid_coord, h);
		double w_curr = w_prev + (W_w * vel(2));
		
		// Accumulate velocities
		MACGrid_->set_w(i, j, k, w_curr);
		
		// Accumulate weights
		double W_w_prev = MACGrid_->get_weights_w(i, j, k);
		double W_w_curr = W_w_prev + W_w;
		MACGrid_->set_weights_w(i, j, k, W_w_curr);
	}
}


void FLIP::normalize_accumulated_u( bool* const visited_u ){
	
	// Get total number of cells on each axis
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned L = MACGrid_->get_num_cells_z();
	
	// Iterate over all horizontal grid-velocities and divide for the 
	// corresponding weight (if non-zero). Also set the flags of 
	// visited_u to 1 if the weight is non-zero (-> visited)
	for( unsigned k = 0; k < L; ++k ){
		for( unsigned j = 0; j < M; ++j ){
			for( unsigned i = 0; i < N+1; ++i ){
				
				// Get weight of current grid-velocity
				double W_u = MACGrid_->get_weights_u(i, j, k);
				
				if ( W_u != 0 ){
					double u_prev = MACGrid_->get_u(i, j, k);
					double u_curr = u_prev/W_u;
					MACGrid_->set_u(i, j, k, u_curr);
					*(visited_u + (N+1)*j + i + (N+1)*M*k) = 1;
				}
			}	
		}
	}	
}


void FLIP::normalize_accumulated_v( bool* const visited_v ){
	
	// Get total number of cells on each axis
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned L = MACGrid_->get_num_cells_z();
	
	// Iterate over all vertical grid-velocities and divide for the 
	// corresponding weight (if non-zero). Also set the flags of 
	// visited_v to 1 if the weight is non-zero (-> visited)
	for( unsigned k = 0; k < L; ++k ){	
		for( unsigned j = 0; j < M+1; ++j ){
			for( unsigned i = 0; i < N; ++i ){
				
				// Get weight of current grid-velocity
				double W_v = MACGrid_->get_weights_v(i, j, k);
				
				if ( W_v != 0 ){
					double v_prev = MACGrid_->get_v(i, j, k);
					double v_curr = v_prev/W_v;
					MACGrid_->set_v(i, j, k, v_curr);
					*(visited_v + N*j + i + N*(M+1)*k) = 1;
				}
			}	
		}	
	}
}


void FLIP::normalize_accumulated_w( bool* const visited_w ){
	
	// Get total number of cells on each axis
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned L = MACGrid_->get_num_cells_z();
	
	// Iterate over all outgoing grid-velocities and divide for the 
	// corresponding weight (if non-zero). Also set the flags of 
	// visited_w to 1 if the weight is non-zero (-> visited)
	for( unsigned k = 0; k < L+1; ++k ){
		for( unsigned j = 0; j < M; ++j ){
			for( unsigned i = 0; i < N; ++i ){
				
				// Get weight of current grid-velocity
				double W_w = MACGrid_->get_weights_w(i, j, k);
				
				if ( W_w != 0 ){
					double w_prev = MACGrid_->get_w(i, j, k);
					double w_curr = w_prev/W_w;
					MACGrid_->set_w(i, j, k, w_curr);
					*(visited_w + N*j + i + N*M*k) = 1;
				}
			}	
		}	
	}
}


void FLIP::extrapolate_u( const bool* const visited_u ){
	
	// Get total number of cells on each axis
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned L = MACGrid_->get_num_cells_z();
	
	// Count the times a specific grid-velocity is accessed (to average 
	// the neighboring velocities)
	unsigned* counter = new unsigned[M*L*(N+1)];
	std::fill(counter, counter + M*L*(N+1), 0);
	
	// Iterate over all horizontal grid-velocities and extrapolate into 
	// the air cells (not visited) the average velocities of the 
	// neighboring fluid/visited cells
	for( unsigned k = 0; k < L; ++k ){	
		for( unsigned j = 0; j < M; ++j ){
			for( unsigned i = 0; i < N+1; ++i ){
				
				if ( *(visited_u + (N+1)*j + i + (N+1)*M*k) ){
					
					if ( i != 0 and !(*(visited_u + (N+1)*j + (i-1) + (N+1)*M*k)) ){
						double tmp = MACGrid_->get_u(i-1, j, k) * *(counter + (N+1)*j + (i-1) + (N+1)*M*k);
						*(counter + (N+1)*j + (i-1) + (N+1)*M*k) += 1;
						MACGrid_->set_u(i-1, j, k, (tmp + MACGrid_->get_u(i, j, k))/(*(counter + (N+1)*j + (i-1) + (N+1)*M*k)));
					}
					
					if ( j != 0 and !(*(visited_u + (N+1)*(j-1) + i + (N+1)*M*k)) ){
						double tmp = MACGrid_->get_u(i, j-1, k) * *(counter + (N+1)*(j-1) + i + (N+1)*M*k);
						*(counter + (N+1)*(j-1) + i + (N+1)*M*k) += 1;
						MACGrid_->set_u(i, j-1, k, (tmp + MACGrid_->get_u(i, j, k))/(*(counter + (N+1)*(j-1) + i + (N+1)*M*k)));
					}
					
					if ( k != 0 and !(*(visited_u + (N+1)*j + i + (N+1)*M*(k-1))) ){
						double tmp = MACGrid_->get_u(i, j, k-1) * *(counter + (N+1)*j + i + (N+1)*M*(k-1));
						*(counter + (N+1)*j + i + (N+1)*M*(k-1)) += 1;
						MACGrid_->set_u(i, j, k-1, (tmp + MACGrid_->get_u(i, j, k))/(*(counter + (N+1)*j + i + (N+1)*M*(k-1))));
					}
					
					if ( i != N and !(*(visited_u + (N+1)*j + (i+1) + (N+1)*M*k)) ){
						double tmp = MACGrid_->get_u(i+1, j, k) * *(counter + (N+1)*j + (i+1) + (N+1)*M*k);
						*(counter + (N+1)*j + (i+1) + (N+1)*M*k) += 1;
						MACGrid_->set_u(i+1, j, k, (tmp + MACGrid_->get_u(i, j, k))/(*(counter + (N+1)*j + (i+1) + (N+1)*M*k)));
					}
					
					if ( j != M-1 and !(*(visited_u + (N+1)*(j+1) + i + (N+1)*M*k)) ){
						double tmp = MACGrid_->get_u(i, j+1, k) * *(counter + (N+1)*(j+1) + i + (N+1)*M*k);
						*(counter + (N+1)*(j+1) + i + (N+1)*M*k) += 1;
						MACGrid_->set_u(i, j+1, k, (tmp + MACGrid_->get_u(i, j, k))/(*(counter + (N+1)*(j+1) + i + (N+1)*M*k)));
					}
					
					if ( k != L-1 and !(*(visited_u + (N+1)*j + i + (N+1)*M*(k+1))) ){
						double tmp = MACGrid_->get_u(i, j, k+1) * *(counter + (N+1)*j + i + (N+1)*M*(k+1));
						*(counter + (N+1)*j + i + (N+1)*M*(k+1)) += 1;
						MACGrid_->set_u(i, j, k+1, (tmp + MACGrid_->get_u(i, j, k))/(*(counter + (N+1)*j + i + (N+1)*M*(k+1))));
					}
				}
			}
		}
	}
	
	// Clear the Heap
	delete[] counter;
}

void FLIP::extrapolate_v( const bool* const visited_v ){
	
	// Get total number of cells on each axis
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned L = MACGrid_->get_num_cells_z();
	
	// Count the times a specific grid-velocity is accessed (to average 
	// the neighboring velocities)
	unsigned* counter = new unsigned[N*L*(M+1)];
	std::fill(counter, counter + N*L*(M+1), 0);
	
	// Iterate over all vertical grid-velocities and extrapolate into 
	// the air cells (not visited) the average velocities of the 
	// neighboring fluid/visited cells
	for( unsigned k = 0; k < L; ++k ){	
		for( unsigned j = 0; j < M+1; ++j ){
			for( unsigned i = 0; i < N; ++i ){
				
				if ( *(visited_v + N*j + i + N*(M+1)*k) ){
					
					if ( i != 0 and !(*(visited_v + N*j + (i-1) + N*(M+1)*k)) ){
						double tmp = MACGrid_->get_v(i-1, j, k) * *(counter + N*j + (i-1) + N*(M+1)*k);
						*(counter + N*j + (i-1) + N*(M+1)*k) += 1;
						MACGrid_->set_v(i-1, j, k, (tmp + MACGrid_->get_v(i, j, k))/(*(counter + N*j + (i-1) + N*(M+1)*k)));
					}
					
					if ( j != 0 and !(*(visited_v + N*(j-1) + i + N*(M+1)*k)) ){
						double tmp = MACGrid_->get_v(i, j-1, k) * *(counter + N*(j-1) + i + N*(M+1)*k);
						*(counter + N*(j-1) + i + N*(M+1)*k) += 1;
						MACGrid_->set_v(i, j-1, k, (tmp + MACGrid_->get_v(i, j, k))/(*(counter + N*(j-1) + i + N*(M+1)*k)));
					}
					
					if ( k != 0 and !(*(visited_v + N*j + i + N*(M+1)*(k-1))) ){
						double tmp = MACGrid_->get_v(i, j, k-1) * *(counter + N*j + i + N*(M+1)*(k-1));
						*(counter + N*j + i + N*(M+1)*(k-1)) += 1;
						MACGrid_->set_v(i, j, k-1, (tmp + MACGrid_->get_v(i, j, k))/(*(counter + N*j + i + N*(M+1)*(k-1))));
					}
					
					if ( i != N-1 and !(*(visited_v + N*j + (i+1) + N*(M+1)*k)) ){
						double tmp = MACGrid_->get_v(i+1, j, k) * *(counter + N*j + (i+1) + N*(M+1)*k);
						*(counter + N*j + (i+1) + N*(M+1)*k) += 1;
						MACGrid_->set_v(i+1, j, k, (tmp + MACGrid_->get_v(i, j, k))/(*(counter + N*j + (i+1) + N*(M+1)*k)));
					}
					
					if ( j != M and !(*(visited_v + N*(j+1) + i + N*(M+1)*k)) ){
						double tmp = MACGrid_->get_v(i, j+1, k) * *(counter + N*(j+1) + i + N*(M+1)*k);
						*(counter + N*(j+1) + i + N*(M+1)*k) += 1;
						MACGrid_->set_v(i, j+1, k, (tmp + MACGrid_->get_v(i, j, k))/(*(counter + N*(j+1) + i + N*(M+1)*k)));
					}
					
					if ( k != L-1 and !(*(visited_v + N*j + i + N*(M+1)*(k+1))) ){
						double tmp = MACGrid_->get_v(i, j, k+1) * *(counter + N*j + i + N*(M+1)*(k+1));
						*(counter + N*j + i + N*(M+1)*(k+1)) += 1;
						MACGrid_->set_v(i, j, k+1, (tmp + MACGrid_->get_v(i, j, k))/(*(counter + N*j + i + N*(M+1)*(k+1))));
					}
				}
			}
		}
	}
	
	// Clear the Heap
	delete[] counter;
}

void FLIP::extrapolate_w( const bool* const visited_w ){
	
	// Get total number of cells on each axis
	unsigned M = MACGrid_->get_num_cells_y();
	unsigned N = MACGrid_->get_num_cells_x();
	unsigned L = MACGrid_->get_num_cells_z();
	
	// Count the times a specific grid-velocity is accessed (to average 
	// the neighboring velocities)
	unsigned* counter = new unsigned[N*M*(L+1)];
	std::fill(counter, counter + N*M*(L+1), 0);
	
	// Iterate over all outgoing grid-velocities and extrapolate into 
	// the air cells (not visited) the average velocities of the 
	// neighboring fluid/visited cells
	for( unsigned k = 0; k < L+1; ++k ){	
		for( unsigned j = 0; j < M; ++j ){
			for( unsigned i = 0; i < N; ++i ){
				
				if ( *(visited_w + N*j + i + N*M*k) ){
					
					if ( i != 0 and !(*(visited_w + N*j + (i-1) + N*M*k)) ){
						double tmp = MACGrid_->get_w(i-1, j, k) * *(counter + N*j + (i-1) + N*M*k);
						*(counter + N*j + (i-1) + N*M*k) += 1;
						MACGrid_->set_w(i-1, j, k, (tmp + MACGrid_->get_w(i, j, k))/(*(counter + N*j + (i-1) + N*M*k)));
					}
					
					if ( j != 0 and !(*(visited_w + N*(j-1) + i + N*M*k)) ){
						double tmp = MACGrid_->get_w(i, j-1, k) * *(counter + N*(j-1) + i + N*M*k);
						*(counter + N*(j-1) + i + N*M*k) += 1;
						MACGrid_->set_w(i, j-1, k, (tmp + MACGrid_->get_w(i, j, k))/(*(counter + N*(j-1) + i + N*M*k)));
					}
					
					if ( k != 0 and !(*(visited_w + N*j + i + N*M*(k-1))) ){
						double tmp = MACGrid_->get_w(i, j, k-1) * *(counter + N*j + i + N*M*(k-1));
						*(counter + N*j + i + N*M*(k-1)) += 1;
						MACGrid_->set_w(i, j, k-1, (tmp + MACGrid_->get_w(i, j, k))/(*(counter + N*j + i + N*M*(k-1))));
					}
					
					if ( i != N-1 and !(*(visited_w + N*j + (i+1) + N*M*k)) ){
						double tmp = MACGrid_->get_w(i+1, j, k) * *(counter + N*j + (i+1) + N*M*k);
						*(counter + N*j + (i+1) + N*M*k) += 1;
						MACGrid_->set_w(i+1, j, k, (tmp + MACGrid_->get_w(i, j, k))/(*(counter + N*j + (i+1) + N*M*k)));
					}
					
					if ( j != M-1 and !(*(visited_w + N*(j+1) + i + N*M*k)) ){
						double tmp = MACGrid_->get_w(i, j+1, k) * *(counter + N*(j+1) + i + N*M*k);
						*(counter + N*(j+1) + i + N*M*k) += 1;
						MACGrid_->set_w(i, j+1, k, (tmp + MACGrid_->get_w(i, j, k))/(*(counter + N*(j+1) + i + N*M*k)));
					}
					
					if ( k != L and !(*(visited_w + N*j + i + N*M*(k+1))) ){
						double tmp = MACGrid_->get_w(i, j, k+1) * *(counter + N*j + i + N*M*(k+1));
						*(counter + N*j + i + N*M*(k+1)) += 1;
						MACGrid_->set_w(i, j, k+1, (tmp + MACGrid_->get_w(i, j, k))/(*(counter + N*j + i + N*M*(k+1))));
					}
				}
			}
		}
	}
	
	// Clear the Heap
	delete[] counter;
}


/*** APPLY EXTERNAL FORCES ***/
void FLIP::apply_forces(const double dt) {
	
	// Compute&apply external forces (gravity, vorticity confinement, ...) 
	// Apply them to the velocity field via forward euler
	// Only worry about gravity for now
	
	// Alias for MAC Grid
	auto& g = MACGrid_;
	
	// Get total number of cells on each axis
	const unsigned N = g->get_num_cells_x();
	const unsigned M = g->get_num_cells_y();
	const unsigned L = g->get_num_cells_z();
	
	// Iterate over cells & update: dv = dt*g
	for(unsigned k = 0; k < L; ++k){
		for (unsigned j = 0; j <= M; ++j){
			for (unsigned i = 0; i < N; ++i){
				g->set_v(i, j, k, g->get_v(i,j,k) - dt*gravity_mag_);
			}
		}
	}
}


void FLIP::explode(const double dt, const unsigned long step, const int x, const int y, const int z, const double r, const double value){
	
	// Apply external forces to simulate meteorite crash
	const double force = value;
	
	// Set the slope of the fall
	const double slope_x = 1;
	const double slope_y = 1;
	const double slope_z = 1;
	
	// Move the forces based on the current step
	const double sty = slope_y*step;
	const double stx = slope_x*sty;
	const double stz = slope_z*sty;
	
	// Get total number of cells on each axis
	const int nx = MACGrid_->get_num_cells_x();
	const int ny = MACGrid_->get_num_cells_y();
	const int nz = MACGrid_->get_num_cells_z();
	
	// Indices of the grid-cell at the center of the "meteorite"
	const int x_center = x + slope_x*nx/2 - stx;
	const int y_center = y + slope_y*ny/2 - sty;
	const int z_center = z + slope_z*nz/2 - stz;
	
	// Set the radius of the "meteorite"
	const double radius = r;
	
	// Iterate over all grid-velocities. If they are within the 
	// radius r from the center, update the velocties with the given 
	// force value (the force always points outwards)
	for( int k = 0; k < nz; ++k ){
		for( int j = 0; j < ny; ++j ){
			for( int i = 0; i < nx; ++i ){
				
				if (std::abs(i-x_center) <= radius and std::abs(j-y_center) <= radius and std::abs(k-z_center) <= radius){
					
					if (i-x_center == 0 or std::signbit(i-x_center)){
						MACGrid_->set_u(i, j, k, MACGrid_->get_u(i, j, k) - dt*force*(2.+(i-x_center)/radius));
					} else {
						MACGrid_->set_u(i, j, k, MACGrid_->get_u(i, j, k) + dt*force*(2.-(i-x_center)/radius));
					}
					
					if (j-y_center == 0 or std::signbit(j-y_center)){
						MACGrid_->set_v(i, j, k, MACGrid_->get_v(i, j, k) - dt*force*(2.+(j-y_center)/radius));
					} else {
						MACGrid_->set_v(i, j, k, MACGrid_->get_v(i, j, k) + dt*force*(2.-(j-y_center)/radius));
					}
					
					if (k-z_center == 0 or std::signbit(k-z_center)){
						MACGrid_->set_w(i, j, k, MACGrid_->get_w(i, j, k) - dt*force*(2.+(k-z_center)/radius));
					} else {
						MACGrid_->set_w(i, j, k, MACGrid_->get_w(i, j, k) + dt*force*(2.-(k-z_center)/radius));
					}
				}
			}
		}
	}
}


/*** BOUNDARY CONDITIONS ***/
void FLIP::apply_boundary_conditions() {
	
	// Get total number of cells on each axis
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	unsigned nz = MACGrid_->get_num_cells_z();
	
	// Enforce boundary conditions for outer (system) boundaries
	for (unsigned k = 0; k < nz; k++) {
		for (unsigned i = 0; i < nx; i++) {
			MACGrid_->set_v(i, 0, k, 0);
			MACGrid_->set_v(i, ny, k, 0);
		}
	}
	
	for (unsigned k = 0; k < nz; k++) {
		for (unsigned j = 0; j < ny; j++) {
			MACGrid_->set_u(0, j, k, 0);
			MACGrid_->set_u(nx, j, k, 0);
		}
	}
	
	for (unsigned j = 0; j < ny; j++) {
		for (unsigned i = 0; i < nx; i++) {
			MACGrid_->set_w(i, j, 0, 0);
			MACGrid_->set_w(i, j, nz, 0);
		}
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
	using solver_t = ConjugateGradient< SparseMatrix<double>, Lower|Upper, IncompleteCholesky<double> >;

	solver_t solver;
	solver.setMaxIterations(100);
	solver.compute(A_);
	VectorXd p = solver.solve(d_);

	// Copy pressures to MAC Grid
	MACGrid_->set_pressure(p);

	// Apply pressure gradients to velocity field
	//     -> see SIGGRAPH §4
	apply_pressure_gradients(dt);
}


void FLIP::compute_pressure_matrix() {
	
	// Compute matrix for pressure solve and store in A_
	// See eq. (4.19) and (4.24) in SIGGRAPH notes
	
	// Vector of triplets to construct pressure matrix with
	std::vector< Mac3d::Triplet_t > triplets;
	
	// Get total number of cells on each axis
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	unsigned nz = MACGrid_->get_num_cells_z();

	// Index of the current grid-cell [0, nx*ny*nz[ (in the matrix)
	unsigned cellidx = 0;
	
	// Iterate over all grid-cells
	for (unsigned k = 0; k < nz; ++k){	
		for (unsigned j = 0; j < ny; ++j){
			for (unsigned i = 0; i < nx; ++i, ++cellidx){
				
				// Copy diagonal entry
				auto& diag_e = MACGrid_->get_a_diag()[i + j*nx + nx*ny*k];
				triplets.push_back(diag_e);

				// Compute off-diagonal entries
				if (MACGrid_->is_fluid(i, j, k)){
					
					// x-adjacent cells
					if (i+1 < nx && MACGrid_->is_fluid(i+1, j, k)){
						
						// Compute (i+1,j,k)
						triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx+1, -1));
						
						// Use symmetry to avoid computing (i-1,j,k) separately
						triplets.push_back(Mac3d::Triplet_t(cellidx+1, cellidx, -1));
					}
					
					// y-adjacent cells
					if (j+1 < ny && MACGrid_->is_fluid(i, j+1, k)){
						
						// Compute (i,j+1,k)
						triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + nx, -1));
						
						// Use symmetry to avoid computing (i,j-1,k) separately
						triplets.push_back(Mac3d::Triplet_t(cellidx + nx, cellidx, -1));
					}
					
					// z-adjacent cells
					if (k+1 < nz && MACGrid_->is_fluid(i, j, k+1)){
						
						// Compute (i,j,k+1)
						triplets.push_back(Mac3d::Triplet_t(cellidx, cellidx + nx*ny, -1));
						
						// Use symmetry to avoid computing (i,j-1) separately
						triplets.push_back(Mac3d::Triplet_t(cellidx + nx*ny, cellidx, -1));
					}
				} // if is_fluid(i,j,k)
			}
		} // Outer for
	}
	
	// Set A_ to zero. A_ could be resized only at the start of the sim
	A_.resize(nx*ny*nz, nx*ny*nz);
	A_.setZero();
	A_.setFromTriplets(triplets.begin(), triplets.end());
}

void FLIP::compute_pressure_rhs(const double dt) {
	
	// Compute right-hand side of the pressure equations and store in d_
	// See eq. (4.19) and (4.24) in SIGGRAPH notes
	// Note: u_{solid} = 0
	
	// Get total number of cells on each axis
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	unsigned nz = MACGrid_->get_num_cells_z();

	// Alias for MAC Grid
	auto& g = MACGrid_;
	
	// Set d_ to zero. d_ could be resized only at the start of the sim
	d_.resize(nx*ny*nz);
	d_.setZero();
	
	// Index of the current grid-cell [0, nx*ny*nz[ (in the matrix)
	unsigned cellidx = 0;
	
	// Iterate over all grid-cells
	for (unsigned k = 0; k < nz; ++k) {
		for (unsigned j = 0; j < ny; ++j) {
			for (unsigned i = 0; i < nx; ++i, ++cellidx) {
				
				if (g->is_fluid(i,j,k)) {
					
					// Apply the formulas of the SIGGRAPH notes
					double d_ij = -(g->get_u(i+1,j,k) - g->get_u(i,j,k));
					d_ij -= g->get_v(i,j+1,k) - g->get_v(i,j,k);
					d_ij -= g->get_w(i,j,k+1) - g->get_w(i,j,k);
					
					// Note: u_{solid} = 0
		
					// Check each adjacent cell. If solid, alter term as in (4.24)
					// Consider cells outside of the boundary as solid
					
					// (i+1, j, k)
					if ((i < (nx-1) && g->is_solid(i+1,j,k)) || i == nx-1) {
						d_ij += g->get_u(i+1,j,k);
					}
					
					// (i-1, j, k)
					if ((i > 0 && g->is_solid(i-1,j,k)) || i == 0) {
						d_ij += g->get_u(i,j,k);
					}
		
					// (i, j+1, k)
					if ((j < (ny-1) && g->is_solid(i,j+1,k)) || j == ny-1) {
						d_ij += g->get_v(i,j+1,k);
					}
					
					// (i, j-1, k)
					if ((j > 0 && g->is_solid(i,j-1,k)) || j == 0) {
						d_ij += g->get_v(i,j,k);
					}
					
					// (i, j, k+1)
					if ((k < (nz-1) && g->is_solid(i,j,k+1)) || k == nz-1) {
						d_ij += g->get_w(i,j,k+1);
					}
					
					// (i, j, k-1)
					if ((k > 0 && g->is_solid(i,j,k-1)) || k == 0) {
						d_ij += g->get_w(i,j,k);
					}
					
					d_(cellidx) = fluid_density_ * g->get_cell_sizex() * d_ij / dt;
					
				} else { // if is_fluid(i,j,k)
					
					// Set the entry to zero if the current cell is a 
					// fluid-cell
					d_(cellidx) = 0;
				}
			}
		}
	}
}


void FLIP::apply_pressure_gradients(const double dt) {
	
	// Apply pressure gradients to velocity field
	
	// Get total number of cells on each axis
	unsigned nx = MACGrid_->get_num_cells_x();
	unsigned ny = MACGrid_->get_num_cells_y();
	unsigned nz = MACGrid_->get_num_cells_z();
	
	// Alias for MAC Grid
	auto& g = MACGrid_;
	
	// Get the length of an edge
	double dx = g->get_cell_sizex();
	
	// Iterate over all grid-cells
	for (unsigned k = 0; k < nz; ++k) {	
		for (unsigned j = 0; j < ny; ++j) {
			for (unsigned i = 0; i < nx; ++i) {
				
				// Update grid-velocities with new velocities induced by
				// pressures
				
				if (i != 0) {
					
					// get_u(i,j,k) = u_{ (i-1/2, j, k) }
					// See SIGGRAPH eq. (4.6)
					double du = (g->get_pressure(i,j,k) - g->get_pressure(i-1,j,k));
					du *= (dt/(dx*fluid_density_));
					
					g->set_u(i,j,k, g->get_u(i,j,k) - du);
				}
				
				if (j != 0) {
					
					// get_v(i,j,k) = v_{ (i, j-1/2, k) }
					// See SIGGRAPH eq. (4.7)
					double dv = (g->get_pressure(i,j,k) - g->get_pressure(i,j-1,k));
					dv *= (dt/(dx*fluid_density_));
					
					g->set_v(i,j,k, g->get_v(i,j,k) - dv);
				}
				
				if (k != 0) {
					
					// get_w(i,j,k) = w_{ (i, j, k-1/2) }
					// See SIGGRAPH eq. (4.8)
					
					double dw = (g->get_pressure(i,j,k) - g->get_pressure(i,j,k-1));
					dw *= (dt/(dx*fluid_density_));
					
					g->set_w(i,j,k, g->get_w(i,j,k) - dw);
				}
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
	
	// Get total number of cells on each axis
	int nx = MACGrid_->get_num_cells_x();
	int ny = MACGrid_->get_num_cells_y();
	int nz = MACGrid_->get_num_cells_z();
	
	// Iterate over all particles
	for(unsigned i = 0; i < num_particles_; ++i){
		
		// Store the initial positions and velocities of the particles
		Eigen::Vector3d initial_position = (particles_+i)->get_position();
		Eigen::Vector3d initial_velocity = (particles_+i)->get_velocity();
		
		// Get the index of the grid-cell containing the current 
		// particle
		auto initial_idx = MACGrid_->index_from_coord(initial_position(0), 
													  initial_position(1),
													  initial_position(2));
		
		// Initialization of variables
		Eigen::Vector3d interp_u_star;
		interp_u_star.setZero();
		Eigen::Vector3d interp_u_n1;
		interp_u_n1.setZero();
		Eigen::Vector3d u_update;
		u_update.setZero();

		// Aliases for initial_position components
		double x = initial_position[0];
		double y = initial_position[1];
		double z = initial_position[2];
		
		// With u*, v* and w* we can do the interpolation by means of 
		// interp(u*, x_p)
		// With the new u, v and w we can do the interpolation by means 
		// of interp(u_n1, x_p)
		
		// Update the horizontal velocity (trilinear interpolation)
		interp_u_star[0] = MACGrid_->get_interp_u(x,y,z,true);
		interp_u_n1[0] = MACGrid_->get_interp_u(x,y,z);
		
		// Update the vertical velocity (trilinear interpolation)
		interp_u_star[1] = MACGrid_->get_interp_v(x,y,z,true);
		interp_u_n1[1] = MACGrid_->get_interp_v(x,y,z);
		
		// Update the outgoing velocity (trilinear interpolation)
		interp_u_star[2] = MACGrid_->get_interp_w(x,y,z,true);
		interp_u_n1[2] = MACGrid_->get_interp_w(x,y,z);
	
		// Blend PIC and FLIP, use double the amount of PIC on boundary
		if (initial_idx(0) == 0 or initial_idx(0) == nx-1
		 or initial_idx(1) == 0 or initial_idx(1) == ny-1
		 or initial_idx(2) == 0 or initial_idx(2) == nz-1){
			
			u_update = initial_velocity*(1 - std::min(1., 2*alpha)) + interp_u_n1 - interp_u_star*(1 - std::min(1., 2*alpha));
		
		} else {
			
			u_update = initial_velocity*(1 - alpha) + interp_u_n1 - interp_u_star*(1 - alpha);
		
		}
		
		// Finally, update the velocities of the particles
		(particles_ + i)->set_velocity(u_update);
	}
}


double FLIP::compute_timestep( const double dt ){
	
	// New timestep that satisfies CFL condition
	double dt_new;
	
	// Particle velocity
	Eigen::Vector3d vel;
	
	// Get the maximal particle velocity components
	double u_max = 0;
	double v_max = 0;
	double w_max = 0;
	for( unsigned int n = 0; n < num_particles_; ++n ){
		vel = (particles_ + n)->get_velocity();
		if ( std::abs(vel(0)) > std::abs(u_max) ){
			u_max = vel(0);
		}
		if ( std::abs(vel(1)) > std::abs(v_max) ){
			v_max = vel(1);
		}
		if ( std::abs(vel(2)) > std::abs(w_max) ){
			w_max = vel(2);
		}
	}
	
	// Check if the fastest particles travel a distance larger than the 
	// length of an edge of a cell
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
	
	if ( w_max == 0 ){
		dt_new = dt;
	} else {
		double tmp = std::abs(MACGrid_->get_cell_sizez()/w_max);
		if ( tmp < dt_new){
			dt_new = tmp;
		}
		if ( dt_new > dt ){
			dt_new = dt;
		}
	}
	
	return dt_new;
}


void FLIP::advance_particles(const double dt, const unsigned long step) {
	
	// Iterate over all particles	
	for(unsigned n = 0; n < num_particles_; ++n){
		
		// Get current position and velocity of the particle
		Eigen::Vector3d pos_curr = (particles_ + n)->get_position();
		Eigen::Vector3d vel = (particles_ + n)->get_velocity();
		
		// Coordinates of the future location of the particle
		Eigen::Vector3d pos_next;
		
		// Euler estimate
		Eigen::Vector3d pos_half = pos_curr + 0.5*dt*vel;
		
		// Get the size of the grid (in meters)
		double size_x = MACGrid_->get_grid_size()(0);
		double size_y = MACGrid_->get_grid_size()(1);
		double size_z = MACGrid_->get_grid_size()(2);
		
		// Get the length of the edges of a grid-cell (in meters)
		double cell_sizex = MACGrid_->get_cell_sizex();
		double cell_sizey = MACGrid_->get_cell_sizey();
		double cell_sizez = MACGrid_->get_cell_sizez();
		
		// Coordinates of intermediate position, computed with Euler
		double x_half = pos_half(0);
		double y_half = pos_half(1);
		double z_half = pos_half(2);
		
		// Check if pos_half is out of the grid
		if ((x_half <= -0.5*cell_sizex) or (x_half >= size_x - 0.5*cell_sizex)
		 or (y_half <= -0.5*cell_sizey) or (y_half >= size_y - 0.5*cell_sizey)
		 or (z_half <= -0.5*cell_sizez) or (z_half >= size_z - 0.5*cell_sizez)){
			continue;
		}
		
		// RK2
		pos_next(0) = pos_curr(0) + dt*MACGrid_->get_interp_u(x_half, y_half, z_half);
		pos_next(1) = pos_curr(1) + dt*MACGrid_->get_interp_v(x_half, y_half, z_half);
		pos_next(2) = pos_curr(2) + dt*MACGrid_->get_interp_w(x_half, y_half, z_half);
		
		// Aliases for pos_next coordinates
		double x = pos_next(0);
		double y = pos_next(1);
		double z = pos_next(2);
		
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
		
		if (z <= -0.5*cell_sizez) {
			pos_next(2) = 0.;
		}
		
		if (z >= size_z - 0.5*cell_sizez) {
			pos_next(2) = size_z - cell_sizez;
		}

		// Check if the particle enters in a solid
		
		// Get the indices of the grid-cells containing the particle at 
		// the current time and in the future
		auto prev_indices = MACGrid_->index_from_coord(pos_curr(0), pos_curr(1), pos_curr(2));
		auto new_indices = MACGrid_->index_from_coord(pos_next(0), pos_next(1), pos_curr(2));
		
		// Get the length of the edges of a grid-cell (in meters)
		// [shorter names]
		double sx = MACGrid_->get_cell_sizex();
		double sy = MACGrid_->get_cell_sizey();
		double sz = MACGrid_->get_cell_sizez();
		
		// Shift a particle if it would exit the system 
		// (should not happen)
		if (MACGrid_->is_solid(new_indices(0), new_indices(1), new_indices(2))) {
			
			if (prev_indices(0)  > new_indices(0))
				pos_next(0) = (prev_indices(0) - 0.25) * sx;
			else if (prev_indices(0) < new_indices(0))
				pos_next(0) = (prev_indices(0) + 0.25) * sx;

			if (prev_indices(1) > new_indices(1))
				pos_next(1) = (prev_indices(1) - 0.25) * sy;
			else if (prev_indices(1) < new_indices(1))
				pos_next(1) = (prev_indices(1) + 0.25) * sy;
				
			if (prev_indices(2) > new_indices(2))
				pos_next(2) = (prev_indices(2) - 0.25) * sz;
			else if (prev_indices(2) < new_indices(2))
				pos_next(2) = (prev_indices(2) + 0.25) * sz;
		}
		
		// Update the position of the current particle
		(particles_ + n)->set_position(pos_next);
	}
}
