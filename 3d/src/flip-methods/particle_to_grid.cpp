/**
 * This file contains the implementation for the particle-to-grid FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/*** COMPUTE VELOCITY FIELD ***/
void FLIP::particle_to_grid() {

	// Set all grid velocities to zero and reset all fluid flags
	MACGrid_->set_velocities_to_zero();
	MACGrid_->set_weights_to_zero();
	MACGrid_->reset_fluid();

	// Sizes of the edges of a cell (in meters)
	double cell_size_x = MACGrid_->get_cell_sizex();
	double cell_size_y = MACGrid_->get_cell_sizey();
	double cell_size_z = MACGrid_->get_cell_sizez();

	// Grid indices
	Mac3d::cellIdx_t nx = MACGrid_->get_num_cells_x();
	Mac3d::cellIdx_t ny = MACGrid_->get_num_cells_y();
	Mac3d::cellIdx_t nz = MACGrid_->get_num_cells_z();

	// Threshold h and h_scaled so that it is equal to the distance
	// expressed in number of cells
	//double hx     = 2. * cell_size_x;
	//double hy     = 2. * cell_size_y;
	//double hz     = 2. * cell_size_z;
	//int hx_scaled = std::ceil(hx/cell_size_x);
	//int hy_scaled = std::ceil(hy/cell_size_y);
	//int hz_scaled = std::ceil(hz/cell_size_z);
	double h      = 2. * cell_size_x;
	int hx_scaled = std::ceil(h/cell_size_x);
	int hy_scaled = std::ceil(h/cell_size_y);
	int hz_scaled = std::ceil(h/cell_size_z);

	// Lists of flags for visited grid-velocities: 1 -> visited
	bool* visited_u = (bool*) calloc(ny*nz*(nx+1), sizeof(bool));
	bool* visited_v = (bool*) calloc(nx*nz*(ny+1), sizeof(bool));
	bool* visited_w = (bool*) calloc(nx*ny*(nz+1), sizeof(bool));

	// Position and velocity of the particleent particle
	double x_particle;
	double y_particle;
	double z_particle;
	double u_particle;
	double v_particle;
	double w_particle;

	// Coordinate of the center of a cell face
	double face_coord_x;
	double face_coord_y;
	double face_coord_z;

	// Indices of the cell containing the current particle
	int cell_idx_x;
	int cell_idx_y;
	int cell_idx_z;

	// Iterate over all particles and add weighted particle velocities
	// to grid points within a threshold h (in this case equal to double
	// the length of an edge of a cell)
	for( Particles::particleIdx_t n = 0; n < num_particles_; ++n ){

		// Get the position and velocity of the current particle
		x_particle = particles_.x[n];
		y_particle = particles_.y[n];
		z_particle = particles_.z[n];
		u_particle = particles_.u[n];
		v_particle = particles_.v[n];
		w_particle = particles_.w[n];

		// Get the indices corresponding to the cell containing the
		// current particle
		particles_.get_cell_index(n, cell_idx_x, cell_idx_y, cell_idx_z);

		// Set the cell of the current particle to a fluid-cell
		if( !(MACGrid_->is_fluid(cell_idx_x, cell_idx_y, cell_idx_z) or MACGrid_->is_solid(cell_idx_x, cell_idx_y, cell_idx_z)) ){
			
			MACGrid_->set_fluid(cell_idx_x, cell_idx_y, cell_idx_z);
		}

		// For each particle iterate only over the grid-velocities in a
		// h_scaled neighborhood
		for( int k = cell_idx_z - hz_scaled; k <= cell_idx_z + hz_scaled + 1; ++k ){
		for( int j = cell_idx_y - hy_scaled; j <= cell_idx_y + hy_scaled + 1; ++j ){
		for( int i = cell_idx_x - hx_scaled; i <= cell_idx_x + hx_scaled + 1; ++i ){

					if( i >= 0 and j >= 0 and k >= 0 ){
						
						// Coordinates of the points on the cell centers (in meters)
						// Face-offsets will be added afterwards
						face_coord_x = i * cell_size_x;
						face_coord_y = j * cell_size_y;
						face_coord_z = k * cell_size_z;

						// Left Face
						if ( ( i <= nx and j < ny and k < nz ) ){

							accumulate_vel( MACGrid_->pu_, MACGrid_->pweights_u_, u_particle, 
											x_particle, y_particle, z_particle, 
											face_coord_x - 0.5*cell_size_x, face_coord_y, face_coord_z, 
											h, i + (nx+1) * (j + ny*k) );
						}
						
						// Lower Face
						if ( ( i < nx and j <= ny and k < nz ) ){

							accumulate_vel( MACGrid_->pv_, MACGrid_->pweights_v_, v_particle, 
											x_particle, y_particle, z_particle, 
											face_coord_x, face_coord_y - 0.5*cell_size_y, face_coord_z, 
											h, i + nx * (j + (ny+1)*k) );
						}

						// Farthest Face (the closest to the origin)
						if ( ( i < nx and j < ny and k <= nz ) ){

							accumulate_vel( MACGrid_->pw_, MACGrid_->pweights_w_, w_particle, 
											x_particle, y_particle, z_particle, 
											face_coord_x, face_coord_y, face_coord_z - 0.5*cell_size_z, 
											h, i + nx * (j + ny*k) );
						}
					}
				}
			}
		}
	}

	// Normalize grid-velocities
	normalize_accumulated_vels(visited_u, visited_v, visited_w, nx, ny, nz);
	
	// Extrapolate velocities
	extrapolate_u( visited_u );
	extrapolate_v( visited_v );
	extrapolate_w( visited_w );

	// Clear the Heap
	free(visited_u);
	free(visited_v);
	free(visited_w);
}


double FLIP::compute_weight( const double r2,
							 const double h2,
							 const double h )
{
	
	double h4   = h2 * h2;
	double diff = h2 - r2;

	return (315./(64. * M_PI * h4*h4*h)) * diff*diff*diff;

}


void FLIP::accumulate_vel( double* const vel_grid,
						   double* const weights,
						   const double vel_particle,
						   const double x_particle,
						   const double y_particle,
						   const double z_particle,
						   const double face_coord_x,
						   const double face_coord_y,
						   const double face_coord_z,
						   const double h,
						   const Mac3d::globalCellIdx_t idx )
{

	// Compute (squared) L2 norm of particle position and face coordinate
	double rx = x_particle - face_coord_x;
	double ry = y_particle - face_coord_y;
	double rz = z_particle - face_coord_z;
	double r2 = rx*rx + ry*ry + rz*rz;
	double h2 = h*h;

	// If the current particle is within the given threshold h, update
	// the horizontal grid-velocity at grid_coord with the weighted
	// particle-velocity
	if( r2 <= h2 ){

		double weight = compute_weight(r2, h2, h);

		// Accumulate velocities and weights
		vel_grid[idx] += weight * vel_particle;
		weights[idx]  += weight;
	}

}


void FLIP::normalize_accumulated_vels( bool* const visited_u, 
									   bool* const visited_v, 
									   bool* const visited_w, 
									   Mac3d::cellIdx_t nx, 
									   Mac3d::cellIdx_t ny, 
									   Mac3d::cellIdx_t nz )
{

	Mac3d::globalCellIdx_t u_idx;
	Mac3d::globalCellIdx_t v_idx;
	Mac3d::globalCellIdx_t w_idx;

	double weight_u;
	double weight_v;
	double weight_w;
	
	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
			for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
				
				u_idx = i + (nx+1) * (j +  ny    * k);
				v_idx = i +  nx    * (j + (ny+1) * k);
				w_idx = i +  nx    * (j +  ny    * k);

				weight_u = MACGrid_->pweights_u_[u_idx];
				weight_v = MACGrid_->pweights_v_[v_idx];
				weight_w = MACGrid_->pweights_w_[w_idx];
				
				if( weight_u != 0 ){
					
					MACGrid_->pu_[u_idx] /= weight_u;
					visited_u[u_idx] = true;
				}

				if( weight_v != 0 ){
					
					MACGrid_->pv_[v_idx] /= weight_v;
					visited_v[v_idx] = true;
				}

				if( weight_w != 0 ){
					
					MACGrid_->pw_[w_idx] /= weight_w;
					visited_w[w_idx] = true;
				}
			}	
		}
	}

	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
			
			u_idx = nx + (nx+1) * (j + ny*k);

			weight_u = MACGrid_->pweights_u_[u_idx];

			if ( weight_u != 0 ){
				
				MACGrid_->pu_[u_idx] /= weight_u;
				visited_u[u_idx] = true;
			}
		}
	}

	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
			
			v_idx = i + nx * (ny * (k+1) + k);

			weight_v = MACGrid_->pweights_v_[v_idx];

			if ( weight_v != 0 ){
				
				MACGrid_->pv_[v_idx] /= weight_v;
				visited_v[v_idx] = true;
			}
		}
	}

	for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
		for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
			
			w_idx = i + nx * (j + ny * nz);

			weight_w = MACGrid_->pweights_w_[w_idx];

			if ( weight_w != 0 ){
				
				MACGrid_->pw_[w_idx] /= weight_w;
				visited_w[w_idx] = true;
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
	unsigned* counter = (unsigned*) calloc(M*L*(N+1), sizeof(unsigned));

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
	free(counter);
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

