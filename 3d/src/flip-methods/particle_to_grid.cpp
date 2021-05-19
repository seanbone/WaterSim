/**
 * This file contains the implementation for the particle-to-grid FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"


void FLIP::particle_to_grid() {

	// Grid dimensions
	const Mac3d::cellIdx_t nx = MACGrid_->N_;
	const Mac3d::cellIdx_t ny = MACGrid_->M_;
	const Mac3d::cellIdx_t nz = MACGrid_->L_;

	// Sizes of the edges of a cell (in meters)
	const double cell_size_x = MACGrid_->cell_sizex_;
	const double cell_size_y = MACGrid_->cell_sizey_;
	const double cell_size_z = MACGrid_->cell_sizez_;

	const double cell_size_x_half = 0.5 * cell_size_x;
	const double cell_size_y_half = 0.5 * cell_size_y;
	const double cell_size_z_half = 0.5 * cell_size_z;
	
	// Threshold h and coefficient for weight computation
	const double h     = 2. * cell_size_x;
	const double h2    = h*h;
	const double h4    = h2*h2;
	const double coeff = 315./(64. * M_PI * h4*h4*h);

	// Threshold h expressed in number of cells
	const int hx_scaled = std::ceil(h/cell_size_x);
	const int hy_scaled = std::ceil(h/cell_size_y);
	const int hz_scaled = std::ceil(h/cell_size_z);

	// Lists of flags for visited grid-velocities: 1 -> visited
	bool* visited_u = (bool*) calloc(ny*nz*(nx+1), sizeof(bool));
	bool* visited_v = (bool*) calloc(nx*nz*(ny+1), sizeof(bool));
	bool* visited_w = (bool*) calloc(nx*ny*(nz+1), sizeof(bool));

	// Set all grid velocities to zero and reset all fluid flags
	MACGrid_->set_velocities_to_zero();
	MACGrid_->set_weights_to_zero();
	MACGrid_->reset_fluid();

	// Counters of nearby visited faces (to average the neighboring velocities 
	// during extrapolation)
	short u_counter;
	short v_counter;
	short w_counter;

	// Nearby faces velocities
	double u_left;
	double u_right;
	double u_down;
	double u_up;
	double u_back;
	double u_front;

	double v_left;
	double v_right;
	double v_down;
	double v_up;
	double v_back;
	double v_front;

	double w_left;
	double w_right;
	double w_down;
	double w_up;
	double w_back;
	double w_front;

	// Position and velocity of the current particle
	double x_particle;
	double y_particle;
	double z_particle;
	double u_particle;
	double v_particle;
	double w_particle;

	// Distance from the particle and the cell centers
	double rx;
	double ry;
	double rz;

	// Difference between threshold and distance of the particle from the face 
	// center
	double x_diff;
	double y_diff;
	double z_diff;

	// Temporary variables to store weights and remove aliasing
	double u_weight;
	double v_weight;
	double w_weight;

	// Placeholders for the global indices of u, v and w
	Mac3d::globalCellIdx_t u_idx;
	Mac3d::globalCellIdx_t v_idx;
	Mac3d::globalCellIdx_t w_idx;

	// Indices of the cell containing the current particle
	int cell_idx_x;
	int cell_idx_y;
	int cell_idx_z;

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
		if( !(MACGrid_->pfluid_[cell_idx_x + nx*cell_idx_y + nx*ny*cell_idx_z] or MACGrid_->psolid_[cell_idx_x + nx*cell_idx_y + nx*ny*cell_idx_z]) ){
			
			MACGrid_->pfluid_[cell_idx_x + nx*cell_idx_y + nx*ny*cell_idx_z] = true;
		}

		// For each particle iterate only over the grid-velocities in a
		// h_scaled neighborhood
		if( cell_idx_x >= hx_scaled          and cell_idx_y >= hy_scaled          and cell_idx_z >= hz_scaled          and 
			cell_idx_x <  nx - hx_scaled - 1 and cell_idx_y <  ny - hy_scaled - 1 and cell_idx_z <  nz - hz_scaled - 1 )
		{

			for( Mac3d::cellIdx_t k = cell_idx_z - hz_scaled; k <= cell_idx_z + hz_scaled + 1; ++k ){
			for( Mac3d::cellIdx_t j = cell_idx_y - hy_scaled; j <= cell_idx_y + hy_scaled + 1; ++j ){
			for( Mac3d::cellIdx_t i = cell_idx_x - hx_scaled; i <= cell_idx_x + hx_scaled + 1; ++i ){

						rx = x_particle - i * cell_size_x;
						ry = y_particle - j * cell_size_y;
						rz = z_particle - k * cell_size_z;

						x_diff = h2 - ry*ry - rz*rz - (rx + cell_size_x_half)*(rx + cell_size_x_half);
						y_diff = h2 - rx*rx - rz*rz - (ry + cell_size_y_half)*(ry + cell_size_y_half);
						z_diff = h2 - rx*rx - ry*ry - (rz + cell_size_z_half)*(rz + cell_size_z_half);
						
						// If the current particle is within the given threshold h, update
						// the horizontal grid-velocity at grid_coord with the weighted
						// particle-velocity

						// Left Face
						if( x_diff >= 0. ){

							u_weight = coeff * x_diff * x_diff * x_diff;

							u_idx = i + (nx+1) * (j + ny*k);

							MACGrid_->pu_[u_idx]         += u_weight * u_particle;
							MACGrid_->pweights_u_[u_idx] += u_weight;
						}
						
						// Lower Face
						if( y_diff >= 0. ){

							v_weight = coeff * y_diff * y_diff * y_diff;

							v_idx = i + nx * (j + (ny+1) * k);

							MACGrid_->pv_[v_idx]         += v_weight * v_particle;
							MACGrid_->pweights_v_[v_idx] += v_weight;
						}

						// Farthest Face (the closest to the origin)
						if( z_diff >= 0. ){

							w_weight = coeff * z_diff * z_diff * z_diff;

							w_idx = i + nx * (j + ny*k);

							MACGrid_->pw_[w_idx]         += w_weight * w_particle;
							MACGrid_->pweights_w_[w_idx] += w_weight;
						}
					}
				}
			}

		}
		else{

			for( int k = cell_idx_z - hz_scaled; k <= cell_idx_z + hz_scaled + 1; ++k ){
			for( int j = cell_idx_y - hy_scaled; j <= cell_idx_y + hy_scaled + 1; ++j ){
			for( int i = cell_idx_x - hx_scaled; i <= cell_idx_x + hx_scaled + 1; ++i ){

						if( i >= 0 and j >= 0 and k >= 0 ){
							
							rx = x_particle - i * cell_size_x;
							ry = y_particle - j * cell_size_y;
							rz = z_particle - k * cell_size_z;

							// If the current particle is within the given threshold h, update
							// the horizontal grid-velocity at grid_coord with the weighted
							// particle-velocity

							// Left Face
							if ( ( i <= nx and j < ny and k < nz ) ){

								x_diff = h2 - ry*ry - rz*rz - (rx + cell_size_x_half)*(rx + cell_size_x_half);

								if( x_diff >= 0. ){

									u_weight = coeff * x_diff * x_diff * x_diff;

									u_idx = i + (nx+1) * (j + ny*k);

									MACGrid_->pu_[u_idx]         += u_weight * u_particle;
									MACGrid_->pweights_u_[u_idx] += u_weight;
								}
							}
							
							// Lower Face
							if ( ( i < nx and j <= ny and k < nz ) ){
								
								y_diff = h2 - rx*rx - rz*rz - (ry + cell_size_y_half)*(ry + cell_size_y_half);

								if( y_diff >= 0. ){

									v_weight = coeff * y_diff * y_diff * y_diff;

									v_idx = i + nx * (j + (ny+1)*k);

									MACGrid_->pv_[v_idx]         += v_weight * v_particle;
									MACGrid_->pweights_v_[v_idx] += v_weight;
								}
							}

							// Farthest Face (the closest to the origin)
							if ( ( i < nx and j < ny and k <= nz ) ){

								z_diff = h2 - rx*rx - ry*ry - (rz + cell_size_z_half)*(rz + cell_size_z_half);

								if( z_diff >= 0. ){

									w_weight = coeff * z_diff * z_diff * z_diff;

									w_idx = i + nx * (j + ny*k);

									MACGrid_->pw_[w_idx]         += w_weight * w_particle;
									MACGrid_->pweights_w_[w_idx] += w_weight;
								}
							}
						}
					}
				}
			}

		}
	}
	
	// Normalize the accumulated velocities with the accumulated weights and 
	// set the flags for the visited faces
	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
			for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
				
				u_idx = i + (nx+1) * (j +  ny    * k);
				v_idx = i +  nx    * (j + (ny+1) * k);
				w_idx = i +  nx    * (j +  ny    * k);

				u_weight = MACGrid_->pweights_u_[u_idx];
				v_weight = MACGrid_->pweights_v_[v_idx];
				w_weight = MACGrid_->pweights_w_[w_idx];
				
				if( u_weight != 0. ){
					
					MACGrid_->pu_[u_idx] /= u_weight;
					visited_u[u_idx] = true;
				}

				if( v_weight != 0. ){
					
					MACGrid_->pv_[v_idx] /= v_weight;
					visited_v[v_idx] = true;
				}

				if( w_weight != 0. ){
					
					MACGrid_->pw_[w_idx] /= w_weight;
					visited_w[w_idx] = true;
				}
			}	
		}
	}

	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
			
			u_idx = nx + (nx+1) * (j + ny*k);

			u_weight = MACGrid_->pweights_u_[u_idx];

			if ( u_weight != 0. ){
				
				MACGrid_->pu_[u_idx] /= u_weight;
				visited_u[u_idx] = true;
			}
		}
	}

	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
			
			v_idx = i + nx * (ny * (k+1) + k);

			v_weight = MACGrid_->pweights_v_[v_idx];

			if ( v_weight != 0. ){
				
				MACGrid_->pv_[v_idx] /= v_weight;
				visited_v[v_idx] = true;
			}
		}
	}

	for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
		for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
			
			w_idx = i + nx * (j + ny * nz);

			w_weight = MACGrid_->pweights_w_[w_idx];

			if ( w_weight != 0. ){
				
				MACGrid_->pw_[w_idx] /= w_weight;
				visited_w[w_idx] = true;
			}
		}
	}

	// Iterate over all horizontal grid-velocities and extrapolate into
	// the air cells (not visited) the average velocities of the
	// neighboring fluid/visited cells
	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
			for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){
				
				u_idx = i + (nx+1) * (j +  ny    * k);
				v_idx = i +  nx    * (j + (ny+1) * k);
				w_idx = i +  nx    * (j +  ny    * k);

				if( !visited_u[u_idx] ){

					u_counter = 0;
					
					if( visited_u[u_idx-1        ] and i > 0    ) { u_left  = MACGrid_->pu_[u_idx-1        ]; ++u_counter; } else { u_left  = 0.; } // Left
					if( visited_u[u_idx+1        ]              ) { u_right = MACGrid_->pu_[u_idx+1        ]; ++u_counter; } else { u_right = 0.; } // Right
					if( visited_u[u_idx-(nx+1)   ] and j > 0    ) { u_down  = MACGrid_->pu_[u_idx-(nx+1)   ]; ++u_counter; } else { u_down  = 0.; } // Down
					if( visited_u[u_idx+(nx+1)   ] and j < ny-1 ) { u_up    = MACGrid_->pu_[u_idx+(nx+1)   ]; ++u_counter; } else { u_up    = 0.; } // Up
					if( visited_u[u_idx-(nx+1)*ny] and k > 0    ) { u_back  = MACGrid_->pu_[u_idx-(nx+1)*ny]; ++u_counter; } else { u_back  = 0.; } // Back
					if( visited_u[u_idx+(nx+1)*ny] and k < nz-1 ) { u_front = MACGrid_->pu_[u_idx+(nx+1)*ny]; ++u_counter; } else { u_front = 0.; } // Front

					if(u_counter != 0) MACGrid_->pu_[u_idx] = (u_left + u_right + u_down + u_up + u_back + u_front) / u_counter;
				}

				if( !visited_v[v_idx] ){

					v_counter = 0;
					
					if( visited_v[v_idx-1        ] and i > 0    ) { v_left  = MACGrid_->pv_[v_idx-1        ]; ++v_counter; } else { v_left  = 0.; } // Left
					if( visited_v[v_idx+1        ] and i < nx-1 ) { v_right = MACGrid_->pv_[v_idx+1        ]; ++v_counter; } else { v_right = 0.; } // Right
					if( visited_v[v_idx-nx       ] and j > 0    ) { v_down  = MACGrid_->pv_[v_idx-nx       ]; ++v_counter; } else { v_down  = 0.; } // Down
					if( visited_v[v_idx+nx       ]              ) { v_up    = MACGrid_->pv_[v_idx+nx       ]; ++v_counter; } else { v_up    = 0.; } // Up
					if( visited_v[v_idx-nx*(ny+1)] and k > 0    ) { v_back  = MACGrid_->pv_[v_idx-nx*(ny+1)]; ++v_counter; } else { v_back  = 0.; } // Back
					if( visited_v[v_idx+nx*(ny+1)] and k < nz-1 ) { v_front = MACGrid_->pv_[v_idx+nx*(ny+1)]; ++v_counter; } else { v_front = 0.; } // Front

					if(v_counter != 0) MACGrid_->pv_[v_idx] = (v_left + v_right + v_down + v_up + v_back + v_front) / v_counter;
				}

				if( !visited_w[w_idx] ){

					w_counter = 0;
					
					if( visited_w[w_idx-1    ] and i > 0    ) { w_left  = MACGrid_->pw_[w_idx-1    ]; ++w_counter; } else { w_left  = 0.; } // Left
					if( visited_w[w_idx+1    ] and i < nx-1 ) { w_right = MACGrid_->pw_[w_idx+1    ]; ++w_counter; } else { w_right = 0.; } // Right
					if( visited_w[w_idx-nx   ] and j > 0    ) { w_down  = MACGrid_->pw_[w_idx-nx   ]; ++w_counter; } else { w_down  = 0.; } // Down
					if( visited_w[w_idx+nx   ] and j < ny-1 ) { w_up    = MACGrid_->pw_[w_idx+nx   ]; ++w_counter; } else { w_up    = 0.; } // Up
					if( visited_w[w_idx-nx*ny] and k > 0    ) { w_back  = MACGrid_->pw_[w_idx-nx*ny]; ++w_counter; } else { w_back  = 0.; } // Back
					if( visited_w[w_idx+nx*ny]              ) { w_front = MACGrid_->pw_[w_idx+nx*ny]; ++w_counter; } else { w_front = 0.; } // Front

					if(w_counter != 0) MACGrid_->pw_[w_idx] = (w_left + w_right + w_down + w_up + w_back + w_front) / w_counter;
				}
			}
		}
	}

	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){

			u_idx = nx + (nx+1) * (j + ny * k);

			if( !visited_u[u_idx] ){

				u_counter = 0;
				
				if( visited_u[u_idx-1        ]              ) { u_left  = MACGrid_->pu_[u_idx-1        ]; ++u_counter; } else { u_left  = 0.; } // Left
				if( visited_u[u_idx-(nx+1)   ] and j > 0    ) { u_down  = MACGrid_->pu_[u_idx-(nx+1)   ]; ++u_counter; } else { u_down  = 0.; } // Down
				if( visited_u[u_idx+(nx+1)   ] and j < ny-1 ) { u_up    = MACGrid_->pu_[u_idx+(nx+1)   ]; ++u_counter; } else { u_up    = 0.; } // Up
				if( visited_u[u_idx-(nx+1)*ny] and k > 0    ) { u_back  = MACGrid_->pu_[u_idx-(nx+1)*ny]; ++u_counter; } else { u_back  = 0.; } // Back
				if( visited_u[u_idx+(nx+1)*ny] and k < nz-1 ) { u_front = MACGrid_->pu_[u_idx+(nx+1)*ny]; ++u_counter; } else { u_front = 0.; } // Front

				if(u_counter != 0) MACGrid_->pu_[u_idx] = (u_left + u_down + u_up + u_back + u_front) / u_counter;
			}
		}
	}

	for( Mac3d::cellIdx_t k = 0; k < nz; ++k ){
		for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){

			v_idx = i + nx * (ny * (k+1) + k);

			if( !visited_v[v_idx] ){

				v_counter = 0;
				
				if( visited_v[v_idx-1        ] and i > 0    ) { v_left  = MACGrid_->pv_[v_idx-1        ]; ++v_counter; } else { v_left  = 0.; } // Left
				if( visited_v[v_idx+1        ] and i < nx-1 ) { v_right = MACGrid_->pv_[v_idx+1        ]; ++v_counter; } else { v_right = 0.; } // Right
				if( visited_v[v_idx-nx       ]              ) { v_down  = MACGrid_->pv_[v_idx-nx       ]; ++v_counter; } else { v_down  = 0.; } // Down
				if( visited_v[v_idx-nx*(ny+1)] and k > 0    ) { v_back  = MACGrid_->pv_[v_idx-nx*(ny+1)]; ++v_counter; } else { v_back  = 0.; } // Back
				if( visited_v[v_idx+nx*(ny+1)] and k < nz-1 ) { v_front = MACGrid_->pv_[v_idx+nx*(ny+1)]; ++v_counter; } else { v_front = 0.; } // Front

				if(v_counter != 0) MACGrid_->pv_[v_idx] = (v_left + v_right + v_down + v_back + v_front) / v_counter;
			}
		}
	}

	for( Mac3d::cellIdx_t j = 0; j < ny; ++j ){
		for( Mac3d::cellIdx_t i = 0; i < nx; ++i ){

			w_idx = i + nx * (j + ny * nz);

			if( !visited_w[w_idx] ){

				w_counter = 0;
				
				if( visited_w[w_idx-1    ] and i > 0    ) { w_left  = MACGrid_->pw_[w_idx-1    ]; ++w_counter; } else { w_left  = 0.; } // Left
				if( visited_w[w_idx+1    ] and i < nx-1 ) { w_right = MACGrid_->pw_[w_idx+1    ]; ++w_counter; } else { w_right = 0.; } // Righ
				if( visited_w[w_idx-nx   ] and j > 0    ) { w_down  = MACGrid_->pw_[w_idx-nx   ]; ++w_counter; } else { w_down  = 0.; } // Down
				if( visited_w[w_idx+nx   ] and j < ny-1 ) { w_up    = MACGrid_->pw_[w_idx+nx   ]; ++w_counter; } else { w_up    = 0.; } // Up
				if( visited_w[w_idx-nx*ny]              ) { w_back  = MACGrid_->pw_[w_idx-nx*ny]; ++w_counter; } else { w_back  = 0.; } // Back

				if(w_counter != 0) MACGrid_->pw_[w_idx] = (w_left + w_right + w_down + w_up + w_back) / w_counter;
			}
		}
	}

	free(visited_u);
	free(visited_v);
	free(visited_w);
}