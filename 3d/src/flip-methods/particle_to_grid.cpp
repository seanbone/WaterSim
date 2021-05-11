/**
 * This file contains the implementation for the particle-to-grid FLIP substep and all related methods.
 */

#include "FLIP.h"
#include "tsc_x86.hpp"

/*** COMPUTE VELOCITY FIELD ***/
void FLIP::particle_to_grid() {

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
        pos = (particlesOLD_ + n)->get_position();
        vel = (particlesOLD_ + n)->get_velocity();

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

