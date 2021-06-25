#include <algorithm>
#include "MeshExporter.h"
#include "tsc_x86.hpp"

/**The constructor already initializes the lists point_, x_avrg_num,
 * r_avrg_num and den.
 */
MeshExporter::MeshExporter(Mac3d* Grid, Particles& particles)
	:pMacGrid_(Grid), particles_(particles), num_particles_(particles.get_num_particles()),
	dx{pMacGrid_->get_cell_sizex()},
	dy{pMacGrid_->get_cell_sizey()},
	dz{pMacGrid_->get_cell_sizez()},
	dx_r{1.0/dx},
	dy_r{1.0/dy},
	dz_r{1.0/dz},
	N{(int) pMacGrid_->get_num_cells_x()},
	M{(int) pMacGrid_->get_num_cells_y()},
	L{(int) pMacGrid_->get_num_cells_z()},
	sizex_{pMacGrid_->sizex_},
	sizey_{pMacGrid_->sizey_},
	sizez_{pMacGrid_->sizez_},
	r_avrg{0.87*dx},
	h{dx},
	h_sq{h*h},
	sdf_trivial{0.5*dx},
	plevel_set_map(nullptr, 1, 1),
	points_map(nullptr, 1, 1)
{
	
	//Grid properties
	//Initialization of the list points_
	const int s_sup = (N+2)*(M+2)*(L+2);
	points_d = new double[3*s_sup];
	//x_avrg_num_array is a matrix of N*M*L 3-vectors. it is in row-major format
	x_avrg_num_array = new double[3*N*M*L];
	plevel_set_array = new double[s_sup];
	const double p_level_init = 0.5*dx;

	// initialize the map objects
	new (&plevel_set_map) Eigen::Map<Eigen::VectorXd>(plevel_set_array, s_sup, 1);
	new (&points_map) Eigen::Map<Eigen::MatrixXd>(points_d, s_sup, 3);
	{
		const int s_sm = N + 2;
		const int s_bg = M + 2;
		for(int k = 0; k < L+2; ++k){
			for(int j = 0; j < M+2; ++j){
				for(int i = 0; i < N+2; ++i){
					int index = i + j*s_sm + k*s_sm*s_bg;
					points_d[index] = (i-1) * dx;
					points_d[index+s_sup] = (j-1) * dy;
					points_d[index+2*s_sup] = (k-1) * dz;
					plevel_set_array[index] = p_level_init;
				}
			}
		}
	}
	
	den = new double[N*M*L];

	// Create directory if it doesn't exist
	mkdir(folder_.c_str(), ACCESSPERMS);
}

void MeshExporter::level_set_easy(){
	//Grid properties
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
	
	//Compute the values of level set function
	for(int i = 0;i < N; ++i){
		for(int j = 0; j < M; ++j){
			for(int k = 0; k < L; ++k){
				int index = i + j*N + k*N*M;
				if(pMacGrid_->is_fluid(i,j,k))
					plevel_set_array[index] = -1;
				else
					plevel_set_array[index] = 1;
			}
		}
	}
}

void MeshExporter::level_set(){
	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");
	//Grid properties
	//Initialization of plevel_set_aArray, x_avrg_num and den
	std::fill(x_avrg_num_array, x_avrg_num_array+3*N*M*L, 0);
	std::fill(den, den+N*M*L, 0);

	//Compute the values of x_avrg_num and den
	tsctimer.start_timing("first_part");
	for(unsigned it_particle = 0; it_particle < num_particles_; ++it_particle){
		//const Particle& particle = *(particles_ + it_particle);
		double particle_x = particles_.x[it_particle];
		double particle_y = particles_.y[it_particle];
		double particle_z = particles_.z[it_particle];

		// inlined coords_to_index
		assert(
			(particle_x < sizex_ - 0.5*dx
			 && particle_y < sizey_ - 0.5*dy
			 && particle_z < sizez_ - 0.5*dz
			 && particle_x > -0.5*dx
			 && particle_y > -0.5*dy
			 && particle_z > -0.5*dz
			 ) && "Attention: out of the grid!");

		Mac3d::cellIdx_t init_cell_x;
		Mac3d::cellIdx_t init_cell_y;
		Mac3d::cellIdx_t init_cell_z;
		particles_.get_cell_index(it_particle, init_cell_x, init_cell_y, init_cell_z);

		const int k_max = std::min((int) L, init_cell_z + 2);
		const int j_max = std::min((int) M, init_cell_y + 2);
		const int i_max = std::min((int) N, init_cell_x + 2);
		for(int k = std::max(0, init_cell_z - 2); k < k_max; ++k){
			for(int j = std::max(0, init_cell_y - 2); j < j_max; ++j){
				for(int i = std::max(0, init_cell_x - 2); i < i_max; ++i){
					const int index = i + j*N + k*N*M;
					const double dist_x = i*dx - particle_x;
					const double dist_y = j*dy - particle_y;
					const double dist_z = k*dz - particle_z;
					const double s_sq = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
					if (s_sq < h_sq){
						const double s_sq_inv = h_sq - s_sq;
						const double W_surf = s_sq_inv*s_sq_inv*s_sq_inv;
						x_avrg_num_array[3*index]   += W_surf*particle_x;
						x_avrg_num_array[3*index+1] += W_surf*particle_y;
						x_avrg_num_array[3*index+2] += W_surf*particle_z;
						den[index] += W_surf;
					}
				}
			}
		}
	}


	tsctimer.stop_timing("first_part", true, "");
	//Compute the values of level set function
	for(int k = 0; k < L; ++k){
		for(int j = 0; j < M; ++j){
			for(int i = 0; i < N; ++i){
				int index = (i+1) + (j+1)*(N+2) + (k+1)*(N+2)*(M+2);
					int index2 = i + j*N + k*N*M;
					const double denominator_inv = 1.0/ *(den+index2);
					const double x_avrg_x = x_avrg_num_array[index2*3  ] * denominator_inv;
					const double x_avrg_y = x_avrg_num_array[index2*3+1] * denominator_inv;
					const double x_avrg_z = x_avrg_num_array[index2*3+2] * denominator_inv;
					
					if(*(den+index2) != 0) {
						// eqn (18)
						plevel_set_array[index] = std::sqrt(
							std::pow(i*dx - x_avrg_x, 2)
							+ std::pow(j*dy - x_avrg_y, 2)
							+ std::pow(k*dz - x_avrg_z, 2)
						) - r_avrg;
					}
					else{
						if(pMacGrid_->is_fluid(i,j,k))
							plevel_set_array[index] = -sdf_trivial;
						else
							plevel_set_array[index] = sdf_trivial;
					}
			}
		}
	}
}


void MeshExporter::compute_mesh() {
	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");

	// Perform calculation of mesh
	tsctimer.start_timing("level_set");
	level_set();
	tsctimer.stop_timing("level_set", true, "");
	tsctimer.start_timing("marching_cubes");
	igl::copyleft::marching_cubes(plevel_set_map, points_map, N+2, M+2, L+2, vertices_, faces_);
	tsctimer.stop_timing("marching_cubes", true, "");
}

void MeshExporter::get_mesh(Eigen::MatrixXd &vertices, Eigen::MatrixXi &faces) const {
	vertices = vertices_;
	faces = faces_;
}

void MeshExporter::export_mesh() {
	const unsigned pad_width = 6;

	// Assemble file name
	std::stringstream filename;
	filename << folder_ << file_prefix_ << std::setfill('0') << std::setw(pad_width) << num_exported_;
	filename << ".obj";

	// Export mesh to OBJ file
	igl::writeOBJ(filename.str(), vertices_, faces_);
	std::cout << "Exported " << filename.str() << std::endl;

	num_exported_++;
}
