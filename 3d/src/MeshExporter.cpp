#include <algorithm>
#include "MeshExporter.h"
#include "tsc_x86.hpp"

/**The constructor already initializes the lists point_, x_avrg_num,
 * r_avrg_num and den.
 */
MeshExporter::MeshExporter(Mac3d* Grid, Particle* particles, const int n)
	:pMacGrid_(Grid), pparticles_(particles), num_particles_(n),
	dx{pMacGrid_->get_cell_sizex()},
	dy{pMacGrid_->get_cell_sizey()},
	dz{pMacGrid_->get_cell_sizez()},
	sizex_{pMacGrid_->sizex_},
	sizey_{pMacGrid_->sizey_},
	sizez_{pMacGrid_->sizez_},
	N{pMacGrid_->get_num_cells_x()},
	M{pMacGrid_->get_num_cells_y()},
	L{pMacGrid_->get_num_cells_z()},
	h{dx},
	h_sq_r{1.0/h},
	weight_factor{0.87*dx}
{
	
	//Grid properties
	//Initialization of the list points_
	points_.resize((N+2)*(M+2)*(L+2), 3);
	points_.setZero();
	
	for(int k = -1; k < L+1; ++k){
		for(int j = -1; j < M+1; ++j){
			for(int i = -1; i < N+1; ++i){
				int index = (i+1) + (j+1)*(N+2) + (k+1)*(N+2)*(M+2);
				Eigen::RowVector3d temp = Eigen::RowVector3d(i*dx, j*dy, k*dz);
				points_.row(index) = temp;
			}
		}
	}
	
	//Initialization of x_avrg_num, r_avrg_num and den
	x_avrg_num.resize(N*M*L);
	std::fill(x_avrg_num.begin(), x_avrg_num.end(), Eigen::Vector3d::Zero());
	
	r_avrg_num = new double[N*M*L];
	std::fill(r_avrg_num, r_avrg_num+N*M*L, 0);
	
	den = new double[N*M*L];
	std::fill(den, den+N*M*L, 0);

	// Create directory if it doesn't exist
	mkdir(folder_.c_str(), ACCESSPERMS);
}

void MeshExporter::level_set_easy(){
	//Grid properties
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
	
	//Compute the values of level set function
	plevel_set_.resize(N*M*L);
	for(int i = 0;i < N; ++i){
		for(int j = 0; j < M; ++j){
			for(int k = 0; k < L; ++k){
				int index = i + j*N + k*N*M;
				if(pMacGrid_->is_fluid(i,j,k))
					plevel_set_(index) = -1;
				else
					plevel_set_(index) = 1;
			}
		}
	}
}

void MeshExporter::level_set(){
	//Grid properties
	//Initialization of plevel_set_, x_avrg_num, r_avrg_num and den
	plevel_set_.resize((N+2)*(M+2)*(L+2));
	plevel_set_.setZero();
	std::fill(x_avrg_num.begin(), x_avrg_num.end(), Eigen::Vector3d::Zero());
	std::fill(r_avrg_num, r_avrg_num+N*M*L, 0);
	std::fill(den, den+N*M*L, 0);

	//Compute the values of x_avrg_num, r_avrg_num and den
	for(unsigned it_particle = 0; it_particle < num_particles_; ++it_particle){
		const Particle& particle = *(pparticles_+it_particle);
		const Eigen::Vector3d particle_pos = (pparticles_+it_particle)->get_position();
		assert(
			(particle.x_ < sizex_ - 0.5*dx
			 && particle.y_ < sizey_ - 0.5*dy
			 && particle.z_ < sizez_ - 0.5*dz
			 && particle.x_ > -0.5*dx
			 && particle.y_ > -0.5*dy
			 && particle.z_ > -0.5*dz
			 ) && "Attention: out of the grid!");
		const int init_cell_x = int(particle.x_/dx + 0.5);
		const int init_cell_y = int(particle.y_/dy + 0.5);
		const int init_cell_z = int(particle.z_/dz + 0.5);

		int k_max = std::min((int) L, init_cell_x + 2);
		int j_max = std::min((int) M, init_cell_y + 2);
		int i_max = std::min((int) N, init_cell_z + 2);
		for(int k = std::max(0, init_cell_x - 2); k <= k_max; ++k){
			for(int j = std::max(0, init_cell_y - 2); j <= j_max; ++j){
				for(int i = std::max(0, init_cell_z - 2); i <= i_max; ++i){
					const int index = i + j*N + k*N*M;
					const double dist_x = i*dx - particle.x_;
					const double dist_y = j*dy - particle.y_;
					const double dist_z = k*dz - particle.z_;
					const double s_sq = (dist_x*dist_x + dist_y*dist_y + dist_z*dist_z) * h_sq_r;
					if (s_sq < 1){
						double s_sq_inv = 1 - s_sq;
						double W_surf = s_sq_inv*s_sq_inv*s_sq_inv;
						x_avrg_num[index] += W_surf*particle_pos;
						*(r_avrg_num + index) += W_surf*weight_factor;
						*(den + index) += W_surf;
					}
				}
			}
		}
	}
	
	//Compute the values of level set function
	for(int k = -1; k < L+1; ++k){
		for(int j = -1; j < M+1; ++j){
			for(int i = -1; i < N+1; ++i){
				int index = (i+1) + (j+1)*(N+2) + (k+1)*(N+2)*(M+2);
				
				if(i == -1 || i == N || j == -1 || j == M || k == -1 || k == L){
					plevel_set_[index] = 0.5*dx;
				}
				
				else{
					int index2 = i + j*N + k*N*M;
					double temp = *(den+index2);
					Eigen::Vector3d x_avrg = x_avrg_num[index2]/temp;
					double r_avrg = 0.87*dx;
					Eigen::Vector3d cell_pos = Eigen::Vector3d(i*dx, j*dy, k*dz);
					
					if(*(den+index2) != 0)
						// eqn (18)
						plevel_set_[index] = (cell_pos - x_avrg).norm() - r_avrg;
					else{
						if(pMacGrid_->is_fluid(i,j,k))
							plevel_set_[index] = -1;
						else
							plevel_set_[index] = 1;
					}
				}
			}
		}
	}
}


void MeshExporter::compute_mesh() {
	unsigned nx, ny, nz;
	nx = pMacGrid_->get_num_cells_x();
	ny = pMacGrid_->get_num_cells_y();
	nz = pMacGrid_->get_num_cells_z();
	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");

	// Perform calculation of mesh
	tsctimer.start_timing("level_set");
	level_set();
	tsctimer.stop_timing("level_set", true, "");
	tsctimer.start_timing("marching_cubes");
	igl::copyleft::marching_cubes(plevel_set_, points_, nx+2, ny+2, nz+2, vertices_, faces_);
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
