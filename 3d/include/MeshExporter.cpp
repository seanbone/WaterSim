#include "MeshExporter.h"

MeshExporter::MeshExporter(Mac3d* Grid, Particle* particles, const int n)
	:pMacGrid_(Grid), pparticles_(particles), num_particles_(n){
	
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
	points_.resize((N+2)*(M+2)*(L+2), 3);
	points_.setZero();
	double dx = pMacGrid_->get_cell_sizex();
	double dy = pMacGrid_->get_cell_sizey();
	double dz = pMacGrid_->get_cell_sizez();
	
	for(int k = -1; k < L+1; ++k){
		for(int j = -1; j < M+1; ++j){
			for(int i = -1; i < N+1; ++i){
				int index = (i+1) + (j+1)*(N+2) + (k+1)*(N+2)*(M+2);
				Eigen::RowVector3d temp = Eigen::RowVector3d(i*dx, j*dy, k*dz);
				points_.row(index) = temp;
			}
		}
	}
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
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
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
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
	plevel_set_.resize((N+2)*(M+2)*(L+2));
	plevel_set_.setZero();
	std::fill(x_avrg_num.begin(), x_avrg_num.end(), Eigen::Vector3d::Zero());
	std::fill(r_avrg_num, r_avrg_num+N*M*L, 0);
	std::fill(den, den+N*M*L, 0);
	double dx = pMacGrid_->get_cell_sizex();
	double dy = pMacGrid_->get_cell_sizey();
	double dz = pMacGrid_->get_cell_sizez();
	double h = dx;
	
	for(unsigned it_particle = 0; it_particle < num_particles_; ++it_particle){
		Eigen::Vector3d particle_pos = (pparticles_+it_particle)->get_position();
		Eigen::Vector3d init_cell = pMacGrid_->index_from_coord(particle_pos[0], particle_pos[1], particle_pos[2]);
		for(int k = init_cell(2) - 2; k <= init_cell(2) + 2; ++k){
			for(int j = init_cell(1) - 2; j <= init_cell(1) + 2; ++j){
				for(int i = init_cell(0) - 2; i <= init_cell(0) + 2; ++i){
					if (k < 0 or k >= L or j < 0 or j >= M or i < 0 or i >= N ){
						continue;
					}
					
					Eigen::Vector3d cell = Eigen::Vector3d(i*dx,j*dy,k*dz);
					int index = i + j*N + k*N*M;
					double temp = (cell - particle_pos).norm()/h;
					if (temp < 1){
						//double W_surf = (1-temp)*(1-temp)*(1-temp);
						double W_surf = (1-temp*temp)*(1-temp*temp)*(1-temp*temp);
						//double W_surf = 1-temp*temp*temp;
						x_avrg_num[index] += W_surf*particle_pos;
						*(r_avrg_num + index) += W_surf*0.87*dx;
						*(den + index) += W_surf;
					}
				}
			}
		}
	}
	
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
					//~ double r_avrg = *(r_avrg_num+index)/temp;
					double r_avrg = 0.87*dx;
					Eigen::Vector3d cell_pos = Eigen::Vector3d(i*dx, j*dy, k*dz);
					if(*(den+index2) != 0)
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


void MeshExporter::export_mesh() {
	const unsigned pad_width = 6;

	unsigned nx, ny, nz;
	nx = pMacGrid_->get_num_cells_x();
	ny = pMacGrid_->get_num_cells_y();
	nz = pMacGrid_->get_num_cells_z();

	// Assemble file name
	std::stringstream filename;
	filename << folder_ << file_prefix_ << std::setfill('0') << std::setw(pad_width) << num_exported_;
	filename << ".obj";

	// Perform calculation of mesh
	level_set();
	igl::copyleft::marching_cubes(plevel_set_, points_, nx+2, ny+2, nz+2, vertices_, faces_);

	// Export mesh to OBJ file
	igl::writeOBJ(filename.str(), vertices_, faces_);
	std::cout << "Exported " << filename.str() << std::endl;

	num_exported_++;
}
