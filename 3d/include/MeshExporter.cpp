#include "MeshExporter.h"

MeshExporter::MeshExporter(Mac3d* Grid, Particle* particles, const int n)
	:pMacGrid_(Grid), pparticles_(particles), num_particles_(n){
	
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
	points_.resize((N)*(M)*(L), 3);
	double dx = pMacGrid_->get_cell_sizex();
	double dy = pMacGrid_->get_cell_sizey();
	double dz = pMacGrid_->get_cell_sizez();
	
	//~ for(int k = 0; k < L + 2; ++k){
		//~ for(int j = 0; j < M + 2; ++j){
			//~ for(int i = 0; i < N + 2; ++i){
				//~ int index = i + j*(N+2) + k*(N+2)*(M+2);
				//~ Eigen::RowVector3d temp = Eigen::RowVector3d(i*dx, j*dy, k*dz);
				//~ points_.row(index) = temp;
			//~ }
		//~ }
	//~ }
	
	for(int k = 0; k < L; ++k){
		for(int j = 0; j < M; ++j){
			for(int i = 0; i < N; ++i){
				int index = i + j*(N) + k*(N)*(M);
				Eigen::RowVector3d temp = Eigen::RowVector3d(i*dx, j*dy, k*dz);
				points_.row(index) = temp;
			}
		}
	}
	x_avrg_num.resize(N*M*L);
	r_avrg_num = new double[N*M*L];
	den = new double[N*M*L];

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

void MeshExporter::level_set_open(){
	int N = pMacGrid_->get_num_cells_x();
	int M = pMacGrid_->get_num_cells_y();
	int L = pMacGrid_->get_num_cells_z();
	plevel_set_.resize((N+2)*(M+2)*(L+2));
	double dx = pMacGrid_->get_cell_sizex();
	double dy = pMacGrid_->get_cell_sizey();
	double dz = pMacGrid_->get_cell_sizez();
	double h = 2*dx;
	
	//~ for(int k = 0; k < L; ++k){
		//~ for(int j = 0; j < M; ++j){
			//~ for(int i = 0; i < N; ++i){
				//~ Eigen::Vector3d x_avrg_num = Eigen::Vector3d::Zero();
				//~ double x_avrg_den = 0;
				//~ double r_avrg_num = 0;
				//~ double r_avrg_den = 0;
				//~ int index = i + j*N + k*N*M;
				//~ Eigen::Vector3d cell_pos = Eigen::Vector3d(i*dx, j*dy, k*dz);
				//~ for(unsigned it_particle = 0; it_particle < num_particles_; ++it_particle){
					//~ Eigen::Vector3d particle_pos = (pparticles_+it_particle)->get_position();
					//~ double temp = (cell_pos - particle_pos).norm()/h;
					//~ std::cout << particle_pos << std::endl;
					//~ std::cout << cell_pos << std::endl;
					//~ std::cout << temp << std::endl << std::endl;
					//~ if (temp < 1){
						//~ //double W_surf = (1-temp*temp)*(1-temp*temp)*(1-temp*temp);
						//~ double W_surf = (1-temp)*(1-temp)*(1-temp);
						//~ //double W_surf = 1-temp*temp*temp;
						//~ x_avrg_num += W_surf*particle_pos;
						//~ x_avrg_den += W_surf;
						//~ r_avrg_num += W_surf*temp;
						//~ r_avrg_den += W_surf;
					//~ }	
				//~ }
				//~ Eigen::Vector3d x_avrg = x_avrg_num/x_avrg_den;
				//~ double r_avrg = r_avrg_num/r_avrg_den;
				//~ //*(plevel_set_ + index) = (cell_pos - x_avrg).norm() - r_avrg; 
				//~ if(x_avrg_den != 0)
					//~ plevel_set_[index] = 2*((cell_pos - x_avrg).norm() - r_avrg);
				//~ else{
					//~ if(pMacGrid_->is_fluid(i,j,k))
						//~ plevel_set_[index] = -100;
					//~ else
						//~ plevel_set_[index] = 100;
				//~ }
			//~ }
		//~ }
	//~ }
	for(unsigned it_particle = 0; it_particle < num_particles_; ++it_particle){
		Eigen::Vector3d particle_pos = (pparticles_+it_particle)->get_position();
		Eigen::Vector3d init_cell = pMacGrid_->index_from_coord(particle_pos[0], particle_pos[1], particle_pos[2]);
		for(int i = init_cell(0) - 2; i <= init_cell(0) + 2 && i >= 0 && i < N; ++i){
			for(int j = init_cell(1) - 2; j <= init_cell(1) + 2 && j >= 0 && j < M; ++j){
				for(int k = init_cell(2) - 2; k <= init_cell(2) + 2 && k >= 0 && k < L; ++k){
					Eigen::Vector3d cell = Eigen::Vector3d(i*dx,j*dy,k*dz);
					int index = i + j*N + k*N*M;
					double temp = (cell - particle_pos).norm()/h;
					if (temp < 1){
						double W_surf = (1-temp)*(1-temp)*(1-temp);
						//double W_surf = (1-temp*temp)*(1-temp*temp)*(1-temp*temp);
						//double W_surf = 1-temp*temp*temp;
						x_avrg_num[index] += W_surf*particle_pos;
						*(r_avrg_num + index) += W_surf*temp;
						*(den + index) += W_surf;
					}
				}
			}
		}
	}
	
	for(int k = 0; k < L; ++k){
		for(int j = 0; j < M; ++j){
			for(int i = 0; i < N; ++i){
				int index1 = (i) + (j)*N + (k)*N*M;
				double temp = *(den+index1);
				Eigen::Vector3d x_avrg = x_avrg_num[index1]/temp;
				double r_avrg = *(r_avrg_num+index1)/temp;
				Eigen::Vector3d cell_pos = Eigen::Vector3d((i)*dx, (j)*dy, (k)*dz);
				if(*(den+index1) != 0)
					plevel_set_[index1] = 2*((cell_pos - x_avrg).norm() - r_avrg);
				else{
					if(pMacGrid_->is_fluid(i,j,k))
						plevel_set_[index1] = -100;
					else
						plevel_set_[index1] = 100;
				}
			}
		}
	}
	//~ for(int k = 0; k < L+2; ++k){
		//~ for(int j = 0; j < M+2; ++j){
			//~ for(int i = 0; i < N+2; ++i){
				//~ if(i == 0 || i == N+1 || j == 0 || j == M+1 || k == 0 || k == L+1){
					//~ int index = i + j*(N+2) + k*(N+2)*(M+2);
					//~ plevel_set_[index] = 100;
				//~ }
				//~ else{
					//~ int index1 = (i-1) + (j-1)*N + (k-1)*N*M;
					//~ int index2 = i + j*(N+2) + k*(N+2)*(M+2);
					
					//~ double temp = *(den+index1);
					//~ Eigen::Vector3d x_avrg = x_avrg_num[index1]/temp;
					//~ double r_avrg = *(r_avrg_num+index1)/temp;
					//~ Eigen::Vector3d cell_pos = Eigen::Vector3d((i-1)*dx, (j-1)*dy, (k-1)*dz);
					//~ if(*(den+index1) != 0)
						//~ plevel_set_[index2] = 2*((cell_pos - x_avrg).norm() - r_avrg);
					//~ else{
						//~ if(pMacGrid_->is_fluid(i-1,j-1,k-1))
							//~ plevel_set_[index2] = -100;
						//~ else
							//~ plevel_set_[index2] = 100;
					//~ }
				//~ }
			//~ }
		//~ }
	//~ }
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
	MeshExporter::level_set_easy();
	igl::copyleft::marching_cubes(plevel_set_, points_, nx, ny, nz, vertices_, faces_);

	// Export mesh to OBJ file
	igl::writeOBJ(filename.str(), vertices_, faces_);
	std::cout << "Exported " << filename.str() << std::endl;

	num_exported_++;
}
