#include "MeshExporter.h"

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
	plevel_set_.resize(N*M*L);
	double dx = pMacGrid_->get_cell_sizex();
	double dy = pMacGrid_->get_cell_sizey();
	double dz = pMacGrid_->get_cell_sizez();
	double h = 2*dx;
	double r = h;
	//double r = 0.25*dx;
	//~ Eigen::Vector3d x_avrg_num = Eigen::Vector3d::Zero();
	//~ double x_avrg_den = 0;
	//~ double r_avrg_num = 0;
	//~ double r_avrg_den = 0;
	
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < M; ++j){
			for(int k = 0; k < L; ++k){
				Eigen::Vector3d x_avrg_num = Eigen::Vector3d::Zero();
				double x_avrg_den = 0;
				double r_avrg_num = 0;
				double r_avrg_den = 0;
				int index = i + j*N + k*N*M;
				Eigen::Vector3d cell_pos = Eigen::Vector3d(i*dx, j*dy, k*dz);
				for(unsigned it_particle = 0; it_particle < num_particles_; ++it_particle){
					Eigen::Vector3d particle_pos = (pparticles_+it_particle)->get_position();
					double temp = (cell_pos - particle_pos).norm()/h;
					if (temp < 1){
						//double W_surf = (1-temp*temp)*(1-temp*temp)*(1-temp*temp);
						double W_surf = (1-temp)*(1-temp)*(1-temp);
						//double W_surf = 1-temp*temp*temp;
						x_avrg_num += W_surf*particle_pos;
						x_avrg_den += W_surf;
						r_avrg_num += W_surf*temp;
						r_avrg_den += W_surf;
					}	
				}
				Eigen::Vector3d x_avrg = x_avrg_num/x_avrg_den;
				double r_avrg = r_avrg_num/r_avrg_den;
				//*(plevel_set_ + index) = (cell_pos - x_avrg).norm() - r_avrg; 
				if(x_avrg_den != 0)
					plevel_set_[index] = 2*((cell_pos - x_avrg).norm() - r_avrg);
				else{
					if(pMacGrid_->is_fluid(i,j,k))
						plevel_set_[index] = -100;
					else
						plevel_set_[index] = 100;
				}
			}
		}
	}
	
	//~ for(int i = 0; i < N; ++i){
		//~ for(int j = 0; j < M; ++j){
			//~ for(int k = 0; k < L; ++k){
				//~ int index = i + j*N+ k*N*M;
				//~ if(k == 0)
					//~ std::cout << i << " " << j << " " << k << "   " << plevel_set_[index] << std::endl;
			//~ }
		//~ }
	//~ }
}
