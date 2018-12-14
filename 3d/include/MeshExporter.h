#ifndef MESHEXPORTER_H
#define MESHEXPORTER_H

#include "Mac3d.h"
#include "Particle.h"
#include "Eigen/Dense"
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>

class MeshExporter{
	private:
		Mac3d* pMacGrid_;
		Particle* pparticles_;
		unsigned num_particles_;
	
	public:
		Eigen::VectorXd plevel_set_;
		Eigen::MatrixXd points_;
		Eigen::MatrixXd vertices_;
		Eigen::MatrixXi faces_;
		
		MeshExporter(){}
		
		MeshExporter(Mac3d* Grid, Particle* particles, const int n)
			:pMacGrid_(Grid), pparticles_(particles), num_particles_(n){
			
			int N = pMacGrid_->get_num_cells_x();
			int M = pMacGrid_->get_num_cells_y();
			int L = pMacGrid_->get_num_cells_z();
			points_.resize(N*M*L, 3);
			double dx = pMacGrid_->get_cell_sizex();
			double dy = pMacGrid_->get_cell_sizey();
			double dz = pMacGrid_->get_cell_sizez();
			
			for(int k = 0; k < L; ++k){
				for(int j = 0; j < M; ++j){
					for(int i = 0; i < N; ++i){
						int index = i + j*N + k*N*M;
						Eigen::RowVector3d temp = Eigen::RowVector3d(i*dx, j*dy, k*dz);
						points_.row(index) = temp;
					}
				}
			}	
		}
		
		~MeshExporter(){
			num_particles_ = 0;
		}
		
		void level_set();
		void level_set_easy();
};

#endif //MESHEXPORTER_H
