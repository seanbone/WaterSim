#include "MeshExporter.h"

MeshExporter::MeshExporter(Mac3d* Grid, Particle* particles, const int n)
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
