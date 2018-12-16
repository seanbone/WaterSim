#ifndef MESHEXPORTER_H
#define MESHEXPORTER_H

#include "Mac3d.h"
#include "Particle.h"
#include "Eigen/Dense"
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <sstream>
#include <iomanip>  
#include <sys/stat.h>

#include <igl/copyleft/marching_cubes.h>
#include <igl/writeOBJ.h>

class MeshExporter{
	private:
		Mac3d* pMacGrid_;
		Particle* pparticles_;
		unsigned num_particles_;
		unsigned num_exported_ = 0;

		// Name of folder to export file to
		std::string folder_ = "../out_meshes/";
		std::string file_prefix_ = "mesh_";
	
	public:
		Eigen::VectorXd plevel_set_;
		Eigen::MatrixXd points_;
		Eigen::MatrixXd vertices_;
		Eigen::MatrixXi faces_;
		std::vector<Eigen::Vector3d> x_avrg_num;
		double* r_avrg_num;
		double* den;
		
		MeshExporter(Mac3d* Grid, Particle* particles, const int n);
		
		~MeshExporter(){
			num_particles_ = 0;
			delete[] r_avrg_num;
			delete[] den;
		}
		
		void level_set_open();
		void level_set_easy();

		void export_mesh();
};

#endif //MESHEXPORTER_H
