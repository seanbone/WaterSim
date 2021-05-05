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
		//pointer to the MacGrid
		Mac3d* pMacGrid_;
		
		//pointer to the fluid particle
		Particle* pparticles_;
		
		//number of fluid particles in the list pparticles_
		unsigned num_particles_;
		
		//counter for the number of exported meshes
		unsigned num_exported_ = 0;

		// Name of folder to export file to
		std::string folder_ = "out_meshes/";
		std::string file_prefix_ = "mesh_";
		const double dx, dy, dz;
		const unsigned int N, M, L;
		const double sizex_, sizey_, sizez_;
		const double h, h_sq_r;
		const double weight_factor;
	public:
		//pointer to the values of the level set function
		Eigen::VectorXd plevel_set_;
		
		//vector which contains the numerator of the value x_average
		//used by the level set function
		std::vector<Eigen::Vector3d> x_avrg_num;
		
		//pointer to the list which contains the numerator of the value
		//r_average used by the level set function
		double* r_avrg_num;
		
		//pointer to the list of the denominator used by both x_average and 
		//r_average (used by the level set function)
		double* den;
		
		//matrix with the coordinate of the points in which the level set 
		//function is computed. Every row represents the coordinate of a point
		Eigen::MatrixXd points_;
		
		//matrix to store the vertices of the mesh (used by marching_cubes)
		Eigen::MatrixXd vertices_;
		
		//matrix to store the faces of the mesh (used by marching_cubes)
		Eigen::MatrixXi faces_;
		
		/**Constructor
		 * Params:
		 * - Grid is a pointer to a Mac3d grid
		 * - particles is a pointer to a list of Particle
		 * - n is the number of particles contained in particles
		 */
		MeshExporter(Mac3d* Grid, Particle* particles, const int n);
		
		/**Destructor
		 */
		~MeshExporter(){
			num_particles_ = 0;
			delete[] r_avrg_num;
			delete[] den;
		}
		
		/**Level set function:
		 * computes the value of the level set function at every point in 
		 * the points_ list. The value of the level set function is bigger 
		 * than 0 if the point in which is evaluated lies out of the fluid 
		 * and smaller than 0 else. At the point where the level set function 
		 * is exactly 0 lies the surface of the fluid.
		 */
		void level_set();
		
		/**Level set function "easy version":
		 * computes the value of the level set function at every point in 
		 * the points_ list. The value of the level set function is 1+
		 * if the point in which is evaluated lies out of the fluid 
		 * and -1 else. At the point where the level set function 
		 * is exactly 0 lies the surface of the fluid.
		 */
		void level_set_easy();

		/**
		 * Compute the mesh without exporting it to a file.
		 */
		void compute_mesh();
		
		/**
		 * Export the mesh in an OBJ file. Note that in order to update the mesh,
		 * compute_mesh must be called separately.
		 */
		void export_mesh();

		/**
		 * Get the computed mesh.
		 */
		 void get_mesh(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces) const;
};

#endif //MESHEXPORTER_H
