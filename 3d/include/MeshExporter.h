#ifndef MESHEXPORTER_H
#define MESHEXPORTER_H

#include "Mac3d.h"
#include "Eigen/Dense"
#include "Particles.h"
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
		Particles& particles_;
		
		//number of fluid particles in the list particles_
		unsigned num_particles_;
		
		//counter for the number of exported meshes
		unsigned num_exported_ = 0;

		// Name of folder to export file to
		std::string folder_ = "out_meshes/";
		std::string file_prefix_ = "mesh_";

		// grid spacing (initialized from pMacGrid_->get_cell_size{x,y,z}())
		const double dx, dy, dz;
		const double dx_r, dy_r, dz_r;

		// number of grid cells, per direction
		// (initialized from pMacGrid_->get_num_cells_{x,y,z}_())
		const int N, M, L;
		const double sizex_, sizey_, sizez_;
		const double r_avrg;


		// cell spacing (dx), reciprocal of square of the same number
		const double h, h_sq;
		const double sdf_trivial;
		double* points_d;
	public:
		//pointer to the values of the level set function
		double* plevel_set_array;

		using level_set_type = Eigen::Map<Eigen::VectorXd>;
		level_set_type plevel_set_map;
		using sdf_points_type = Eigen::Map<Eigen::MatrixXd>;
		sdf_points_type points_map;
		
		//vector which contains the numerator of the value x_average
		//used by the level set function
		std::vector<Eigen::Vector3d> x_avrg_num;
		double* x_avrg_num_array;
		
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
		MeshExporter(Mac3d* Grid, Particles& particles);
		
		/**Destructor
		 */
		~MeshExporter(){
			delete[] den;
			delete[] points_d;
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
