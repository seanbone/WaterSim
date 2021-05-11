#ifndef MAC3D_H
#define MAC3D_H

#include <Eigen/Dense>	//used for the Eigen dense matrix
#include <iostream> 	//used for input and output
#include <iomanip>		//used for input and output
#include <Eigen/Sparse>	//used for the matrix A
#include <vector>		//used for std::vector
#include <algorithm>	//std::fill
#include <cassert>		//assertions

class Mac3d{
	public:
		using Triplet_t = Eigen::Triplet<double>; 
		using Pair_t = std::pair<int, int>;

		using cellIdx_t = unsigned int;
		using globalCellIdx_t = unsigned int;

		//------------------- GRID Properties --------------------------
		//number of cells respectively in x-direction,
		//y-direction, z-direction
		const unsigned N_; 
		const unsigned M_;
		const unsigned L_;
		
		//size (in meter) of the grid respectively in x-direction,
		//y-direction, z-direction
		const double sizex_;
		const double sizey_;
		const double sizez_;
		
		//size (in meter) of one cell respectively in x-direction,
		//y-direction, z-direction
		const double cell_sizex_; 
		const double cell_sizey_;
		const double cell_sizez_;
		
	private:
		//---------------- PHYSICAL VALUES objects ---------------------
		//pointer to array for the pressure
		double* ppressure_;
		
		//pointer to array for the velocities respecitvely in x-direction
		//y-direction, z-direction
		double* pu_;
		double* pv_;
		double* pw_;
		
		//pointer to the temporary copy of velocity field u*, v*, w*
		//respectively in x-direction, y-direction, z-direction
		double* pu_star_;
		double* pv_star_;
		double* pw_star_;
		
		//pointer to array for specifing if a cell is solid (1) or not(0)
		bool* psolid_;
		
		//pointer to array for specifing if a cell contains fluid (1) or not(0)
		bool* pfluid_;
		
		//pointer to a std::vector which contians the triplets for
		//the diagonal of the matrix A, used to solve the pressures
		std::vector<Triplet_t> A_diag_;
		
		//pointer to the weights for particle-to-grid respectively for u, v, w
		double* pweights_u_;
		double* pweights_v_;
		double* pweights_w_;

	public:
		
		/** Default Constructor
		*/
		Mac3d()
			: N_(0), M_(0), L_(0), sizex_(0), sizey_(0), sizez_(0), cell_sizex_(0), cell_sizey_(0), cell_sizez_(0){}
		
		/**Constructor
		 * Params:
		 * - n, m and l represents the number of grid cells in each direction
		 * - dx, dy and dz are the real lengths of the grid in each direction (in meter)
		 */
		Mac3d(const unsigned n, const unsigned m, const unsigned l, const double dx, const double dy, const double dz);

		/** Destructor
		 */
		~Mac3d(){
			delete[] ppressure_;
			delete[] pu_;
			delete[] pv_;
			delete[] pw_;
			delete[] pu_star_;
			delete[] pv_star_;
			delete[] pw_star_;
			delete[] psolid_;
			delete[] pfluid_;
			delete[] pweights_u_;
			delete[] pweights_v_;
			delete[] pweights_w_;
		}

		/**Initialize the arrays of the class to zero, in particular:
		 * ppressure_, pu_, pv_, pw_, pu_star_, pv_star_, pw_star_,
		 * psolid_, pfluid_, pweights_u_,  pweights_v_, pweights_w_.
		 */
		void initArrays();

		/**Initialize the diagonal of the pressure matrix A 
		 */
		void initAdiag();
		
		//------------------------- GETTERS ----------------------------
		/**Return the size in meter of the grid in every direction 
		 */
		Eigen::Vector3d get_grid_size() const;
		
		/**Return the x-velocity in the mathematical point (i-1/2, j, k)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		 double get_u(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return the y-velocity in the mathematical point (i, j-1/2, k)
		 * * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_v(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/** Return the z-velocity in the mathematical point (i, j, k-1/2)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_w(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return the intermediate x-velocity u* in the mathematical point (i-1/2, j, k)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_u_star(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return the intermediate y-velocity v* in the mathematical point (i, j-1/2, k)
		 * * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_v_star(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/** Return the intermediate z-velocity w* in the mathematical point (i, j, k-1/2)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_w_star(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return the pressure in the mathematical point(i,j,k)
		 * Params:
		 * - i, j, k indicate in which position of the grid the pressure is asked.
		 */
		double get_pressure(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return if the cell with center (i,j,k) is a solid cell
		 * Params:
		 * - i, j, k indicate for which cell the physical property is asked.
		 */
		bool is_solid(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return if the cell with center (i,j,k) is a fluid cell
		 * Params:
		 * - i, j, k indicate for which cell the physical property is asked.
		 */
		bool is_fluid(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return if the cell with center (i,j,k) is empty (not solid && not fluid)
		 * Params:
		 * - i, j, k indicate for which cell the physical property is asked.
		 */
		bool is_empty(const unsigned i, const unsigned j, const unsigned k = 0);

		/**Return the number of cells of the grid in x-direction
		 */
		unsigned get_num_cells_x() const;
		
		/**Return the number of cells of the grid in y-direction
		 */
		unsigned get_num_cells_y() const;
		
		/**Return the number of cells of the grid in z-direction
		 */
		unsigned get_num_cells_z() const;
		
		/**Return the total number of cells of the grid
		 */
		unsigned get_num_cells() const;
		
		/**Return the dimension of one cell in x-direction in meter
		 */
		double get_cell_sizex() const;
		
		/**Return the dimension of one cell in y-direction in meter
		 */
		double get_cell_sizey() const;
		
		/**Return the dimension of one cell in z-direction in meter
		 */
		double get_cell_sizez() const;
		
		/** Return a const-reference to the diagonal of the pressure's matrix A
		 */
		const std::vector< Triplet_t >& get_a_diag();
		
		/**Return the weights for u in the mathematical point (i-1/2, j, k)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_weights_u(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return the weights for v in the mathematical point (i, j-1/2, k)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_weights_v(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Return the weights for v in the mathematical point (i, j, k-1/2)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is asked.
		 */
		double get_weights_w(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Assign the right x-indices for the interpolation of the velocities.
		 * Params:
		 * - x is the x coordinate of the point in which the interpolation is required
		 * - indices_x is the index in x-driection which denotes the cell
		 * 	 in which the point for the interpolation lies
		 * - ix0 is assigned by the function and refers to the indices in 
		 *   https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - ix1 is assigned by the function and refers to the indices in
		 *   https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - x0 is assigned by the function and is (ix0-0.5)*cell_sizex_
		 * - x1 is assigned by the function and is (ix1-0.5)*cell_sizex_
		 */
		void assign_x(double x, int indices_x, int& ix0, int& ix1, double& x0, double& x1);
		
		/**Assign the right y-indices for the interpolation of the velocities.
		 * Params:
		 * - y is the y coordinate of the point in which the interpolation is required
		 * - indices_y is the index in y-driection which denotes the cell
		 * 	 in which the point for the interpolation lies
		 * - iy0 is assigned by the function and refers to the indices in 
		 *   https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - iy1 is assigned by the function and refers to the indices in
		 *   https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - y0 is assigned by the function and is (ix0-0.5)*cell_sizex_
		 * - y1 is assigned by the function and is (ix1-0.5)*cell_sizex_
		 */
		void assign_y(double y, int indices_y, int& iy0, int& iy1, double& y0, double& y1);
		
		/**Assign the right z-indices for the interpolation of the velocities.
		 * Params:
		 * - z is the z coordinate of the point in which the interpolation is required
		 * - indices_z is the index in z-driection which denotes the cell
		 * 	 in which the point for the interpolation lies
		 * - iz0 is assigned by the function and refers to the indices in 
		 *   https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - iz1 is assigned by the function and refers to the indices in
		 *   https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - z0 is assigned by the function and is (ix0-0.5)*cell_sizex_
		 * - z1 is assigned by the function and is (ix1-0.5)*cell_sizex_
		 */
		void assign_z(double z, int indices_z, int& iz0, int& iz1, double& z0, double& z1);
		
		/**Assign the values xd, yd, zd explained in https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * Params:
		 * - xd, yd, zd are assigned by the function
		 * - x0, x1, y0, y1, z0, z1 describe the points in which the values 
		 *   of the velocities are taken in order to interpolate the velocity.
		 *   For better description see https://en.wikipedia.org/wiki/Trilinear_interpolation
		 * - x, y, z are the coordinate of the point in which the interpolation has to be done
		 */
		void assign_d(double& xd, double& yd, double& zd, double x0, double x1, double y0, double y1, double z0, double z1, double x, double y, double z);
		
		/**Return the interpolation of u or u* in the point (x,y,z)
		 * Params:
		 * - x, y, z are the coordinate in which the interpolation is made
		 * - use_u_star indicate if the velocity u (use_u_star == 0) or the 
		 * 	 velocity u* (use_u_star == 1) is used for the interpolation
		 */
		double get_interp_u(double x, double y, double z, const bool use_u_star = false);
		
		/**Return the interpolation of v or v* in the point (x,y,z)
		 * Params:
		 * - x, y, z are the coordinate in which the interpolation is made
		 * - use_v_star indicate if the velocity v (use_v_star == 0) or the 
		 * 	 velocity v* (use_v_star == 1) is used for the interpolation
		 */
		double get_interp_v(double x, double y, double z, const bool use_v_star = false);
		
		/**Return the interpolation of w or w* in the point (x,y,z)
		 * Params:
		 * - x, y, z are the coordinate in which the interpolation is made
		 * - use_w_star indicate if the velocity w (use_w_star == 0) or the 
		 * 	 velocity w* (use_w_star == 1) is used for the interpolation
		 */
		double get_interp_w(double x, double y, double z, const bool use_w_star = false);
		
		//------------------------ SETTERS -----------------------------
		/**Set the x-velocity in the mathematical point (i-1/2, j, k) at value
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is to be set
		 * - value is the value to which the velocity is set
		 */
		void set_u(const unsigned i, const unsigned j, const unsigned k, double value); 
		
		/**Set the y-velocity in the mathematical point (i,j-1/2,k) at value
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is to be set
		 * - value is the value to which the velocity is set
		 */
		void set_v(const unsigned i, const unsigned j, const unsigned k, double value);
		
		/**Set the z-velocity in the mathematical point (i,j,k-1/2) at value
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is to be set
		 * - value is the value to which the velocity is set
		 */
		void set_w(const unsigned i, const unsigned j, const unsigned k, double value);
		
		/**Set pu_star_ (used in validation)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is to be set
		 * - value is the value to which the velocity is set
		 */
		void set_u_star(const unsigned i, const unsigned j, const unsigned k, double value); 
		
		/**Set pv_star_ (used in validation)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is to be set
		 * - value is the value to which the velocity is set
		 */
		void set_v_star(const unsigned i, const unsigned j, const unsigned k, double value);
		
		/**Set pw_star_ (used in validation)
		 * Params:
		 * - i, j, k indicate in which position of the grid the velocity 
		 * 	 is to be set
		 * - value is the value to which the velocity is set
		 */
		void set_w_star(const unsigned i, const unsigned j, const unsigned k, double value);

		/**Copy velocity field to temporary copy arrays
		 */		
		void set_uvw_star();

		/**Set the pressure in the mathematical point(i,j,k)
		 * Params:
		 * - i, j, k indicate in which position of the grid the pressure 
		 * 	 is to be set
		 * - value is the value to which the pressure is set
		 */
		void set_pressure(const unsigned i, const unsigned j, const unsigned k, double value);
		
		/**Reset all the pressure with the values contained in p
		 * Params:
		 * - p is the vector which cointains the new value of the pressures
		 */
		void set_pressure(const Eigen::VectorXd& p);
		
		/**Set the weight for u in the mathematical point (i-1/2, j, k) at value
		 * Params:
		 * - i, j, k indicate in which position of the grid the weight_u 
		 * 	 is set
		 * - value is the value at which the weight is set
		 */
		void set_weights_u(const unsigned i, const unsigned j, const unsigned k, double value);
		
		/**Set the weight for v in the mathematical point (i,j-1/2,k) at value
		 * Params:
		 * - i, j, k indicate in which position of the grid the weight_v
		 * 	 is set
		 * - value is the value at which the weight is set
		 */
		void set_weights_v(const unsigned i, const unsigned j, const unsigned k, double value);
		
		/**Set the weight for w in the mathematical point (i,j,k-1/2) at value
		 * Params:
		 * - i, j, k indicate in which position of the grid the weight_w
		 * 	 is set
		 * - value is the value at which the weight is set
		 */
		void set_weights_w(const unsigned i, const unsigned j, const unsigned k, double value);
				
		/**Set all grid velocities to zero for particle-to-grid transfer
		 */
		void set_velocities_to_zero();
		
		/**Set all weights to zero for particle-to-grid transfer
		 */
		void set_weights_to_zero();
		
		/**Set the cell with center (i,j,k) as a solid cell
		 * Prams:
		 * - i, j, k are the indices of the cell to be set as solid
		 */
		void set_solid(const unsigned i, const unsigned j, const unsigned k = 0);

		/**Set the cell with center (i,j,k) as a solid cell
		 * Prams:
		 * - i, j, k are the indices of the cell to be set as solid
		 */
		void set_fluid(const unsigned i, const unsigned j, const unsigned k = 0);
		
		/**Reset all cells to not contain any fluid
		 */
		void reset_fluid();

		/** Return the indices in which the point with coordinate (x,y,z)
		 * lies as a Eigen::Vector3d
		 * Params:
		 * - x,y,z are the coordinate of the point 
		 */
		Eigen::Vector3d index_from_coord(const double x, const double y, const double z);
};
#endif
