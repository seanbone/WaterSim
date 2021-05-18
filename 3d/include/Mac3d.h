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

		/**
		 * Index type for per-axis cell indices.
		 * Signed type to allow for negative offsets.
		 */
		using cellIdx_t = int;

		/**
		 * Index type for global cell indices.
		 */
		using globalCellIdx_t = int;

		/**
		 * Used to tell methods which grid to operate on.
		 */
		enum GRID { GRID_P, GRID_U, GRID_V, GRID_W, GRID_U_STAR, GRID_V_STAR, GRID_W_STAR };

		/**
		 * Offsets of the grids from the origin.
		 * The first entry of grid GRID will be at spatial coordinates GRID_OFFSETS[GRID].
		 */
		double GRID_OFFSETS[7][3];

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
		const double rcell_sizex_;
		const double rcell_sizey_;
		const double rcell_sizez_;

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
		
		/** Default Constructor
		*/
		Mac3d()
			: N_(0), M_(0), L_(0), 
			  sizex_(0), sizey_(0), sizez_(0), 
			  cell_sizex_(0), cell_sizey_(0), cell_sizez_(0), 
			  rcell_sizex_(0), rcell_sizey_(0), rcell_sizez_(0){}
		
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

		void index_from_coord( const double x, 
							   const double y, 
							   const double z, 
							   Mac3d::cellIdx_t &cell_idx_x, 
							   Mac3d::cellIdx_t &cell_idx_y, 
							   Mac3d::cellIdx_t &cell_idx_z );



		/********** WIP INTERPOLATION **********/


		/**
		 * Perform linear interpolation on the normalized [0, 1] domain.
		 * @param value0 	Value at position x=0
		 * @param value1 	Value at position x=1
		 * @param alpha 	Position in the [0, 1] space to interpolate to
		 */
		static inline double linear_interpolation_normalized(double value0, double value1, double alpha) {
			return value0 + alpha * (value1 - value0);
		}

		/**
		 * Perform bilinear interpolation on the normalized [0,1]x[0,1] domain.
		 * @param v00 	Value at position (0, 0)
		 * @param v01 	Value at position (0, 1)
		 * @param v10 	Value at position (1, 0)
		 * @param v11 	Value at position (1, 1)
		 * @param alpha	Position on the x-axis in [0, 1] space to interpolate to
		 * @param beta 	Position on the y-axis in [0, 1] space to interpolate to
		 */
		static inline double bilinear_interpolation_normalized(double v00, double v01,
		                                                       double v10, double v11,
		                                                       double alpha, double beta) {
			double v0 = (v00 + alpha*(v10 - v00));
			double v1 = (v01 + alpha*(v11 - v01));
			return v0 + beta*(v1 - v0);
		}


		/**
		 * Perform trilinear interpolation on the normalized [0,1]x[0,1]x[0,1] domain.
		 * @param v000	Value at position (0, 0, 0)
		 * @param v001	Value at position (0, 0, 1)
		 * @param v010	Value at position (0, 1, 0)
		 * @param v011	Value at position (0, 1, 1)
		 * @param v100	Value at position (1, 0, 0)
		 * @param v101	Value at position (1, 0, 1)
		 * @param v110	Value at position (1, 1, 0)
		 * @param v111	Value at position (1, 1, 1)
		 * @param alpha Position on the x-axis in [0, 1] space to interpolate to
		 * @param beta 	Position on the y-axis in [0, 1] space to interpolate to
		 * @param gamma Position on the z-axis in [0, 1] space to interpolate to
		 */
		static inline double trilinear_interpolation_normalized(double v000, double v001, double v010, double v011,
		                                                        double v100, double v101, double v110, double v111,
		                                                        double alpha, double beta, double gamma) {
			double v00 = v000 + alpha*(v100 - v000);
			double v01 = v001 + alpha*(v101 - v001);
			double v10 = v010 + alpha*(v110 - v010);
			double v11 = v011 + alpha*(v111 - v011);

			double v0 = v00 + beta*(v10 - v00);
			double v1 = v01 + beta*(v11 - v01);

			return v0 + gamma*(v1 - v0);
		}


		/**
		 * Perform linear interpolation on the general [min_pos, min_pos + 1/r_size] domain.
		 * @param v0 		Value at start of the domain, x = min_pos.
		 * @param v1 		Value at end of the domain, x = min_pos + 1/r_size.
		 * @param min_pos 	Start of domain.
		 * @param r_size	Reciprocal of the size of the domain.
		 * @param pos 		Position to interpolate to.
		 */
		static inline double linear_interpolation(double v0, double v1, double min_pos, double r_size, double pos) {
			double alpha = (pos - min_pos) * r_size;
			return linear_interpolation_normalized(v0, v1, alpha);
		}

		/**
		 * Perform bilinear interpolation on the general
		 * [x_min, x_min + 1/r_size_x] x [y_min, y_min + 1/r_size_y] domain.
		 * @param v00 		Value at position (x_min, y_min)
		 * @param v01  		Value at position (x_min, y_min + 1/r_size_y)
		 * @param v10  		Value at position (x_min + 1/r_size_x, y_min)
		 * @param v11 		Value at position (x_min + 1/r_size_x, y_min + 1/r_size_y)
		 * @param min_x 	Start of domain in X axis.
		 * @param min_y 	Start of domain in Y axis.
		 * @param r_size_x 	Reciprocal of the domain size on the X axis.
		 * @param r_size_y 	Reciprocal of the domain size on the Y axis.
		 * @param pos_x 	Position to interpolate to on X axis.
		 * @param pos_y 	Position to interpolate to on Y axis.
		 */
		static inline double bilinear_interpolation(double v00, double v01,
		                                            double v10, double v11,
		                                            double min_x, double min_y,
		                                            double r_size_x, double r_size_y,
		                                            double pos_x, double pos_y) {
			double alpha = (pos_x - min_x) * r_size_x;
			double beta  = (pos_y - min_y) * r_size_y;
			return bilinear_interpolation_normalized(v00, v01, v10, v11, alpha, beta);
		}

		/**
		 * Perform trilinear interpolation on the general domain
		 * [x_min, x_min + 1/r_size_x] x [y_min, y_min + 1/r_size_y] x [z_min, z_min + 1/r_size_z]
		 * @param min_x  	Start of domain in X axis.
		 * @param min_y  	Start of domain in Y axis.
		 * @param min_z   	Start of domain in Z axis.
		 * @param r_size_x 	Reciprocal of the domain size on the X axis.
		 * @param r_size_y 	Reciprocal of the domain size on the Y axis.
		 * @param r_size_z 	Reciprocal of the domain size on the Z axis.
		 * @param pos_x 	Position to interpolate to on X axis.
		 * @param pos_y 	Position to interpolate to on Y axis.
		 * @param pos_z 	Position to interpolate to on Z axis.
		 */
		static inline double trilinear_interpolation(double v000, double v001, double v010, double v011,
		                                             double v100, double v101, double v110, double v111,
		                                             double min_x, double min_y, double min_z,
		                                             double r_size_x, double r_size_y, double r_size_z,
		                                             double pos_x, double pos_y, double pos_z) {
			double alpha = (pos_x - min_x) * r_size_x;
			double beta  = (pos_y - min_y) * r_size_y;
			double gamma = (pos_z - min_z) * r_size_z;
			return trilinear_interpolation_normalized(v000, v001, v010, v011,
											          v100, v101, v110, v111,
											          alpha, beta, gamma);
		}


		template<GRID grid_name>
		double grid_interpolate(double pos_x, double pos_y, double pos_z) {
			// We are working on the staggered grid - that is, the grid where the pressures
			// (or velocities) are not the the center of the cells (cell faces), but at the
			// intersection points of the grid. Therefore we have one fewer cell on each axis.
			unsigned nx;
			unsigned ny;
			unsigned nz;
			double* g;
			// The different arrays have different dimensions
			// get_grid_properties is a funky template hack to have this evaluation done at compile time
			get_grid_properties<grid_name>(nx, ny, nz, g);

			const double& offset_x = GRID_OFFSETS[grid_name][0];
			const double& offset_y = GRID_OFFSETS[grid_name][1];
			const double& offset_z = GRID_OFFSETS[grid_name][2];

			double pos_x_offset = pos_x - offset_x;
			double pos_y_offset = pos_y - offset_y;
			double pos_z_offset = pos_z - offset_z;

			int cell_x = pos_x_offset * rcell_sizex_;
			int cell_y = pos_y_offset * rcell_sizey_;
			int cell_z = pos_z_offset * rcell_sizez_;

			int i000 = (cell_x    ) + nx * (cell_y    ) + ny * nx * (cell_z    );
			double min_x = offset_x + cell_sizex_*cell_x;
			double min_y = offset_y + cell_sizey_*cell_y;
			double min_z = offset_z + cell_sizez_*cell_z;

			bool inside_left = pos_x_offset > 0;
			bool inside_right = pos_x_offset < (nx-1) * cell_sizex_;
			bool inside_bottom = pos_y_offset > 0;
			bool inside_top = pos_y_offset < (ny-1) * cell_sizey_;
			bool inside_front = pos_z_offset > 0;
			bool inside_back = pos_z_offset < (nz-1) * cell_sizez_;

			bool inside = inside_left && inside_right && inside_bottom && inside_top && inside_front && inside_back;
			if (!inside) {
				if (!inside_left) {
					if (!inside_top) {
						// linear interpolation on top-left edge
						int i0 = nx*(ny-1) + nx*ny*cell_z;
						int i1 = i0 + nx*ny;
						return linear_interpolation(g[i0], g[i1], min_z, rcell_sizez_, pos_z);
					} else if (!inside_bottom) {
						// linear interpolation on bottom-left edge
						int i0 = nx*ny*cell_z;
						int i1 = i0 + nx*ny;
						//std::cout << "<<<<\n";
						//std::cout << i0 << std::endl;
						//std::cout << i1 << std::endl;
						//std::cout << g[i0] << std::endl;
						//std::cout << g[i1] << std::endl;
						//std::cout << "min_z " << min_z << std::endl;
						//std::cout << "rcell_size_z " << rcell_sizez_ << std::endl;
						//std::cout << "pos_z " << pos_z << std::endl;
						//std::cout << "<<<<\n";
						return linear_interpolation(g[i0], g[i1], min_z, rcell_sizez_, pos_z);
					} else if (!inside_front) {
						// linear interpolation on front-left edge
						int i0 = nx*cell_y;
						int i1 = i0 + nx;
						return linear_interpolation(g[i0], g[i1], min_y, rcell_sizey_, pos_y);
					} else if (!inside_back) {
						// linear interpolation on back-left edge
						int i0 = nx*cell_y + nx*ny*(nz-1);
						int i1 = i0 + nx;
						return linear_interpolation(g[i0], g[i1], min_y, rcell_sizey_, pos_y);
					} else {
						// Bilinear interpolation on left face
						std::cout << " * Bilinear on left face" << std::endl;
						int i00 = nx*cell_y + nx*ny*cell_z;
						int i01 = i00 + nx*ny;
						int i10 = i00 + nx;
						int i11 = i00 + nx + nx*ny;
						return bilinear_interpolation(g[i00], g[i01], g[i10], g[i11],
									min_z, min_y, rcell_sizez_, rcell_sizey_, pos_z, pos_y);
					}
				} else if (!inside_right) {
					if (!inside_top) {
						// Linear interpolation on right-top edge
						int i0 = nx-1 + nx*(ny-1) + nx*ny*cell_z;
						int i1 = i0 + nx*ny;
						return linear_interpolation(g[i0], g[i1], min_z, rcell_sizez_, pos_z);
					} else if (!inside_bottom) {
						// Linear interpolation on right-bottom edge
						int i0 = nx-1 + nx*ny*cell_z;
						int i1 = i0 + nx*ny;
						return linear_interpolation(g[i0], g[i1], min_z, rcell_sizez_, pos_z);
					} else if (!inside_front) {
						// Linear interpolation on right-front edge
						int i0 = nx-1 + nx*cell_y;
						int i1 = i0 + nx;
						return linear_interpolation(g[i0], g[i1], min_y, rcell_sizey_, pos_y);
					} else if (!inside_back) {
						// Linear interpolation on right-back edge
						int i0 = nx-1 + nx*cell_y + nx*ny*(nz-1);
						int i1 = i0 + nx;
						return linear_interpolation(g[i0], g[i1], min_y, rcell_sizey_, pos_y);
					} else {
						// Bilinear on right face
						std::cout << " * Bilinear on right face" << std::endl;
						int i00 = nx-1 + nx*cell_y + nx*ny*cell_z;
						int i01 = i00 + nx*ny;
						int i10 = i00 + nx;
						int i11 = i00 + nx + nx*ny;
						return bilinear_interpolation(g[i00], g[i01], g[i10], g[i11],
						                              min_z, min_y, rcell_sizez_, rcell_sizey_, pos_z, pos_y);
					}
				} // we now know that inside_left && inside_right
				else if (!inside_top) {
					if (!inside_front) {
						// Linear interpolation on top-front edge
						int i0 = cell_x + nx*(ny-1);
						int i1 = i0+1;
						return linear_interpolation(g[i0], g[i1], min_x, rcell_sizex_, pos_x);
					} else if (!inside_back) {
						// Linear interpolation on top-back edge
						int i0 = cell_x + nx*(ny-1) + nx*ny*(nz-1);
						int i1 = i0+1;
						return linear_interpolation(g[i0], g[i1], min_x, rcell_sizex_, pos_x);
					} else {
						// bilinear on top face
						int i00 = cell_x + nx * (ny-1) + nx*ny*cell_z;
						int i01 = i00 + 1;
						int i10 = i00 + nx*ny;
						int i11 = i10 + 1;
						//std::cout << "bilinear interp" << std::endl;
						//std::cout << "cell_x " << cell_x << std::endl;
						//std::cout << "cell_y " << cell_y << std::endl;
						//std::cout << "cell_z " << cell_z << std::endl;
						//std::cout << "min_x " << min_x << std::endl;
						//std::cout << "min_y " << min_y << std::endl;
						//std::cout << "min_z " << min_z << std::endl;
						//std::cout << "cell_sizex_ " << cell_sizex_ << std::endl;
						//std::cout << "cell_sizey_ " << cell_sizey_ << std::endl;
						//std::cout << "cell_sizez_ " << cell_sizez_ << std::endl;
						//std::cout << "pos_x " << pos_x << std::endl;
						//std::cout << "pos_y " << pos_y << std::endl;
						//std::cout << "pos_z " << pos_z << std::endl;
						//std::cout << i00 << std::endl;
						//std::cout << i01 << std::endl;
						//std::cout << i10 << std::endl;
						//std::cout << i11 << std::endl;
						//std::cout << g[i00] << std::endl;
						//std::cout << g[i01] << std::endl;
						//std::cout << g[i10] << std::endl;
						//std::cout << g[i11] << std::endl;
						return bilinear_interpolation(g[i00], g[i01], g[i10], g[i11],
						                              min_x, min_z, rcell_sizex_, rcell_sizez_, pos_x, pos_z);
					}
				} else if (!inside_bottom) {
					if (!inside_front) {
						// Linear interpolation on bottom-front edge
						int i0 = cell_x;
						int i1 = i0+1;
						return linear_interpolation(g[i0], g[i1], min_x, rcell_sizex_, pos_x);
					} else if (!inside_back) {
						// Linear interpolation on bottom-back
						int i0 = cell_x + nx*ny*(nz-1);
						int i1 = i0+1;
						return linear_interpolation(g[i0], g[i1], min_x, rcell_sizex_, pos_x);
					} else {
						// Bilinear interpolation on bottom face
						int i00 = cell_x + nx*ny*cell_z;
						int i01 = i00 + 1;
						int i10 = i00 + nx*ny;
						int i11 = i10 + 1;
						return bilinear_interpolation(g[i00], g[i01], g[i10], g[i11],
						                              min_x, min_z, rcell_sizex_, rcell_sizez_, pos_x, pos_z);
					}
				} else if (!inside_front) {
					// Bilinear interpolation on front face
					int i00 = cell_x + nx*cell_y;
					int i01 = i00 + nx;
					int i10 = i00 + 1;
					int i11 = i10 + nx;
					return bilinear_interpolation(g[i00], g[i01], g[i10], g[i11],
					                              min_x, min_y, rcell_sizex_, rcell_sizey_, pos_x, pos_y);
				} else {
					// Bilinear interpolation on back face
					int i00 = cell_x + nx*cell_y + nx*ny*(nz-1);
					int i01 = i00 + nx;
					int i10 = i00 + 1;
					int i11 = i10 + nx;
					return bilinear_interpolation(g[i00], g[i01], g[i10], g[i11],
					                              min_x, min_y, rcell_sizex_, rcell_sizey_, pos_x, pos_y);
				}
			}

			int i001 = (cell_x    ) + nx * (cell_y    ) + ny * nx * (cell_z + 1);
			int i010 = (cell_x    ) + nx * (cell_y + 1) + ny * nx * (cell_z    );
			int i011 = (cell_x    ) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);
			int i100 = (cell_x + 1) + nx * (cell_y    ) + ny * nx * (cell_z    );
			int i101 = (cell_x + 1) + nx * (cell_y    ) + ny * nx * (cell_z + 1);
			int i110 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z    );
			int i111 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);

			/*
			std::cout << "cell_x " << cell_x << std::endl;
			std::cout << "cell_y " << cell_y << std::endl;
			std::cout << "cell_z " << cell_z << std::endl;
			std::cout << i000 << std::endl;
			std::cout << i100 << std::endl;
			std::cout << i010 << std::endl;
			std::cout << i110 << std::endl;
			std::cout << i001 << std::endl;
			std::cout << i101 << std::endl;
			std::cout << i011 << std::endl;
			std::cout << i111 << std::endl;
			 */

			return trilinear_interpolation(g[i000], g[i001], g[i010], g[i011],
								           g[i100], g[i101], g[i110], g[i111],
								           min_x, min_y, min_z,
								           rcell_sizex_, rcell_sizey_, rcell_sizez_,
								           pos_x, pos_y, pos_z);
		}

		template<GRID grid_type>
		double grid_interpolate_marginal(double pos_x, double pos_y, double pos_z) {

		}

		template<GRID grid_type>
		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid) {
			get_grid_properties(nx, ny, nz, grid, Identity<grid_type>());
		}

	private:
		template<GRID grid_type>
		struct Identity { };

		template<GRID grid_type>
		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<grid_type>) {
			throw std::runtime_error("[Mac3d::get_grid_properties] Invalid grid type!");
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_P>) {
			nx = M_;
			ny = N_;
			nz = L_;
			grid = ppressure_;
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_U>) {
			nx = M_+1;
			ny = N_;
			nz = L_;
			grid = pu_;
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_V>) {
			nx = M_;
			ny = N_+1;
			nz = L_;
			grid = pv_;
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_W>) {
			nx = M_;
			ny = N_;
			nz = L_+1;
			grid = pw_;
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_U_STAR>) {
			nx = M_+1;
			ny = N_;
			nz = L_;
			grid = pu_star_;
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_V_STAR>) {
			nx = M_;
			ny = N_+1;
			nz = L_;
			grid = pv_star_;
		}

		inline void get_grid_properties(unsigned& nx, unsigned& ny, unsigned& nz, double*& grid, Identity<GRID_W_STAR>) {
			nx = M_;
			ny = N_;
			nz = L_+1;
			grid = pw_star_;
		}

	/********** END WIP INTERPOLATION **********/
};
#endif
