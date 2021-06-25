#ifndef MAC3D_H
#define MAC3D_H

#include <Eigen/Dense>	//used for the Eigen dense matrix
#include <iostream> 	//used for input and output
#include <iomanip>		//used for input and output
#include <Eigen/Sparse>	//used for the matrix A
#include <vector>		//used for std::vector
#include <algorithm>	//std::fill
#include <cassert>		//assertions
#include <cstdlib>
#include <stdlib.h>

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
		enum GRID { GRID_U, GRID_V, GRID_W };

		/**
		 * Used to tell grid_interpolate whether to interpolate just on e.g. U or also U_star
		 */
		enum INTERPOLATION_MODE { INTERPOLATE_ONE, INTERPOLATE_BOTH };

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
		// the same information stored in three vectors
		std::vector<double> A_diag_val;
		
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
		bool is_fluid(const unsigned i, const unsigned j, const unsigned k = 0) const;
		
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
		const std::vector< Triplet_t >& get_a_diag() const;
		

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



		/********** INTERPOLATION METHODS **********/


		/**
		 * Interpolate a velocity component from the MAC grid.
		 * Alias for grid_interpolate<grid_name, INTERPOLATE_BOTH>(...)
		 * @tparam grid_name	The grid to interpolate from.
		 * @param pos_x 		Position x to interpolate to.
		 * @param pos_y			Position y to interpolate to
		 * @param pos_z			Position z to interpolate to
		 * @return 				Interpolated value.
		 */
		template<GRID grid_name>
		inline double grid_interpolate(double pos_x, double pos_y, double pos_z) {
			double ret;
			std::tie(ret, std::ignore) = grid_interpolate<grid_name, INTERPOLATE_ONE>(pos_x, pos_y, pos_z);
			return ret;
		}

		/**
		 * Interpolate a velocity component from the MAC grid.
		 * If interpolation_mode == INTERPOLATE_BOTH, the second component of the return
		 * value will be the interpolated value on the _star field corresponding to the velocity field
		 * defined in grid_name.
		 * @tparam grid_name			The grid to interpolate from.
		 * @tparam interpolation_mode
		 * @param pos_x 				Position x to interpolate to.
		 * @param pos_y					Position y to interpolate to
		 * @param pos_z					Position z to interpolate to
		 * @return 						Interpolated values.
		 */
		template<GRID grid_name, INTERPOLATION_MODE interpolation_mode>
		std::pair<double, double> grid_interpolate(double pos_x, double pos_y, double pos_z) {
			static_assert(grid_name == GRID_U || grid_name == GRID_V || grid_name == GRID_W,
			              "[Mac3d::grid_interpolate] Invalid template arguments!");

			// We are working on the staggered grid - that is, the grid where the
			// velocities are not at the center of cell faces, but at the
			// intersection points of the grid.
			// This means the grid will have different dimensions and offsets depending
			// on which one we are working on, given by template parameter grid_name.
			unsigned nx = N_;
			unsigned ny = M_;
			unsigned nz = L_;
			double* g;
			double* g_star;
			double offset_x = 0;
			double offset_y = 0;
			double offset_z = 0;

			// `if constexpr` is evaluated at compile time
			if constexpr(grid_name == GRID_U) {
				nx += 1;
				g = pu_;
				offset_x = -0.5 * cell_sizex_;
			} else if constexpr(grid_name == GRID_V) {
				ny += 1;
				g = pv_;
				offset_y = -0.5 * cell_sizey_;
			} else if constexpr(grid_name == GRID_W) {
				nz += 1;
				g = pw_;
				offset_z = -0.5 * cell_sizez_;
			}

			if constexpr(interpolation_mode == INTERPOLATE_BOTH) {
				if constexpr(grid_name == GRID_U) {
					g_star = pu_star_;
				} else if constexpr(grid_name == GRID_V) {
					g_star = pv_star_;
				} else if constexpr(grid_name == GRID_W) {
					g_star = pw_star_;
				}
			}

			double pos_x_offset = pos_x - offset_x;
			double pos_y_offset = pos_y - offset_y;
			double pos_z_offset = pos_z - offset_z;

			auto cell_x = unsigned(pos_x_offset * rcell_sizex_);
			auto cell_y = unsigned(pos_y_offset * rcell_sizey_);
			auto cell_z = unsigned(pos_z_offset * rcell_sizez_);

			double min_x = offset_x + cell_sizex_*cell_x;
			double min_y = offset_y + cell_sizey_*cell_y;
			double min_z = offset_z + cell_sizez_*cell_z;

			double ret, ret_star = 0;

			bool inside_left = pos_x_offset > 0;
			bool inside_right = pos_x_offset < (nx-1) * cell_sizex_;
			bool inside_bottom = pos_y_offset > 0;
			bool inside_top = pos_y_offset < (ny-1) * cell_sizey_;
			bool inside_front = pos_z_offset > 0;
			bool inside_back = pos_z_offset < (nz-1) * cell_sizez_;

			bool inside = inside_left && inside_right && inside_bottom && inside_top && inside_front && inside_back;
			if (inside) {
				// If we're inside the domain, perform trilinear interpolation
				unsigned i000 = (cell_x) + nx * (cell_y) + ny * nx * (cell_z);
				unsigned i001 = (cell_x) + nx * (cell_y) + ny * nx * (cell_z + 1);
				unsigned i010 = (cell_x) + nx * (cell_y + 1) + ny * nx * (cell_z);
				unsigned i011 = (cell_x) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);
				unsigned i100 = (cell_x + 1) + nx * (cell_y) + ny * nx * (cell_z);
				unsigned i101 = (cell_x + 1) + nx * (cell_y) + ny * nx * (cell_z + 1);
				unsigned i110 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z);
				unsigned i111 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);

				double alpha = (pos_x - min_x) * rcell_sizex_;
				double beta  = (pos_y - min_y) * rcell_sizey_;
				double gamma = (pos_z - min_z) * rcell_sizez_;
				ret = trilinear_interpolation_normalized(g[i000], g[i001], g[i010], g[i011],
				                                         g[i100], g[i101], g[i110], g[i111],
				                                         alpha, beta, gamma);
				if constexpr(interpolation_mode == INTERPOLATE_BOTH) {
					ret_star = trilinear_interpolation_normalized(g_star[i000], g_star[i001], g_star[i010], g_star[i011],
												   g_star[i100], g_star[i101], g_star[i110], g_star[i111],
												   alpha, beta, gamma);
				}
				return std::make_pair(ret, ret_star);
			}

			// On the boundaries of the domain, handle special cases
			else {
				unsigned i0, i1;
				unsigned i00, i01, i10, i11;

				if (!inside_left) {
					if (!inside_top) {
						// linear interpolation on top-left edge
						i0 = nx*(ny-1) + nx*ny*cell_z;
						i1 = i0 + nx*ny;
						return do_lerp<interpolation_mode>(i0, i1, min_z, rcell_sizez_, pos_z, g, g_star);
					} else if (!inside_bottom) {
						// linear interpolation on bottom-left edge
						i0 = nx*ny*cell_z;
						i1 = i0 + nx*ny;
						return do_lerp<interpolation_mode>(i0, i1, min_z, rcell_sizez_, pos_z, g, g_star);
					} else if (!inside_front) {
						// linear interpolation on front-left edge
						i0 = nx*cell_y;
						i1 = i0 + nx;
						return do_lerp<interpolation_mode>(i0, i1, min_y, rcell_sizey_, pos_y, g, g_star);
					} else if (!inside_back) {
						// linear interpolation on back-left edge
						i0 = nx*cell_y + nx*ny*(nz-1);
						i1 = i0 + nx;
						return do_lerp<interpolation_mode>(i0, i1, min_y, rcell_sizey_, pos_y, g, g_star);
					} else {
						// Bilinear interpolation on left face
						i00 = nx*cell_y + nx*ny*cell_z;
						i01 = i00 + nx*ny;
						i10 = i00 + nx;
						i11 = i00 + nx + nx*ny;
						return do_blerp<interpolation_mode>(i00, i01, i10, i11, min_y, min_z,
										  rcell_sizey_, rcell_sizez_, pos_y, pos_z, g, g_star);
					}
				} else if (!inside_right) {
					if (!inside_top) {
						// Linear interpolation on right-top edge
						i0 = nx-1 + nx*(ny-1) + nx*ny*cell_z;
						i1 = i0 + nx*ny;
						return do_lerp<interpolation_mode>(i0, i1, min_z, rcell_sizez_, pos_z, g, g_star);
					} else if (!inside_bottom) {
						// Linear interpolation on right-bottom edge
						i0 = nx-1 + nx*ny*cell_z;
						i1 = i0 + nx*ny;
						return do_lerp<interpolation_mode>(i0, i1, min_z, rcell_sizez_, pos_z, g, g_star);
					} else if (!inside_front) {
						// Linear interpolation on right-front edge
						i0 = nx-1 + nx*cell_y;
						i1 = i0 + nx;
						return do_lerp<interpolation_mode>(i0, i1, min_y, rcell_sizey_, pos_y, g, g_star);
					} else if (!inside_back) {
						// Linear interpolation on right-back edge
						i0 = nx-1 + nx*cell_y + nx*ny*(nz-1);
						i1 = i0 + nx;
						return do_lerp<interpolation_mode>(i0, i1, min_y, rcell_sizey_, pos_y, g, g_star);
					} else {
						// Bilinear on right face
						i00 = nx-1 + nx*cell_y + nx*ny*cell_z;
						i01 = i00 + nx*ny;
						i10 = i00 + nx;
						i11 = i00 + nx + nx*ny;
						return do_blerp<interpolation_mode>(i00, i01, i10, i11, min_y, min_z,
						                                    rcell_sizey_, rcell_sizez_, pos_y, pos_z, g, g_star);
					}
				} // we now know that inside_left && inside_right
				else if (!inside_top) {
					if (!inside_front) {
						// Linear interpolation on top-front edge
						i0 = cell_x + nx*(ny-1);
						i1 = i0+1;
						return do_lerp<interpolation_mode>(i0, i1, min_x, rcell_sizex_, pos_x, g, g_star);
					} else if (!inside_back) {
						// Linear interpolation on top-back edge
						i0 = cell_x + nx*(ny-1) + nx*ny*(nz-1);
						i1 = i0+1;
						return do_lerp<interpolation_mode>(i0, i1, min_x, rcell_sizex_, pos_x, g, g_star);
					} else {
						// bilinear on top face
						i00 = cell_x + nx * (ny-1) + nx*ny*cell_z;
						i01 = i00 + 1;
						i10 = i00 + nx*ny;
						i11 = i10 + 1;
						return do_blerp<interpolation_mode>(i00, i01, i10, i11, min_z, min_x,
						                                    rcell_sizez_, rcell_sizex_, pos_z, pos_x, g, g_star);
					}
				} else if (!inside_bottom) {
					if (!inside_front) {
						// Linear interpolation on bottom-front edge
						i0 = cell_x;
						i1 = i0+1;
						return do_lerp<interpolation_mode>(i0, i1, min_x, rcell_sizex_, pos_x, g, g_star);
					} else if (!inside_back) {
						// Linear interpolation on bottom-back
						i0 = cell_x + nx*ny*(nz-1);
						i1 = i0+1;
						return do_lerp<interpolation_mode>(i0, i1, min_x, rcell_sizex_, pos_x, g, g_star);
					} else {
						// Bilinear interpolation on bottom face
						i00 = cell_x + nx*ny*cell_z;
						i01 = i00 + 1;
						i10 = i00 + nx*ny;
						i11 = i10 + 1;
						return do_blerp<interpolation_mode>(i00, i01, i10, i11, min_z, min_x,
						                                    rcell_sizez_, rcell_sizex_, pos_z, pos_x, g, g_star);
					}
				} else if (!inside_front) {
					// Bilinear interpolation on front face
					i00 = cell_x + nx*cell_y;
					i01 = i00 + nx;
					i10 = i00 + 1;
					i11 = i10 + nx;
					return do_blerp<interpolation_mode>(i00, i01, i10, i11, min_x, min_y,
					                                    rcell_sizex_, rcell_sizey_, pos_x, pos_y, g, g_star);
				} else {
					// Bilinear interpolation on back face
					i00 = cell_x + nx*cell_y + nx*ny*(nz-1);
					i01 = i00 + nx;
					i10 = i00 + 1;
					i11 = i10 + nx;
					return do_blerp<interpolation_mode>(i00, i01, i10, i11, min_x, min_y,
					                                    rcell_sizex_, rcell_sizey_, pos_x, pos_y, g, g_star);
				}
			}

		}

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


private:

	template<INTERPOLATION_MODE interpolation_mode>
	inline std::pair<double, double> do_lerp(const int i0, const int i1, const double min_pos,
										  const double r_size, const double pos,
										  const double* g, const double* g_star) const {
		double alpha = (pos - min_pos) * r_size;
		double ret = linear_interpolation_normalized(g[i0], g[i1], alpha);
		double ret_star = 0;
		if constexpr(interpolation_mode == INTERPOLATE_BOTH) {
			ret_star = linear_interpolation_normalized(g_star[i0], g_star[i1], alpha);
		}
		return std::make_pair(ret, ret_star);
	}

	template<INTERPOLATION_MODE interpolation_mode>
	inline std::pair<double, double> do_blerp(const int i00, const int i01, int i10, int i11,
										   const double min_x, const double min_y,
										   const double r_size_x, const double r_size_y,
										   const double pos_x, const double pos_y,
										   const double* g, const double* g_star) const {
		double alpha = (pos_x - min_x) * r_size_x;
		double beta = (pos_y - min_y) * r_size_y;
		double ret = bilinear_interpolation_normalized(g[i00], g[i01], g[i10], g[i11], alpha, beta);
		double ret_star = 0;
		if constexpr(interpolation_mode == INTERPOLATE_BOTH) {
			ret_star = bilinear_interpolation_normalized(g_star[i00], g_star[i01], g_star[i10], g_star[i11], alpha, beta);
		}
		return std::make_pair(ret, ret_star);
	}
};
#endif
