#ifndef MAC3D_H
#define MAC3D_H

#include <Eigen/Dense>	//used for the Eigen dense matrix
#include <utility>		//used for std::pair
#include <iostream> 	//used for input and output
#include <iomanip>		//used for input and output
#include <Eigen/Sparse>	//used for the matrix A
#include <vector>		//used for std::vector
#include <algorithm>	// std::fill
#include <cassert>		//assertions

class Mac3d{
	public:
		using Triplet_t = Eigen::Triplet<double>; 
		using Pair_t = std::pair<int, int>; //To access the component of the Pair use the methods "first" and "second"

	private:
		//PARAMETER for the grid
		//number of cells in x
		const unsigned N_; 
		//number of cells in y
		const unsigned M_;
		//number of cells in z
		const unsigned L_;
		//size (in meter) of the grid in x-direction
		const double sizex_;
		//size (in meter) of the grid in y-direction
		const double sizey_;
		//size (in meter) of the grid in z-direction
		const double sizez_;
		//size (in meter) of one cell in x-direction
		const double cell_sizex_; 
		//size (in meter) of one cell in y-direction
		const double cell_sizey_;
		//size (in meter) of one cell in z-direction
		const double cell_sizez_;
		
		//PHYSICAL VALUES objects
		//pointer to array for the pressure
		double* ppressure_;
		//pointer to array for the velocities in x-direction
		double* pu_;
		//pointer to array for the velocities in y-direction
		double* pv_;
		//pointer to array for the velocities in z-direction
		double* pw_;
		// Temporary copy of velocity field -> u^*, v^*
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
		// Weights for particle-to-grid for u
		double* pweights_u_;
		// Weights for particle-to-grid for v
		double* pweights_v_;
		// Weights for particle-to-grid for w
		double* pweights_w_;

	public:
		//CONSTRUCTORS
		//Default constructor
		Mac3d()
			: N_(0), M_(0), L_(0), sizex_(0), sizey_(0), sizez_(0), cell_sizex_(0), cell_sizey_(0), cell_sizez_(0){}
		
		//Constructor with the number of cells (in all the three directions) and the length of the box (in metres)
		Mac3d(const unsigned n, const unsigned m, const unsigned l, const double dx, const double dy, const double dz);

		// Initialize arrays to zero (pressure, u, v, w, etc)
		void initArrays();

		// Initialize A diagonal
		void initAdiag();
		
		//Destructor
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
		}
		
		//GETS
		// Get size of the full grid
		Eigen::Vector3d get_grid_size() const;
		//Get the x-velocity in the mathematical point (i-1/2,j, k) 
		double get_u(const unsigned i, const unsigned j, const unsigned k);
		//Get the y-velocity in the mathematical point (i,j-1/2, k)
		double get_v(const unsigned i, const unsigned j, const unsigned k);
		//Get the z-velocity in the mathematical point (i,j-1/2, k-1/2)
		double get_w(const unsigned i, const unsigned j, const unsigned k);
		// Equivalent for intermediate field
		double get_u_star(const unsigned i, const unsigned j, const unsigned k);
		double get_v_star(const unsigned i, const unsigned j, const unsigned k);
		double get_w_star(const unsigned i, const unsigned j, const unsigned k);
		//Get the pressure in the mathematical point(i,j,k)
		double get_pressure(const unsigned i, const unsigned j, const unsigned k = 0);
		//Return if the cell with center (i,j,k) is a solid cell
		bool is_solid(const unsigned i, const unsigned j, const unsigned k);
		//Return if the cell with center (i,j,k) is a fluid cell
		bool is_fluid(const unsigned i, const unsigned j, const unsigned k);
		//Return if the cell with center (i,j,k) is empty (not solid && not fluid)
		bool is_empty(const unsigned i, const unsigned j, const unsigned k);

		// Grid dimensions in #cells
		unsigned get_num_cells_x();
		unsigned get_num_cells_y();
		unsigned get_num_cells_z();
		unsigned get_num_cells();
		// Cell dimensions in metres
		double get_cell_sizex();
		double get_cell_sizey();
		double get_cell_sizez();
		// Get const-reference to A diagonal
		const std::vector< Triplet_t >& get_a_diag();
		
		//Get the weights for u in the mathematical point (i-1/2,j,k) 
		double get_weights_u(const unsigned i, const unsigned j, const unsigned k);
		//Get the weights for v in the mathematical point (i,j-1/2,k)
		double get_weights_v(const unsigned i, const unsigned j, const unsigned k);
		//Get the weights for v in the mathematical point (i,j,k-1/2)
		double get_weights_w(const unsigned i, const unsigned j, const unsigned k);
		
		//Getters for the trilinear interpolation of u, v, w, u*, v*, w* in the point x,y,z
		void assign_x(double x, int indices_x, int& ix0, int& ix1, double& x0, double& x1);
		void assign_y(double y, int indices_y, int& iy0, int& iy1, double& y0, double& y1);
		void assign_z(double z, int indices_z, int& iz0, int& iz1, double& z0, double& z1);
		void assign_d(double& xd, double& yd, double& zd, double x0, double x1, double y0, double y1, double z0, double z1, double x, double y, double z);
		double get_interp_u(double x, double y, double z, bool use_u_star);	// if use_u_star is true, then is the interpolation of u*!
		double get_interp_v(double x, double y, double z, bool use_v_star); // if use_v_star is true, then is the interpolation of v*!
		double get_interp_w(double x, double y, double z, bool use_w_star);	// if use_w_star is true, then is the interpolation of w*!	


		//SETTERS
		//Set the x-velocity in the mathematical point (i-1/2,j,k)
		void set_u(const unsigned i, const unsigned j, const unsigned k, double value); 
		//Set the y-velocity in the mathematical point (i,j-1/2,k)
		void set_v(const unsigned i, const unsigned j, const unsigned k, double value);
		//Set the z-velocity in the mathematical point (i,j,k-1/2)
		void set_w(const unsigned i, const unsigned j, const unsigned k, double value);
		// Copy velocity field to temporary copy arrays
		void set_uvw_star();

		//Set the pressure in the mathematical point(i,j,k)
		void set_pressure(const unsigned i, const unsigned j, const unsigned k, double value);
		// Reset all pressures
		void set_pressure(const Eigen::VectorXd& p);
		
		//Set the weights for u in the mathematical point (i-1/2,j,k) 
		void set_weights_u(const unsigned i, const unsigned j, const unsigned k, double value);
		//Set the weights for v in the mathematical point (i,j-1/2,k)
		void set_weights_v(const unsigned i, const unsigned j, const unsigned k, double value);
		//Set the weights for w in the mathematical point (i,j-1/2,k)
		void set_weights_w(const unsigned i, const unsigned j, const unsigned k, double value);
				
		// Set all grid velocities to zero for particle-to-grid transfer
		void set_velocities_to_zero();
		
		// Set all weights to zero for particle-to-grid transfer
		void set_weights_to_zero();
		
		//Set the cell with center (i,j,k) as a solid cell
		void set_solid(const unsigned i, const unsigned j, const unsigned k);
		//Set the cell with center (i,j,k) as a fluid cell
		void set_fluid(const unsigned i, const unsigned j, const unsigned k);
		// Reset all cells to not contain any fluid
		void reset_fluid();

		//USEFUL FUNCTIONS
		//Return a pair with the grid-coordinates (i,j,k) given a spatial coordinate(x,y,z)
		Eigen::Vector3d index_from_coord(const double x, const double y, const double z);
};
#endif
