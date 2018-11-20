#ifndef MAC2D_H
#define MAC2D_H

#include <Eigen/Sparse>	//used for the Eigen sparse matrix
#include <utility>		//used for std::pair
#include <iostream> 	//used for input and output
#include <iomanip>		//used for input and output
#include <Eigen/Sparse>	//used for the matrix A
#include <vector>		//used for std::vector
#include <algorithm>	// std::fill
#include <cassert>

class Mac2d{
	public:
		using Triplet_t = Eigen::Triplet<double>; 
		using Pair_t = std::pair<int, int>; //To access the component of the Pair use the methods "first" and "second"

	private:
		//PARAMETER for the grid
		//number of cells in x
		const int N_; 
		//number of cells in y
		const int M_;
		//size (in meter) of the grid in x-direction
		const double sizex_;
		//size (in meter) of the grid in y-direction
		const double sizey_;
		//size (in meter) of one cell in x-direction
		const double cell_sizex_; 
		//size (in meter) of one cell in x-direction
		const double cell_sizey_;
		
		//PHYSICAL VALUES objects
		//pointer to array for the pressure
		double* ppressure_;
		//pointer to array for the velocities in x-direction
		double* pu_;
		//pointer to array for the velocities in y-direction
		double* pv_;
		// Temporary copy of velocity field -> u^*, v^*
		double* pu_star_;
		double* pv_star_;
		//pointer to array for specifing if a cell is solid (1) or not(0)
		bool* psolid_;
		//pointer to array for specifing if a cell contains fluid (1) or not(0)
		bool* pfluid_;
		//pointer to a std::vector which contians the triplets for
		//   the diagonal of the matrix A, used to solve the pressures
		std::vector<Triplet_t> A_diag_;
		
		// Weights for particle-to-grid for u
		double* pweights_u_;
		
		// Weights for particle-to-grid for v
		double* pweights_v_;

	public:
		//CONSTRUCTORS
		//Default constructor
		Mac2d()
			: N_(0), M_(0), sizex_(0), sizey_(0), cell_sizex_(0), cell_sizey_(0){}
		
		//Constructor with the number of cells (in both directions) and the length of the box (in metres)
		Mac2d(const int n, const int m, const double dx, const double dy)
			: N_(n), M_(m), sizex_(dx), sizey_(dy), cell_sizex_(sizex_/(1.*N_)), cell_sizey_(sizey_/(1.*M_)){
			ppressure_ = new double[N_*M_];
			pu_ = new double[(N_+1)*M_];
			pu_star_ = new double[(N_+1)*M_];
			pv_ = new double[N_*(M_+1)];
			pv_star_ = new double[N_*(M_+1)];
			psolid_ = new bool[N_*M_];
			pfluid_ = new bool[N_*M_];
			pweights_u_ = new double[(N_+1)*M_];
			pweights_v_ = new double[N_*(M_+1)];

			// Initialize solid cells as a "box"
			//  Top & bottom
			for (int i = 0; i < N_; i++) {
				psolid_[i] = true;
				psolid_[i + (M_ - 1)*N_] = true;
			}
			// Sides
			for (int i = 1; i < M_ - 1; i++) {
				psolid_[N_*i] = true;
				for (int j = 1; j < N_ - 1; j++) {
					psolid_[j + N_*i] = false;
				}
				psolid_[N_-1 + N_*i] = true;
			}

			//Initialization of the diagonal of A
			for(int j = 0; j< M_; ++j){
				for(int i = 0; i < N_; ++i){
					int index = N_ * j + i;
					int count = 0;
					if (i == 0){
						if (j == 0){
							count = !is_solid(i+1,j) + !is_solid(i,j+1);
						}
						else if (j == M_ - 1){
							count = !is_solid(i+1,j) + !is_solid(i,j-1);
						}
						else{
							count = !is_solid(i+1,j) + !is_solid(i,j-1) + !is_solid(i,j+1);
						}
					}
					else if (i == N_ - 1){
						if (j == 0){
							count = !is_solid(i-1,j) + !is_solid(i,j+1);
						}
						else if (j == M_ - 1){
							count = !is_solid(i-1,j) + !is_solid(i,j-1);
						}
						else{
							count = !is_solid(i-1,j) + !is_solid(i,j-1) + !is_solid(i,j+1);
						}
					}
					else {
						if (j == 0){
							count = !is_solid(i+1,j) + !is_solid(i-1,j) + !is_solid(i,j+1);
						}
						else if (j == M_ - 1){
							count = !is_solid(i+1,j) + !is_solid(i-1,j) + !is_solid(i,j-1);
						}
						else{
							count = !is_solid(i-1,j) + !is_solid(i+1,j) + !is_solid(i,j-1) + !is_solid(i,j+1);
						}
					}
					A_diag_.push_back(Triplet_t(index, index, count));
				}
			}
		}
		
		//Destructor
		~Mac2d(){
			delete[] ppressure_;
			delete[] pu_;
			delete[] pv_;
			delete[] pu_star_;
			delete[] pv_star_;
			delete[] psolid_;
			delete[] pfluid_;
		}
		
		//GETS
		//Get the x-velocity in the mathematical point (i-1/2,j) 
		double get_u(const int i, const int j);
		//Get the y-velocity in the mathematical point (i,j-1/2)
		double get_v(const int i, const int j);
		// Equivalent for intermediate field
		double get_u_star(const int i, const int j);
		double get_v_star(const int i, const int j);
		//Get the pressure in the mathematical point(i,j)
		double get_pressure(const int i, const int j);
		//Return if the cell with center (i,j) is a solid cell
		bool is_solid(const int i, const int j);
		//Return if the cell with center (i,j) is a fluid cell
		bool is_fluid(const int i, const int j);
		//Return if the cell with center (i,j) is empty (not solid && not fluid)
		bool is_empty(const int i, const int j);

		// Grid dimensions in #cells
		unsigned get_num_cells_x();
		unsigned get_num_cells_y();
		unsigned get_num_cells();
		// Cell dimensions in metres
		double get_cell_sizex();
		double get_cell_sizey();
		// Get const-reference to A diagonal
		const std::vector< Triplet_t >& get_a_diag();
		
		//Get the weights for u in the mathematical point (i-1/2,j) 
		double get_weights_u(const int i, const int j);
		//Get the weights for v in the mathematical point (i,j-1/2)
		double get_weights_v(const int i, const int j);

		//SETTERS
		//Set the x-velocity in the mathematical point (i-1/2,j)
		void set_u(const int i, const int j, double value); 
		//Set the y-velocity in the mathematical point (i,j-1/2)
		void set_v(const int i, const int j, double value);
		// Copy velocity field to temporary copy arrays
		void set_uv_star();

		//Set the pressure in the mathematical point(i,j)
		void set_pressure(const int i, const int j, double value);
		// Reset all pressures
		void set_pressure(const Eigen::VectorXd& p);
		
		//Set the weights for u in the mathematical point (i-1/2,j) 
		void set_weights_u(const int i, const int j, double value);
		//Set the weights for v in the mathematical point (i,j-1/2)
		void set_weights_v(const int i, const int j, double value);
		
		// Set all grid velocities to zero for particle-to-grid transfer
		void set_velocities_to_zero();
		
		// Set all weights to zero for particle-to-grid transfer
		void set_weights_to_zero();
		
		//Set the cell with center (i,j) as a solid cell
		void set_solid(const int i, const int j);
		//Set the cell with center (i,j) as a fluid cell
		void set_fluid(const int i, const int j);
		// Reset all cells to not contain any fluid
		void reset_fluid();

		//USEFUL FUNCTIONS
		//Return a pair with the grid-coordinates (i,j) given a spatial coordinate(x,y)
		Pair_t index_from_coord(const double x, const double y);
};
#endif
