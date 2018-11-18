#ifndef MAC2D_H
#define MAC2D_H

#include <Eigen/Sparse>	//used for the Eigen sparse matrix
#include <utility>		//used for std::pair
#include <iostream> 	//used for input and output
#include <iomanip>		//used for input and output
#include <Eigen/Sparse>	//used for the matrix A
#include <vector>		//used for std::vector

using Triplet_t = Eigen::Triplet<double>; 
using Pair_t = std::pair<int, int>; //To access the component of the Pair use the methods "first" and "second"

class Mac2d{
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
		//pointer to array for specifing if a cell is solid (1) or not(0)
		bool* psolid_;
		//pointer to a std::vector which contians the triplets for
		//   the diagonal of the matrix A, used to solve the pressures
		std::vector<Triplet_t> A_diag_; 
		

	public:
		//CONSTRUCTORS
		//Default constructor
		Mac2d()
			: N_(0), M_(0), sizex_(0), sizey_(0), cell_sizex_(0), cell_sizey_(0){}
		
		//Constructor with the number of cells (in both directions) and the length of the cells (in both directions)
		Mac2d(const int n, const int m, const double dx, const double dy)
			: N_(n), M_(m), sizex_(dx), sizey_(dy), cell_sizex_(sizex_/(1.*N_)), cell_sizey_(sizex_/(1.*N_)){
			ppressure_ = new double[N_*M_];
			pu_ = new double[(N_+1)*M_];
			pv_ = new double[N_*(M_+1)];
			psolid_ = new bool[N_*M_];
			
			//Initialization of the diagonal of A
			for(int i = 0; i < N_; ++i){
				for(int j = 0; j < M_; ++j){
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
				}
				A_diag_.push_back(Triplet_t(index, index, count));
			}
		}
		
		//Destructor
		~Mac2d(){
			delete ppressure_;
			delete pu_;
			delete pv_;
			delete psolid_;
		}
		
		//GETS
		//Get the x-velocity in the mathematical point (i-1/2,j) 
		double get_u(const int i, const int j);
		//Get the y-velocity in the mathematical point (i,j-1/2)
		double get_v(const int i, const int j);
		//Get the pressure in the mathematical point(i,j)
		double get_pressure(const int i, const int j);
		//Return if the cell with center (i,j) is a solid cell
		bool is_solid(const int i, const int j);
		
		//SETTERS
		//Set the x-velocity in the mathematical point (i-1/2,j)
		void set_u(const int i, const int j, double value); 
		//Set the y-velocity in the mathematical point (i,j-1/2)
		void set_v(const int i, const int j, double value);
		//Set the pressure in the mathematical point(i,j)
		void set_pressure(const int i, const int j, double value);
		
		//Set the cell with center (i,j) as a solid cell
		void set_solid(const int i, const int j);

		//USEFUL FUNCTIONS
		//Return a pair with the grid-coordinates (i,j) given a spatial coordinate(x,y)
		Pair_t index_from_coord(const double x, const double y);
};
#endif
