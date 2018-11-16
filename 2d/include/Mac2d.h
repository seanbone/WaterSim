#ifndef Mac2d_h
#define Mac2d_h

#include <Eigen/Sparse>	//used for the Eigen sparse matrix
#include <utility>		//used for std::pair
#include <iostream> 	//used for input and output
#include <iomanip>		//used for input and output

using Triplet_t = Eigen::Triplet<double>; 
using Pair_t = std::pair<int, int>; //To access the component of the Pair use the methods "first" and "second"

class Mac2d{
	private:
		int N_; 				//number of cells in x
		int M_; 				//number of cells in y
		double sizex_; 			//size (in meter) of the grid in x-direction
		double sizey_; 			//size (in meter) of the grid in y-direction
		double cell_sizex_; 	//size (in meter) of one cell in x-direction
		double cell_sizey_; 	//size (in meter) of one cell in x-direction
		double* ppressure_;		//pointer to array for the pressure
		double* pu_;			//pointer to array for the velocities in x-direction
		double* pv_;			//pointer to array for the velocities in y-direction
		Triplet_t* pA_diag_;	//sparse matrix which contains the data for the pressure equations
		bool* psolid_; 			//pointer to array for specifing if a cell is solid (1) or not(0)
		
	public:
		Mac2d() {}															//Default constructur
		Mac2d(const int n, const int m, const double dx, const double dy);	//Constructor with the number of cells (in both directions) and the length of the cells (in both directions)
		double get_u(const int i, const int j);								//Get the x-velocity in the mathematical point (i-1/2,j) 
		double get_v(const int i, const int j);								//Get the y-velocity in the mathematical point (i,j-1/2)
		double get_pressure(const int i, const int j);						//Get the pressure in the mathematical point(i,j)
		bool is_solid(const int i, const int j);							//Return if the cell with center (i,j) is a solid cell
		Pair_t index_from_coord(const double x, const double y);			//Return a pair with the grid-coordinates (i,j) given a spatial coordinate(x,y)
		void set_u(const int i, const int j);								//Set the x-velocity in the mathematical point (i-1/2,j) 
		void set_v(const int i, const int j);								//Set the y-velocity in the mathematical point (i,j-1/2)
		void set_pressure(const int i, const int j);						//Set the pressure in the mathematical point(i,j)
		void set_solid(const int i, const int j);							//Set the cell with center (i,j) as a solid cell
		std::ostream& operator<< (std::ostream& out, double* list);			//Operator<< overloading to print an array pointed by list
};
#endif
