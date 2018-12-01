#include "Mac2d.h"

// Constructor

Mac2d::Mac2d(const unsigned n, const unsigned m, const double dx, const double dy)
	: N_(n), M_(m), sizex_(dx), sizey_(dy), cell_sizex_(sizex_/(1.*N_)), cell_sizey_(sizey_/(1.*M_)){

	initArrays();

	// Initialize solid cells as a "box"
	//  Top & bottom
	//for (unsigned i = 0; i < N_; i++) {
	//	psolid_[i] = true;
	//	psolid_[i + (M_ - 1)*N_] = true;
	//}
	//// Sides
	//for (unsigned i = 1; i < M_ - 1; i++) {
	//	psolid_[N_*i] = true;
	//	for (unsigned j = 1; j < N_ - 1; j++) {
	//		psolid_[j + N_*i] = false;
	//	}
	//	psolid_[N_-1 + N_*i] = true;
	//}

	//Initialization of the diagonal of A
	initAdiag();
}


void Mac2d::initArrays() {
	ppressure_ = new double[N_*M_];
	std::fill(ppressure_, ppressure_+N_*M_, 0.);

	pu_ = new double[(N_+1)*M_];
	std::fill(pu_, pu_+(N_+1)*M_, 0.);

	pu_star_ = new double[(N_+1)*M_];
	std::fill(pu_star_, pu_star_+(N_+1)*M_, 0.);

	pv_ = new double[N_*(M_+1)];
	std::fill(pv_, pv_+N_*(M_+1), 0.);

	pv_star_ = new double[N_*(M_+1)];
	std::fill(pv_star_, pv_star_+N_*(M_+1), 0.);

	psolid_ = new bool[N_*M_];
	std::fill(psolid_, psolid_+N_*M_, 0.);

	pfluid_ = new bool[N_*M_];
	std::fill(pfluid_, pfluid_+N_*M_, 0.);

	pweights_u_ = new double[(N_+1)*M_];
	std::fill(pweights_u_, pweights_u_+(N_+1)*M_, 0.);

	pweights_v_ = new double[N_*(M_+1)];
	std::fill(pweights_v_, pweights_v_+N_*(M_+1), 0.);
}


void Mac2d::initAdiag() {
	A_diag_.clear();

	for(unsigned j = 0; j < M_; ++j){
		for(unsigned i = 0; i < N_; ++i){
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

//************************************************************************************
//***********************************GETTERS******************************************
//************************************************************************************
//1) getters for the layout properties of the grid (cell sizes and sizes of the grid);
//2) getters for the velocities (u, v, u*, v* and their interpolation);
//3) getters for the pressures;
//4) getters for the physical properties of the cells (solid, liquid, empty);
//5) getters for the weights for the particle to grid;
//6) getters for the diagonal of the pressure matrix A;
//7) getters for the indices of the cell to which a particle belong.

//1. Layout properties of the grid -------------------------------------
Eigen::Vector3d Mac2d::get_grid_size() const {
	return Eigen::Vector3d(sizex_, sizey_, 0.);
}
unsigned Mac2d::get_num_cells_x() {
	return N_;
}

unsigned Mac2d::get_num_cells_y() {
	return M_;
}

unsigned Mac2d::get_num_cells() {
	return M_*N_;
}

double Mac2d::get_cell_sizex() {
	return cell_sizex_;
}

double Mac2d::get_cell_sizey() {
	return cell_sizey_;
}

//2. Velocities --------------------------------------------------------
double Mac2d::get_u(const unsigned i, const unsigned j){
	if (i < (N_+1) && j < (M_) && "Index out of bounds!")
		return *(pu_ + (N_+1)*j + i);
	else 
		return 0;
	
}

double Mac2d::get_v(const unsigned i, const unsigned j){
	if (i < (N_) && j < (M_+1) && "Index out of bounds!")
		return *(pv_ + N_*j + i);
	else 
		return 0;
}

double Mac2d::get_u_star(const unsigned i, const unsigned j) {
	if (i < (N_+1) && j < (M_) && "Index out of bounds!")
		return *(pu_star_ + (N_+1)*j + i);
	else 
		return 0;
}

double Mac2d::get_v_star(const unsigned i, const unsigned j) {
	if (i < (N_) && j < (M_+1) && "Index out of bounds!")
		return *(pv_star_ + N_*j + i);
	else
		return 0;
}

double Mac2d::get_interp_u(double x, double y){
	Pair_t indices = index_from_coord(x,y);
	double x1, x2, y1, y2;
	int ix1, ix2, iy1, iy2;
	double u11, u12, u21, u22;
	
	if(y >= 0 && y <= sizey_ - cell_sizey_){
		//Update the u-velocity (bilinear interpolation)
		ix1 = indices.first;
		ix2 = ix1 + 1;
		if(y > indices.second * cell_sizey_){
			iy1 = indices.second;
			iy2 = iy1 + 1;
		}
		else{
			iy2 = indices.second;
			iy1 = iy2 - 1;
		}
		x1 = (ix1-0.5) * cell_sizex_;
		x2 = (ix2-0.5) * cell_sizex_;
		y1 = (iy1) * cell_sizey_;
		y2 = (iy2) * cell_sizey_;
		
		u11 = get_u(ix1,iy1);
		u12 = get_u(ix1,iy2);
		u21 = get_u(ix2,iy1);
		u22 = get_u(ix2,iy2);
		return 1/((x2 - x1)*(y2-y1))*(u11*(x2 - x)*(y2-y) 
				+ u21*(x - x1)*(y2-y) + u12*(x2 - x)*(y-y1) 
				+ u22*(x - x1)*(y-y1));
	}
	else if (y < 0){
		ix1 = indices.first;
		ix2 = ix1 + 1;
		x1 = (ix1-0.5) * cell_sizex_;
		x2 = (ix2-0.5) * cell_sizex_;
		u11 = get_u(ix1, 0);
		u21 = get_u(ix2, 0);
		return u11*(1-(x-x1)/(x2-x1)) + u21*((x-x1)/(x2-x1));
	}
	//else if (y > sizey_ - cell_sizey_){
	else {
		ix1 = indices.first;
		ix2 = ix1 + 1;
		x1 = (ix1-0.5) * cell_sizex_;
		x2 = (ix2-0.5) * cell_sizex_;
		u11 = get_u(ix1, M_-1);
		u21 = get_u(ix2, M_-1);
		return u11*(1- (x-x1)/(x2-x1)) + u21*((x-x1)/(x2-x1));
	}
}
		
double Mac2d::get_interp_v(double x, double y){
	Pair_t indices = index_from_coord(x,y);
	double x1, x2, y1, y2;
	int ix1, ix2, iy1, iy2;
	double v11, v12, v21, v22;
	
	if( x >= 0 && x <= sizex_ - cell_sizex_){
		//Update the v-velocity (bilinear interpolation)
		iy1 = indices.second;
		iy2 = iy1 + 1;
		if(x > indices.first * cell_sizex_){
			ix1 = indices.first;
			ix2 = ix1 + 1;
		}
		else{
			ix2 = indices.first;
			ix1 = ix2 - 1;
		}
		x1 = (ix1) * cell_sizex_;
		x2 = (ix2) * cell_sizex_;
		y1 = (iy1-0.5) * cell_sizey_;
		y2 = (iy2-0.5) * cell_sizey_;
		
		v11 = get_v(ix1,iy1);
		v12 = get_v(ix1,iy2);
		v21 = get_v(ix2,iy1);
		v22 = get_v(ix2,iy2);	
		return 1/((x2 - x1)*(y2-y1))*(v11*(x2 - x)*(y2-y) 
			   + v21*(x - x1)*(y2-y) + v12*(x2 - x)*(y-y1) 
			   + v22*(x - x1)*(y-y1));
	}
	else if (x < 0){
		iy1 = indices.second;
		iy2 = iy1 + 1;
		y1 = (iy1-0.5) * cell_sizex_;
		y2 = (iy2-0.5) * cell_sizex_;
		v11 = get_v(0, iy1);
		v21 = get_v(0, iy2);
		return v11*(1-(y-y1)/(y2-y1)) + v21*((y-y1)/(y2-y1));		
	}
	//else if (x > sizex_ - cell_sizex_){
	else {
		iy1 = indices.second;
		iy2 = iy1 + 1;
		y1 = (iy1-0.5) * cell_sizex_;
		y2 = (iy2-0.5) * cell_sizex_;
		v11 = get_v(N_-1, iy1);
		v21 = get_v(N_-1, iy2);
		return v11*(1-(y-y1)/(y2-y1)) + v21*((y-y1)/(y2-y1));
	}
	
}
	
double Mac2d::get_interp_u_star(double x, double y){
	Pair_t indices = index_from_coord(x,y);
	double x1, x2, y1, y2;
	int ix1, ix2, iy1, iy2;
	double u11, u12, u21, u22;
	
	if(y >= 0 && y <= sizey_ - cell_sizey_){
		//Update the u*-velocity (bilinear interpolation)
		ix1 = indices.first;
		ix2 = ix1 + 1;
		if(y > indices.second * cell_sizey_){
			iy1 = indices.second;
			iy2 = iy1 + 1;
		}
		else{
			iy2 = indices.second;
			iy1 = iy2 - 1;
		}
		x1 = (ix1-0.5) * cell_sizex_;
		x2 = (ix2-0.5) * cell_sizex_;
		y1 = (iy1) * cell_sizey_;
		y2 = (iy2) * cell_sizey_;
		
		u11 = get_u_star(ix1,iy1);
		u12 = get_u_star(ix1,iy2);
		u21 = get_u_star(ix2,iy1);
		u22 = get_u_star(ix2,iy2);
		return 1/((x2 - x1)*(y2-y1))*(u11*(x2 - x)*(y2-y) 
				+ u21*(x - x1)*(y2-y) + u12*(x2 - x)*(y-y1) 
				+ u22*(x - x1)*(y-y1));
	}
	else if (y < 0){
		ix1 = indices.first;
		ix2 = ix1 + 1;
		x1 = (ix1-0.5) * cell_sizex_;
		x2 = (ix2-0.5) * cell_sizex_;
		u11 = get_u_star(ix1, 0);
		u21 = get_u_star(ix2, 0);
		return u11*(1- (x-x1)/(x2-x1)) + u21*((x-x1)/(x2-x1));
	}
	else{
		ix1 = indices.first;
		ix2 = ix1 + 1;
		x1 = (ix1-0.5) * cell_sizex_;
		x2 = (ix2-0.5) * cell_sizex_;
		u11 = get_u_star(ix1, M_-1);
		u21 = get_u_star(ix2, M_-1);
		return u11*(1- (x-x1)/(x2-x1)) + u21*((x-x1)/(x2-x1));
	}
}
	
double Mac2d::get_interp_v_star(double x, double y){
	Pair_t indices = index_from_coord(x,y);
	double x1, x2, y1, y2;
	int ix1, ix2, iy1, iy2;
	double v11, v12, v21, v22;
	
	if( x >= 0 && x <= sizex_ - cell_sizex_){
		//Update the v*-velocity (bilinear interpolation)
		iy1 = indices.second;
		iy2 = iy1 + 1;
		if(x > (indices.first) * cell_sizex_){
			ix1 = indices.first;
			ix2 = ix1 + 1;
		}
		else{
			ix2 = indices.first;
			ix1 = ix2 - 1;
		}
		x1 = ix1 * cell_sizex_;
		x2 = ix2 * cell_sizex_;
		y1 = (iy1-0.5) * cell_sizey_;
		y2 = (iy2-0.5) * cell_sizey_;
		
		v11 = get_v_star(ix1,iy1);
		v12 = get_v_star(ix1,iy2);
		v21 = get_v_star(ix2,iy1);
		v22 = get_v_star(ix2,iy2);
		return 1/((x2 - x1)*(y2-y1))*(v11*(x2 - x)*(y2-y) 
			   + v21*(x - x1)*(y2-y) + v12*(x2 - x)*(y-y1) 
			   + v22*(x - x1)*(y-y1));
	}
	else if (x < 0){
		iy1 = indices.second;
		iy2 = iy1 + 1;
		y1 = (iy1-0.5) * cell_sizex_;
		y2 = (iy2-0.5) * cell_sizex_;
		v11 = get_v_star(0, iy1);
		v21 = get_v_star(0, iy2);
		return v11*(1-(y-y1)/(y2-y1)) + v21*((y-y1)/(y2-y1));		
	}
	else{
		iy1 = indices.second;
		iy2 = iy1 + 1;
		y1 = (iy1-0.5) * cell_sizex_;
		y2 = (iy2-0.5) * cell_sizex_;
		v11 = get_v_star(N_-1, iy1);
		v21 = get_v_star(N_-1, iy2);
		return v11*(1-(y-y1)/(y2-y1)) + v21*((y-y1)/(y2-y1));
	}
}

//3. Pressures ---------------------------------------------------------
double Mac2d::get_pressure(const unsigned i, const unsigned j){
	if (i < N_ && j < M_ && "Index out of bounds!")
		return *(ppressure_ + N_*j + i);
	else 
		return 0;
}

//4. Physical properties -----------------------------------------------
bool Mac2d::is_solid(const unsigned i, const unsigned j){
	if (i < N_ && j < M_ && "Index out of bounds!")
		return *(psolid_ + N_*j + i);
	else 
		return 0;
}

bool Mac2d::is_fluid(const unsigned i, const unsigned j){
	if (i < N_ && j < M_ && "Index out of bounds!")
		return *(pfluid_ + N_*j + i);
	else
		return 0;
}

bool Mac2d::is_empty(const unsigned i, const unsigned j){
	if (i < N_ && j < M_ && "Index out of bounds!")
		return (!is_fluid(i,j) && !is_solid(i,j));
	else 
		return 0;
}

//5. Weights -----------------------------------------------------------
double Mac2d::get_weights_u(const unsigned i, const unsigned j){
	if (i < (N_+1) && j < (M_) && "Index out of bounds!")
		return *(pweights_u_ + (N_+1)*j + i);
	else 
		return 0;
}

double Mac2d::get_weights_v(const unsigned i, const unsigned j){
	if (i < (N_) && j < (M_+1) && "Index out of bounds!")
		return *(pweights_v_ + N_*j + i);
	else
		return 0;
}

//6. Diagonal of A -----------------------------------------------------
const std::vector< Mac2d::Triplet_t >& Mac2d::get_a_diag() {
	return A_diag_;
}

//7. Indices from coordinate -------------------------------------------
Mac2d::Pair_t Mac2d::index_from_coord(const double x, const double y){
	assert((x < sizex_ - 0.5*cell_sizex_ && y < sizey_ - 0.5*cell_sizey_ && x > -0.5*cell_sizex_ && y > -0.5*cell_sizey_)
			&& "Attention: out of the grid!");
	return Pair_t(int(x/cell_sizex_ + 0.5), int(y/cell_sizey_ + 0.5));
}



//************************************************************************************
//***********************************SETTERS******************************************
//************************************************************************************
//1) setters for the velocities (u, v, u*, v*);
//2) setters for the pressures;
//3) setters for the physical properties of the cells (solid, liquid, empty);
//4) setters for the weights for the particle to grid;

//1. Velocities --------------------------------------------------------
void Mac2d::set_u(const unsigned i, const unsigned j, double value){
	if (i < (N_+1) && j < (M_) && "Index out of bounds!")
		*(pu_ + (N_+1)*j + i) = value;
}

void Mac2d::set_v(const unsigned i, const unsigned j, double value){
	if (i < (N_) && j < (M_+1) && "Index out of bounds!")
		*(pv_ + N_*j + i) = value;
}

void Mac2d::set_uv_star() {
	std::copy(pu_, pu_ + (N_+1)*M_, pu_star_);
	std::copy(pv_, pv_ + N_*(M_+1), pv_star_);
}

void Mac2d::set_velocities_to_zero(){
	std::fill(pu_, pu_ + (N_+1)*M_, 0);
	std::fill(pv_, pv_ + (M_+1)*N_, 0);
}

//2. Pressures ---------------------------------------------------------
void Mac2d::set_pressure(const unsigned i, const unsigned j, double value){
	if (i < N_ && j < M_ && "Index out of bounds!")
		*(ppressure_ + N_*j + i) = value;
}

void Mac2d::set_pressure(const Eigen::VectorXd& p) {
	std::copy(p.data(), p.data()+p.size(), ppressure_);
}

//3. Physical properties -----------------------------------------------
void Mac2d::set_solid(const unsigned i, const unsigned j){
	if (i < N_ && j < M_ && "Index out of bounds!")
		*(psolid_ + N_*j + i) = true;
}

void Mac2d::set_fluid(const unsigned i, const unsigned j){
	if (i < N_ && j < M_ && "Index out of bounds!")
		*(pfluid_ + N_*j + i) = true;
}

void Mac2d::reset_fluid() {
	std::fill(pfluid_, pfluid_ + get_num_cells(), false);
}

//4. Weights -----------------------------------------------------------
void Mac2d::set_weights_u(const unsigned i, const unsigned j, double value){
	if (i < (N_+1) && j < (M_) && "Index out of bounds!")
		*(pweights_u_ + (N_+1)*j + i) = value;
}

void Mac2d::set_weights_v(const unsigned i, const unsigned j, double value){
	if (i < (N_) && j < (M_+1) && "Index out of bounds!")
		*(pweights_v_ + N_*j + i) = value;
}

void Mac2d::set_weights_to_zero(){
	std::fill(pweights_u_, pweights_u_ + (N_+1)*M_, 0);
	std::fill(pweights_v_, pweights_v_ + (M_+1)*N_, 0);
}

