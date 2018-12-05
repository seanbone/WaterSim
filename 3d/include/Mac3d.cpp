#include "Mac3d.h"

// Constructor

Mac3d::Mac3d(const unsigned n, const unsigned m, const unsigned l, const double dx, const double dy, const double dz)
	: N_(n), M_(m), L_(l), sizex_(dx), sizey_(dy), sizez_(dz), cell_sizex_(sizex_/(1.*N_)), cell_sizey_(sizey_/(1.*M_)), cell_sizez_(sizez_/(1.*L_)){

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


void Mac3d::initArrays() {
	ppressure_ = new double[N_*M_*L_];
	std::fill(ppressure_, ppressure_+N_*M_*L_, 0.);

	pu_ = new double[(N_+1)*M_*L_];
	std::fill(pu_, pu_+(N_+1)*M_*L_, 0.);

	pu_star_ = new double[(N_+1)*M_*L_];
	std::fill(pu_star_, pu_star_+(N_+1)*M_*L_, 0.);

	pv_ = new double[N_*(M_+1)*L_];
	std::fill(pv_, pv_+N_*(M_+1)*L_, 0.);

	pv_star_ = new double[N_*(M_+1)*L_];
	std::fill(pv_star_, pv_star_+N_*(M_+1)*L_, 0.);
	
	pw_ = new double[N_*M_*(L_+1)];
	std::fill(pw_, pw_+N_*M_*(L_+1), 0.);

	pw_star_ = new double[N_*M_*(L_+1)];
	std::fill(pw_star_, pw_star_+N_*M_*(L_+1), 0.);

	psolid_ = new bool[N_*M_*L_];
	std::fill(psolid_, psolid_+N_*M_*L_, 0.);

	pfluid_ = new bool[N_*M_*L_];
	std::fill(pfluid_, pfluid_+N_*M_*L_, 0.);

	pweights_u_ = new double[(N_+1)*M_*L_];
	std::fill(pweights_u_, pweights_u_+(N_+1)*M_*L_, 0.);

	pweights_v_ = new double[N_*(M_+1)*L_];
	std::fill(pweights_v_, pweights_v_+N_*(M_+1)*L_, 0.);
	
	pweights_w_ = new double[N_*M_*(L_+1)];
	std::fill(pweights_w_, pweights_w_+N_*M_*(L_+1), 0.);
}

void Mac3d::initAdiag() {
	A_diag_.clear();

	for(unsigned k = 0; k < L_; ++k){
		for(unsigned j = 0; j < M_; ++j){
			for(unsigned i = 0; i < N_; ++i){
				int index = N_ * j + i + N_*M_*k;
				int count = 0;
				if (i == 0){
					if (j == 0){
						if (k == 0)
							count = !is_solid(i+1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i+1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i+1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
					else if (j == M_-1){
						if(k == 0)
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
					else{
						if (k == 0)
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
				}
				else if (i == N_ - 1){
					if (j == 0){
						if (k == 0)
							count = !is_solid(i-1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i-1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i-1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
					else if (j == M_-1){
						if(k == 0)
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
					else{
						if(k == 0)
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
				}
				else {
					if (j == 0){
						if(k == 0)
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
					else if (j == M_ - 1){
						if(k == 0)
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
					else{
						if(k == 0)
							count = !is_solid(i-1,j,k) + !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						else if (k == L_-1)
							count = !is_solid(i-1,j,k) + !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						else
							count = !is_solid(i-1,j,k) + !is_solid(i+1,j,k) + !is_solid(i,j-1,k) + !is_solid(i,j+1,k) + !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
					}
				}
				A_diag_.push_back(Triplet_t(index, index, count));
			}
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
Eigen::Vector3d Mac3d::get_grid_size() const {
	return Eigen::Vector3d(sizex_, sizey_, sizez_);
}

unsigned Mac3d::get_num_cells_x() {
	return N_;
}

unsigned Mac3d::get_num_cells_y() {
	return M_;
}

unsigned Mac3d::get_num_cells_z() {
	return L_;
}

unsigned Mac3d::get_num_cells() {
	return M_*N_*L_;
}

double Mac3d::get_cell_sizex() {
	return cell_sizex_;
}

double Mac3d::get_cell_sizey() {
	return cell_sizey_;
}

double Mac3d::get_cell_sizez() {
	return cell_sizez_;
}

//2. Velocities --------------------------------------------------------
double Mac3d::get_u(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < (N_+1) && j < M_ && k < L_ && "Index out of bounds!")
		return *(pu_ + (N_+1)*j + i + (N_+1)*M_*k);
	else 
		return 0;
	
}

double Mac3d::get_v(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < (M_+1) && k < L_ && "Index out of bounds!")
		return *(pv_ + N_*j + i + N_*(M_+1)*k);
	else 
		return 0;
}

double Mac3d::get_w(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < (L_+1) && "Index out of bounds!")
		return *(pw_ + N_*j + i + N_*M_*k);
	else 
		return 0;
}

double Mac3d::get_u_star(const unsigned i, const unsigned j, const unsigned k = 0) {
	if (i < (N_+1) && j < M_ && k < L_ && "Index out of bounds!")
		return *(pu_star_ + (N_+1)*j + i + (N_+1)*M_*k);
	else 
		return 0;
}

double Mac3d::get_v_star(const unsigned i, const unsigned j, const unsigned k = 0) {
	if (i < N_ && j < (M_+1) && k < L_ && "Index out of bounds!")
		return *(pv_star_ + N_*j + i + N_*(M_+1)*k);
	else
		return 0;
}

double Mac3d::get_w_star(const unsigned i, const unsigned j, const unsigned k = 0) {
	if (i < N_ && j < M_ && k < (L_+1) && "Index out of bounds!")
		return *(pw_star_ + N_*j + i + N_*M_*k);
	else
		return 0;
}

//TO DO!!!
void Mac3d::assign_y(double y, int indices_y, int& iy0, int& iy1, double& y0, double& y1){
	if(y > indices_y * cell_sizey_){
		iy0 = indices_y;
		iy1 = iy0 + 1;
	}
	else{
		iy1 = indices_y;
		iy0 = iy1 - 1;
	}
	y0 = iy0*cell_sizey_;
	y1 = iy1*cell_sizey_;
	return;
}

void Mac3d::assign_z(double z, int indices_z, int& iz0, int& iz1, double& z0, double& z1){
	if(z > indices_z * cell_sizez_){
		iz0 = indices_z;
		iz1 = iz0 + 1;
	}
	else{
		iz1 = indices_z;
		iz0 = iz1 - 1;
	}
	z0 = iz0*cell_sizez_;
	z1 = iz1*cell_sizez_;
	return;
}
void Mac3d::assign_d(double& xd, double& yd, double& zd, double x0, double x1, 
			  double y0, double y1, double z0, double z1, double x, double y, double z){
	xd = (x-x0)/(x1-x0);
	yd = (y-y0)/(y1-y0);
	zd = (z-z0)/(z1-z0);
	return;
}

double Mac3d::get_interp_u(double x, double y, double z, bool use_u_star = false){
	//u or u*?
	double (Mac3d::*get_vel)(unsigned, unsigned, unsigned);
	if (use_u_star)
		get_vel = &Mac3d::get_u_star;
	else
		get_vel = &Mac3d::get_u;
	
	Eigen::Vector3d indices =  index_from_coord(x,y,z);
	double x0, x1, y0, y1, z0, z1;
	int ix0, ix1, iy0, iy1, iz0, iz1;
	double u000, u100, u110, u010, u001, u101, u111, u011;
	double xd, yd, zd;
	double u00, u10, u01, u11;
	double u0, u1;
	
	ix0 = indices[0];
	ix1 = ix0 + 1;
	x0 = (ix0-0.5)*cell_sizex_;
	x1 = (ix1-0.5)*cell_sizex_;
	
	if(y >= 0 && y <= sizey_ - cell_sizey_){
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			
			u000 = (this->*get_vel)(ix0, iy0, iz0);
			u100 = (this->*get_vel)(ix1, iy0, iz0);
			u010 = (this->*get_vel)(ix0, iy1, iz0);
			u001 = (this->*get_vel)(ix0, iy0, iz1);
			u110 = (this->*get_vel)(ix1, iy1, iz0);
			u101 = (this->*get_vel)(ix1, iy0, iz1);
			u011 = (this->*get_vel)(ix0, iy1, iz1);
			u111 = (this->*get_vel)(ix1, iy1, iz1);
			
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			
			u00 = u000*(1-xd) + u100*xd;
			u01 = u001*(1-xd) + u101*xd;
			u10 = u010*(1-xd) + u110*xd;
			u11 = u011*(1-xd) + u111*xd;
			
			u0 = u00*(1-yd) + u10*yd;
			u1 = u01*(1-yd) + u11*yd;
			
			return u0*(1-zd) + u1*zd;
		}
		else if (z < 0){
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			
			u00 = (this->*get_vel)(ix0,iy0, 0);
			u01 = (this->*get_vel)(ix0,iy1, 0);
			u10 = (this->*get_vel)(ix1,iy0, 0);
			u11 = (this->*get_vel)(ix1,iy1, 0);
			
			u0 = u00*(1-xd) + u10*xd;
			u1 = u01*(1-xd) + u11*xd;
			
			return u0*(1-yd) + u1*yd;
		}
		else{
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			
			u00 = (this->*get_vel)(ix0,iy0, L_-1);
			u01 = (this->*get_vel)(ix0,iy1, L_-1);
			u10 = (this->*get_vel)(ix1,iy0, L_-1);
			u11 = (this->*get_vel)(ix1,iy1, L_-1);
			
			u0 = u00*(1-xd) + u10*xd;
			u1 = u01*(1-xd) + u11*xd;
			
			return u0*(1-yd) + u1*yd;
		}
	}
	else if (y < 0){
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			
			u00 = (this->*get_vel)(ix0, 0, iz0);
			u01 = (this->*get_vel)(ix0, 0, iz1);
			u10 = (this->*get_vel)(ix1, 0, iz0);
			u11 = (this->*get_vel)(ix1, 0, iz1);
			
			u0 = u00*(1-xd) + u10*xd;
			u1 = u01*(1-xd) + u11*xd;
			
			return u0*(1-zd) + u1*zd;
		}
		else if (z < 0){
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			u0 = (this->*get_vel)(x0, 0, 0);
			u1 = (this->*get_vel)(x1, 0, 0);
			return u0*(1-xd) + u1*xd;				
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			u0 = (this->*get_vel)(x0, 0, L_-1);
			u1 = (this->*get_vel)(x1, 0, L_-1);
			return u0*(1-xd) + u1*xd;
		}
	}
	else{
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			
			u00 = (this->*get_vel)(ix0, M_-1, iz0);
			u01 = (this->*get_vel)(ix0, M_-1, iz1);
			u10 = (this->*get_vel)(ix1, M_-1, iz0);
			u11 = (this->*get_vel)(ix1, M_-1, iz1);
			
			u0 = u00*(1-xd) + u10*xd;
			u1 = u01*(1-xd) + u11*xd;
			
			return u0*(1-zd) + u1*zd;
		}
		else if (z < 0){
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			u0 = (this->*get_vel)(x0, M_-1, 0);
			u1 = (this->*get_vel)(x1, M_-1, 0);
			return u0*(1-xd) + u1*xd;
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			u0 = (this->*get_vel)(x0, M_-1, L_-1);
			u1 = (this->*get_vel)(x1, M_-1, L_-1);
			return u0*(1-xd) + u1*xd;
		}
	}
}
		
/*double Mac3d::get_interp_v(double x, double y){
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
	
double Mac3d::get_interp_u_star(double x, double y){
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
	
double Mac3d::get_interp_v_star(double x, double y){
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
}*/

//3. Pressures ---------------------------------------------------------
double Mac3d::get_pressure(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		return *(ppressure_ + N_*j + i + N_*M_*k);
	else 
		return 0;
}

//4. Physical properties -----------------------------------------------
bool Mac3d::is_solid(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		return *(psolid_ + N_*j + i + N_*M_*k);
	else 
		return 0;
}

bool Mac3d::is_fluid(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		return *(pfluid_ + N_*j + i + N_*M_*k);
	else
		return 0;
}

bool Mac3d::is_empty(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		return (!is_fluid(i,j,k) && !is_solid(i,j,k));
	else 
		return 0;
}

//5. Weights -----------------------------------------------------------
double Mac3d::get_weights_u(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < (N_+1) && j < M_ && k < L_ && "Index out of bounds!")
		return *(pweights_u_ + (N_+1)*j + i + (N_+1)*M_*k);
	else 
		return 0;
}

double Mac3d::get_weights_v(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < (M_+1) && k < L_ && "Index out of bounds!")
		return *(pweights_v_ + N_*j + i + N_*(M_+1)*k);
	else
		return 0;
}

double Mac3d::get_weights_w(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < (L_+1) && "Index out of bounds!")
		return *(pweights_w_ + N_*j + i + N_*M_*k);
	else
		return 0;
}
//6. Diagonal of A -----------------------------------------------------
const std::vector< Mac3d::Triplet_t >& Mac3d::get_a_diag() {
	return A_diag_;
}

//7. Indices from coordinate -------------------------------------------
Eigen::Vector3d Mac3d::index_from_coord(const double x, const double y, const double z){
	assert((x < sizex_ - 0.5*cell_sizex_ && y < sizey_ - 0.5*cell_sizey_ && z < sizez_ - 0.5*cell_sizez_
		    && x > -0.5*cell_sizex_ && y > -0.5*cell_sizey_ && z > -0.5*cell_sizez_)
			&& "Attention: out of the grid!");
	Eigen::Vector3d result(int(x/cell_sizex_ + 0.5), int(y/cell_sizey_ + 0.5), int(z/cell_sizez_ + 0.5));
	return result;
}



//************************************************************************************
//***********************************SETTERS******************************************
//************************************************************************************
//1) setters for the velocities (u, v, u*, v*);
//2) setters for the pressures;
//3) setters for the physical properties of the cells (solid, liquid, empty);
//4) setters for the weights for the particle to grid;

//1. Velocities --------------------------------------------------------
void Mac3d::set_u(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_+1) && j < M_ && k < L_ && "Index out of bounds!")
		*(pu_ + (N_+1)*j + i + (N_+1)*M_*k) = value;
}

void Mac3d::set_v(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < (M_+1) && k < L_ && "Index out of bounds!")
		*(pv_ + N_*j + i + N_*(M_+1)*k) = value;
}

void Mac3d::set_w(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < M_ && k < (L_+1) && "Index out of bounds!")
		*(pw_ + N_*j + i + N_*M_*k) = value;
}

void Mac3d::set_uvw_star() {
	std::copy(pu_, pu_ + (N_+1)*M_*L_, pu_star_);
	std::copy(pv_, pv_ + N_*(M_+1)*L_, pv_star_);
	std::copy(pw_, pw_ + N_*M_*(L_+1), pw_star_);
}

void Mac3d::set_velocities_to_zero(){
	std::fill(pu_, pu_ + (N_+1)*M_*L_, 0);
	std::fill(pv_, pv_ + N_*(M_+1)*L_, 0);
	std::fill(pw_, pw_ + N_*M_*(L_+1), 0);
}

//2. Pressures ---------------------------------------------------------
void Mac3d::set_pressure(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		*(ppressure_ + N_*j + i + N_*M_*k) = value;
}

void Mac3d::set_pressure(const Eigen::VectorXd& p) {
	std::copy(p.data(), p.data()+p.size(), ppressure_);
}

//3. Physical properties -----------------------------------------------
void Mac3d::set_solid(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		*(psolid_ + N_*j + i + N_*M_*k) = true;
}

void Mac3d::set_fluid(const unsigned i, const unsigned j, const unsigned k = 0){
	if (i < N_ && j < M_ && k < L_ && "Index out of bounds!")
		*(pfluid_ + N_*j + i + N_*M_*k) = true;
}


void Mac3d::reset_fluid() {
	std::fill(pfluid_, pfluid_ + get_num_cells(), false);
}

//4. Weights -----------------------------------------------------------
void Mac3d::set_weights_u(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_+1) && j < (M_) && k < L_ && "Index out of bounds!")
		*(pweights_u_ + (N_+1)*j + i + (N_+1)*M_*k) = value;
}

void Mac3d::set_weights_v(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_) && j < (M_+1) && k < L_ && "Index out of bounds!")
		*(pweights_v_ + N_*j + i + N_*(M_+1)*k) = value;
}

void Mac3d::set_weights_w(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_) && j < M_ && k < (L_+1) && "Index out of bounds!")
		*(pweights_w_ + N_*j + i + N_*M_*k) = value;
}

void Mac3d::set_weights_to_zero(){
	std::fill(pweights_u_, pweights_u_ + (N_+1)*M_*L_, 0);
	std::fill(pweights_v_, pweights_v_ + N_*(M_+1)*L_, 0);
	std::fill(pweights_w_, pweights_w_ + N_*M_*(L_+1), 0);
}

