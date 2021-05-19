#include "Mac3d.h"


Mac3d::Mac3d(const unsigned n, const unsigned m, const unsigned l, 
			 const double dx, const double dy, const double dz)
	: N_(n), M_(m), L_(l), sizex_(dx), sizey_(dy), sizez_(dz), 
	  cell_sizex_(sizex_/(1.*N_)), cell_sizey_(sizey_/(1.*M_)), 
	  cell_sizez_(sizez_/(1.*L_)), rcell_sizex_(1./cell_sizex_), 
	  rcell_sizey_(1./cell_sizey_), rcell_sizez_(1./cell_sizez_){
	
	//Initialization of the arrays 
	initArrays();

	//Initialization of the diagonal of A
	initAdiag();

	GRID_OFFSETS[GRID_P][0] = 0;
	GRID_OFFSETS[GRID_P][1] = 0;
	GRID_OFFSETS[GRID_P][2] = 0;

	GRID_OFFSETS[GRID_U][0] = -0.5 * cell_sizex_;
	GRID_OFFSETS[GRID_U][1] = 0;
	GRID_OFFSETS[GRID_U][2] = 0;

	GRID_OFFSETS[GRID_U_STAR][0] = -0.5 * cell_sizex_;
	GRID_OFFSETS[GRID_U_STAR][1] = 0;
	GRID_OFFSETS[GRID_U_STAR][2] = 0;

	GRID_OFFSETS[GRID_V][0] = 0;
	GRID_OFFSETS[GRID_V][1] = -0.5 * cell_sizey_;
	GRID_OFFSETS[GRID_V][2] = 0;

	GRID_OFFSETS[GRID_V_STAR][0] = 0;
	GRID_OFFSETS[GRID_V_STAR][1] = -0.5 * cell_sizey_;
	GRID_OFFSETS[GRID_V_STAR][2] = 0;

	GRID_OFFSETS[GRID_W][0] = 0;
	GRID_OFFSETS[GRID_W][1] = 0;
	GRID_OFFSETS[GRID_W][2] = -0.5 * cell_sizez_;

	GRID_OFFSETS[GRID_W_STAR][0] = 0;
	GRID_OFFSETS[GRID_W_STAR][1] = 0;
	GRID_OFFSETS[GRID_W_STAR][2] = -0.5 * cell_sizez_;
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
						if (k == 0){
							count = !is_solid(i+1,j,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i+1,j,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i+1,j,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
						}
					}
					else if (j == M_-1){
						if(k == 0){
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
						}
					}
					else{
						if (k == 0){
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i+1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k-1) 
									+ !is_solid(i,j,k+1);
						}
					}
				}
				else if (i == N_ - 1){
					if (j == 0){
						if (k == 0){
							count = !is_solid(i-1,j,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i-1,j,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i-1,j,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
						}
					}
					else if (j == M_-1){
						if(k == 0){
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
						}
					}
					else{
						if(k == 0){
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i-1,j,k) + !is_solid(i,j-1,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k-1) 
									+ !is_solid(i,j,k+1);
						}
					}
				}
				else {
					if (j == 0){
						if(k == 0){
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) 
									+ !is_solid(i,j+1,k) + !is_solid(i,j,k-1) 
									+ !is_solid(i,j,k+1);
						}
					}
					else if (j == M_ - 1){
						if(k == 0){
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) 
									+ !is_solid(i,j-1,k) + !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) 
									+ !is_solid(i,j-1,k) + !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i+1,j,k) + !is_solid(i-1,j,k) 
									+ !is_solid(i,j-1,k) + !is_solid(i,j,k-1) 
									+ !is_solid(i,j,k+1);
						}
					}
					else{
						if(k == 0){
							count = !is_solid(i-1,j,k) + !is_solid(i+1,j,k) 
									+ !is_solid(i,j-1,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k+1);
						}
						else if (k == L_-1){
							count = !is_solid(i-1,j,k) + !is_solid(i+1,j,k) 
									+ !is_solid(i,j-1,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k-1);
						}
						else{
							count = !is_solid(i-1,j,k) + !is_solid(i+1,j,k) 
									+ !is_solid(i,j-1,k) + !is_solid(i,j+1,k) 
									+ !is_solid(i,j,k-1) + !is_solid(i,j,k+1);
						}
					}
				}
				A_diag_.push_back(Triplet_t(index, index, count));
			}
		}
	}
}

/************************************************************************************
************************************GETTERS******************************************
*************************************************************************************
//1) getters for the layout properties of the grid (cell sizes and sizes of the grid);
//2) getters for the velocities (u, v, w, u*, v*, w* and their interpolation);
//3) getters for the pressures;
//4) getters for the physical properties of the cells (solid, liquid, empty);
//5) getters for the weights for the particle to grid;
//6) getters for the diagonal of the pressure matrix A;
//7) getters for the indices of the cell to which a particle belong.
*/

//1. Layout properties of the grid -------------------------------------
Eigen::Vector3d Mac3d::get_grid_size() const {
	return Eigen::Vector3d(sizex_, sizey_, sizez_);
}

unsigned Mac3d::get_num_cells_x() const {
	return N_;
}

unsigned Mac3d::get_num_cells_y() const {
	return M_;
}

unsigned Mac3d::get_num_cells_z() const {
	return L_;
}

unsigned Mac3d::get_num_cells() const {
	return M_*N_*L_;
}

double Mac3d::get_cell_sizex() const {
	return cell_sizex_;
}

double Mac3d::get_cell_sizey() const {
	return cell_sizey_;
}

double Mac3d::get_cell_sizez() const {
	return cell_sizez_;
}

//2. Velocities --------------------------------------------------------
double Mac3d::get_u(const unsigned i, const unsigned j, const unsigned k){
	if (i < (N_+1) && j < M_ && k < L_)
		return *(pu_ + (N_+1)*j + i + (N_+1)*M_*k);
	else{ 
		std::cout << "Calling get_u: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
	
}

double Mac3d::get_v(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < (M_+1) && k < L_) {
		//if (i == 0 && j == 1 && k == 0)
		//	std::cout << "-> " << N_ * j + i + N_ * (M_ + 1) * k << std::endl;
		return *(pv_ + N_ * j + i + N_ * (M_ + 1) * k);
	} else {
		std::cout << "Calling get_v: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

double Mac3d::get_w(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < (L_+1))
		return *(pw_ + N_*j + i + N_*M_*k);
	else{ 
		std::cout << "Calling get_w: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

double Mac3d::get_u_star(const unsigned i, const unsigned j, const unsigned k) {
	if (i < (N_+1) && j < M_ && k < L_)
		return *(pu_star_ + (N_+1)*j + i + (N_+1)*M_*k);
	else{ 
		std::cout << "Calling get_u_star: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

double Mac3d::get_v_star(const unsigned i, const unsigned j, const unsigned k) {
	if (i < N_ && j < (M_+1) && k < L_)
		return *(pv_star_ + N_*j + i + N_*(M_+1)*k);
	else{ 
		std::cout << "Calling get_v_star: Index out (" << i << ", " << j << ", " << k << ") of bounds!" << std::endl;
		return 0;
	}
}

double Mac3d::get_w_star(const unsigned i, const unsigned j, const unsigned k) {
	if (i < N_ && j < M_ && k < (L_+1))
		return *(pw_star_ + N_*j + i + N_*M_*k);
	else{ 
		std::cout << "Calling get_w_star: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

void Mac3d::assign_x(double x, int indices_x, int& ix0, int& ix1, double& x0, double& x1){
	if(x > indices_x * cell_sizex_ or  x == 0){
		ix0 = indices_x;
		ix1 = ix0 + 1;
	}
	else{
		ix1 = indices_x;
		ix0 = ix1 - 1;
	}
	x0 = ix0*cell_sizex_;
	x1 = ix1*cell_sizex_;
	return;
}

void Mac3d::assign_y(double y, int indices_y, int& iy0, int& iy1, double& y0, double& y1){
	if(y > indices_y * cell_sizey_ or  y == 0){
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
	if(z > indices_z * cell_sizez_ or  z == 0){
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

double Mac3d::get_interp_u(double x, double y, double z, const bool use_u_star){
	//Setting u or u*
	double (Mac3d::*get_vel)(unsigned, unsigned, unsigned);
	if (use_u_star)
		get_vel = &Mac3d::get_u_star;
	else
		get_vel = &Mac3d::get_u;
	
	//Initialization of variables for the interpolation
	Eigen::Vector3d indices =  index_from_coord(x,y,z);
	double x0, x1, y0, y1, z0, z1, xd, yd, zd;
	int ix0, ix1, iy0, iy1, iz0, iz1;
	double u000, u100, u110, u010, u001, u101, u111, u011, u00, u10, u01, u11, u0, u1;

	//Assigning of the indices and the values in x
	ix0 = indices[0];
	ix1 = ix0 + 1;
	x0 = (ix0-0.5)*cell_sizex_;
	x1 = (ix1-0.5)*cell_sizex_;
	
	//Assigning of the indices and the values in y and z
	//and trilinear interpolation
	if(y >= 0 && y <= sizey_ - cell_sizey_){
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);			
			u000 = (this->*get_vel)(ix0, iy0, iz0);
			u100 = (this->*get_vel)(ix1, iy0, iz0);
			u010 = (this->*get_vel)(ix0, iy1, iz0);
			u001 = (this->*get_vel)(ix0, iy0, iz1);
			u110 = (this->*get_vel)(ix1, iy1, iz0);
			u101 = (this->*get_vel)(ix1, iy0, iz1);
			u011 = (this->*get_vel)(ix0, iy1, iz1);
			u111 = (this->*get_vel)(ix1, iy1, iz1);
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
			u0 = (this->*get_vel)(ix0, 0, 0);
			u1 = (this->*get_vel)(ix1, 0, 0);
			return u0*(1-xd) + u1*xd;				
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			u0 = (this->*get_vel)(ix0, 0, L_-1);
			u1 = (this->*get_vel)(ix1, 0, L_-1);
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
			u0 = (this->*get_vel)(ix0, M_-1, 0);
			u1 = (this->*get_vel)(ix1, M_-1, 0);
			return u0*(1-xd) + u1*xd;
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			u0 = (this->*get_vel)(ix0, M_-1, L_-1);
			u1 = (this->*get_vel)(ix1, M_-1, L_-1);
			return u0*(1-xd) + u1*xd;
		}
	}
}

double Mac3d::get_interp_v(double x, double y, double z, const bool use_v_star){
	//Setting v or v*
	double (Mac3d::*get_vel)(unsigned, unsigned, unsigned);
	if (use_v_star)
		get_vel = &Mac3d::get_v_star;
	else
		get_vel = &Mac3d::get_v;
	
	//Initialization of variables for the interpolation
	Eigen::Vector3d indices =  index_from_coord(x,y,z);
	double x0, x1, y0, y1, z0, z1, xd, yd, zd;
	int ix0, ix1, iy0, iy1, iz0, iz1;
	double v000, v100, v110, v010, v001, v101, v111, v011, v00, v10, v01, v11, v0, v1;

	//Assigning of the indices and the values in y
	iy0 = indices[1];
	iy1 = iy0 + 1;
	y0 = (iy0-0.5)*cell_sizey_;
	y1 = (iy1-0.5)*cell_sizey_;

	//if (verbose) {
	//	std::cout << " + " << x << std::endl;
	//	std::cout << " + " << y << std::endl;
	//	std::cout << " + " << z << std::endl;
	//}

	//Assigning of the indices and the values in x and z
	//and trilinear interpolation
	if(x >= 0 && x <= sizex_ - cell_sizex_){
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			//if (verbose)
			//	std::cout << " + Trilinear" << std::endl;
			assign_x(x, indices[0], ix0, ix1, x0, x1);
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);			
			v000 = (this->*get_vel)(ix0, iy0, iz0);
			v100 = (this->*get_vel)(ix1, iy0, iz0);
			v010 = (this->*get_vel)(ix0, iy1, iz0);
			v001 = (this->*get_vel)(ix0, iy0, iz1);
			v110 = (this->*get_vel)(ix1, iy1, iz0);
			v101 = (this->*get_vel)(ix1, iy0, iz1);
			v011 = (this->*get_vel)(ix0, iy1, iz1);
			v111 = (this->*get_vel)(ix1, iy1, iz1);
			v00 = v000*(1-xd) + v100*xd;
			v01 = v001*(1-xd) + v101*xd;
			v10 = v010*(1-xd) + v110*xd;
			v11 = v011*(1-xd) + v111*xd;
			v0 = v00*(1-yd) + v10*yd;
			v1 = v01*(1-yd) + v11*yd;			
			return v0*(1-zd) + v1*zd;
		}
		else if (z < 0){
			assign_x(x, indices[0], ix0, ix1, x0, x1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);			
			v00 = (this->*get_vel)(ix0,iy0, 0);
			v01 = (this->*get_vel)(ix0,iy1, 0);
			v10 = (this->*get_vel)(ix1,iy0, 0);
			v11 = (this->*get_vel)(ix1,iy1, 0);
			v0 = v00*(1-xd) + v10*xd;
			v1 = v01*(1-xd) + v11*xd;			
			return v0*(1-yd) + v1*yd;
		}
		else{
			assign_x(x, indices[0], ix0, ix1, x0, x1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			v00 = (this->*get_vel)(ix0,iy0, L_-1);
			v01 = (this->*get_vel)(ix0,iy1, L_-1);
			v10 = (this->*get_vel)(ix1,iy0, L_-1);
			v11 = (this->*get_vel)(ix1,iy1, L_-1);
			v0 = v00*(1-xd) + v10*xd;
			v1 = v01*(1-xd) + v11*xd;
			return v0*(1-yd) + v1*yd;
		}
	}
	else if (x < 0){
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			//if (verbose)
			//	std::cout << " + Bilinear on left face" << std::endl;
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			v00 = (this->*get_vel)(0, iy0, iz0);
			v01 = (this->*get_vel)(0, iy0, iz1);
			v10 = (this->*get_vel)(0, iy1, iz0);
			v11 = (this->*get_vel)(0, iy1, iz1);
			//if (verbose) {
			//	std::cout << "+++ Values:" << std::endl;
			//	std::cout << " + " << v00 << std::endl;
			//	std::cout << " + " << v01 << std::endl;
			//	std::cout << " + " << v10 << std::endl;
			//	std::cout << " + " << v11 << std::endl;
			//	std::cout << "+++ Indices:" << std::endl;
			//	std::cout << " + " << iy0 << " " << iz0 << std::endl;
			//	std::cout << " + " << iy0 << " " << iz1 << std::endl;
			//	std::cout << " + " << iy1 << " " << iz0 << std::endl;
			//	std::cout << " + " << iy1 << " " << iz1 << std::endl;
			//	std::cout << " + yd: " << yd << std::endl;
			//	std::cout << " + zd: " << zd << std::endl;
			//}
			v0 = v00*(1-yd) + v10*yd;
			v1 = v01*(1-yd) + v11*yd;
			return v0*(1-zd) + v1*zd;
		}
		else if (z < 0){
			//if (verbose)
			//	std::cout << " + Linear on front-left edge" << std::endl;
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			v0 = (this->*get_vel)(0, iy0, 0);
			v1 = (this->*get_vel)(0, iy1, 0);
			return v0*(1-yd) + v1*yd;				
		}
		else{
			//if (verbose)
			//	std::cout << " + Linear on back-left edge" << std::endl;
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			v0 = (this->*get_vel)(0, iy0, L_-1);
			v1 = (this->*get_vel)(0, iy1, L_-1);
			return v0*(1-yd) + v1*yd;
		}
	}
	else{
		if(z >= 0 && z <= sizez_ - cell_sizez_){
			assign_z(z, indices[2], iz0, iz1, z0, z1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);			
			v00 = (this->*get_vel)(N_-1, iy0, iz0);
			v01 = (this->*get_vel)(N_-1, iy0, iz1);
			v10 = (this->*get_vel)(N_-1, iy1, iz0);
			v11 = (this->*get_vel)(N_-1, iy1, iz1);
			v0 = v00*(1-yd) + v10*yd;
			v1 = v01*(1-yd) + v11*yd;
			return v0*(1-zd) + v1*zd;
		}
		else if (z < 0){
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			v0 = (this->*get_vel)(N_-1, iy0, 0);
			v1 = (this->*get_vel)(N_-1, iy1, 0);
			return v0*(1-yd) + v1*yd;
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			v0 = (this->*get_vel)(N_-1, iy0, L_-1);
			v1 = (this->*get_vel)(N_-1, iy1, L_-1);
			return v0*(1-yd) + v1*yd;
		}
	}
}

double Mac3d::get_interp_w(double x, double y, double z, const bool use_w_star){
	//Setting w or w*
	double (Mac3d::*get_vel)(unsigned, unsigned, unsigned);
	if (use_w_star)
		get_vel = &Mac3d::get_w_star;
	else
		get_vel = &Mac3d::get_w;
	
	//Initialization of variables for the interpolation
	Eigen::Vector3d indices =  index_from_coord(x,y,z);
	double x0, x1, y0, y1, z0, z1, xd, yd, zd;
	int ix0, ix1, iy0, iy1, iz0, iz1;
	double w000, w100, w110, w010, w001, w101, w111, w011, w00, w10, w01, w11, w0, w1;

	//Assigning of the indices and the values in z
	iz0 = indices[2];
	iz1 = iz0 + 1;
	z0 = (iz0-0.5)*cell_sizez_;
	z1 = (iz1-0.5)*cell_sizez_;
	
	//Assigning of the indices and the values in x and y
	//and trilinear interpolation
	if(x >= 0 && x <= sizex_ - cell_sizex_){
		if(y >= 0 && y <= sizey_ - cell_sizey_){
			assign_x(x, indices[0], ix0, ix1, x0, x1);
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w000 = (this->*get_vel)(ix0, iy0, iz0);
			w100 = (this->*get_vel)(ix1, iy0, iz0);
			w010 = (this->*get_vel)(ix0, iy1, iz0);
			w001 = (this->*get_vel)(ix0, iy0, iz1);
			w110 = (this->*get_vel)(ix1, iy1, iz0);
			w101 = (this->*get_vel)(ix1, iy0, iz1);
			w011 = (this->*get_vel)(ix0, iy1, iz1);
			w111 = (this->*get_vel)(ix1, iy1, iz1);
			w00 = w000*(1-xd) + w100*xd;
			w01 = w001*(1-xd) + w101*xd;
			w10 = w010*(1-xd) + w110*xd;
			w11 = w011*(1-xd) + w111*xd;
			w0 = w00*(1-yd) + w10*yd;
			w1 = w01*(1-yd) + w11*yd;			
			return w0*(1-zd) + w1*zd;
		}
		else if (y < 0){
			assign_x(x, indices[0], ix0, ix1, x0, x1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);			
			w00 = (this->*get_vel)(ix0, 0, iz0);
			w01 = (this->*get_vel)(ix0, 0, iz1);
			w10 = (this->*get_vel)(ix1, 0, iz0);
			w11 = (this->*get_vel)(ix1, 0, iz1);
			w0 = w00*(1-xd) + w10*xd;
			w1 = w01*(1-xd) + w11*xd;			
			return w0*(1-zd) + w1*zd;
		}
		else{
			assign_x(x, indices[0], ix0, ix1, x0, x1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w00 = (this->*get_vel)(ix0, M_-1, iz0);
			w01 = (this->*get_vel)(ix0, M_-1, iz1);
			w10 = (this->*get_vel)(ix1, M_-1, iz0);
			w11 = (this->*get_vel)(ix1, M_-1, iz1);
			w0 = w00*(1-xd) + w10*xd;
			w1 = w01*(1-xd) + w11*xd;			
			return w0*(1-zd) + w1*zd;
		}
	}
	else if (x < 0){
		if(y >= 0 && y <= sizey_ - cell_sizey_){
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w00 = (this->*get_vel)(0, iy0, iz0);
			w01 = (this->*get_vel)(0, iy0, iz1);
			w10 = (this->*get_vel)(0, iy1, iz0);
			w11 = (this->*get_vel)(0, iy1, iz1);
			w0 = w00*(1-yd) + w10*yd;
			w1 = w01*(1-yd) + w11*yd;
			return w0*(1-zd) + w1*zd;
		}
		else if (y < 0){
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w0 = (this->*get_vel)(0, 0, iz0);
			w1 = (this->*get_vel)(0, 0, iz1);
			return w0*(1-zd) + w1*zd;				
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w0 = (this->*get_vel)(0, M_-1, iz0);
			w1 = (this->*get_vel)(0, M_-1, iz1);
			return w0*(1-zd) + w1*zd;
		}
	}
	else{
		if(y >= 0 && y <= sizey_ - cell_sizey_){
			assign_y(y, indices[1], iy0, iy1, y0, y1);
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);			
			w00 = (this->*get_vel)(N_-1, iy0, iz0);
			w01 = (this->*get_vel)(N_-1, iy0, iz1);
			w10 = (this->*get_vel)(N_-1, iy1, iz0);
			w11 = (this->*get_vel)(N_-1, iy1, iz1);
			w0 = w00*(1-yd) + w10*yd;
			w1 = w01*(1-yd) + w11*yd;
			return w0*(1-zd) + w1*zd;
		}
		else if (y < 0){
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w0 = (this->*get_vel)(N_-1, 0, iz0);
			w1 = (this->*get_vel)(N_-1, 0, iz1);
			return w0*(1-zd) + w1*zd;
		}
		else{
			assign_d(xd, yd, zd, x0, x1, y0, y1, z0, z1, x, y, z);
			w0 = (this->*get_vel)(N_-1, M_-1, iz0);
			w1 = (this->*get_vel)(N_-1, M_-1, iz1);
			return w0*(1-zd) + w1*zd;
		}
	}
}

//3. Pressures ---------------------------------------------------------
double Mac3d::get_pressure(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < L_)
		return *(ppressure_ + N_*j + i + N_*M_*k);
	else{ 
		std::cout << "Calling get_pressure: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

//4. Physical properties -----------------------------------------------
bool Mac3d::is_solid(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < L_)
		return *(psolid_ + N_*j + i + N_*M_*k);
	else{ 
		std::cout << "Calling is_solid: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 1;
	}
}

bool Mac3d::is_fluid(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < L_)
		return *(pfluid_ + N_*j + i + N_*M_*k);
	else{ 
		std::cout << "Calling is_fluid: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

bool Mac3d::is_empty(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < L_)
		return (!is_fluid(i,j,k) && !is_solid(i,j,k));
	else{ 
		std::cout << "Calling is_empty Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

//5. Weights -----------------------------------------------------------
double Mac3d::get_weights_u(const unsigned i, const unsigned j, const unsigned k){
	if (i < (N_+1) && j < M_ && k < L_)
		return *(pweights_u_ + (N_+1)*j + i + (N_+1)*M_*k);
	else{ 
		std::cout << "Calling get_weights_u: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

double Mac3d::get_weights_v(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < (M_+1) && k < L_)
		return *(pweights_v_ + N_*j + i + N_*(M_+1)*k);
	else{ 
		std::cout << "Calling get_weights_v: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
}

double Mac3d::get_weights_w(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < (L_+1))
		return *(pweights_w_ + N_*j + i + N_*M_*k);
	else{ 
		std::cout << "Calling get_weights_w: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
		return 0;
	}
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

 void Mac3d::index_from_coord( const double x, 
							   const double y, 
							   const double z, 
							   Mac3d::cellIdx_t &cell_idx_x, 
							   Mac3d::cellIdx_t &cell_idx_y, 
							   Mac3d::cellIdx_t &cell_idx_z )
{
	cell_idx_x = (Mac3d::cellIdx_t) (x * rcell_sizex_ + 0.5);
	cell_idx_y = (Mac3d::cellIdx_t) (y * rcell_sizey_ + 0.5);
	cell_idx_z = (Mac3d::cellIdx_t) (z * rcell_sizez_ + 0.5);
}



/************************************************************************************
*************************************SETTERS*****************************************
*************************************************************************************
//1) setters for the velocities (u, v, w, u*, v*, w*);
//2) setters for the pressures;
//3) setters for the physical properties of the cells (solid, liquid, empty);
//4) setters for the weights for the particle to grid;
*/

//1. Velocities --------------------------------------------------------
void Mac3d::set_u(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_+1) && j < M_ && k < L_)
		*(pu_ + (N_+1)*j + i + (N_+1)*M_*k) = value;
	else
		std::cout << "Calling set_u: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_v(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < (M_+1) && k < L_)
		*(pv_ + N_*j + i + N_*(M_+1)*k) = value;
	else
		std::cout << "Calling set_v: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_w(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < M_ && k < (L_+1))
		*(pw_ + N_*j + i + N_*M_*k) = value;
	else
		std::cout << "Calling set_w: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_u_star(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_+1) && j < M_ && k < L_)
		*(pu_star_ + (N_+1)*j + i + (N_+1)*M_*k) = value;
	else
		std::cout << "Calling set_u_star: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_v_star(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < (M_+1) && k < L_)
		*(pv_star_ + N_*j + i + N_*(M_+1)*k) = value;
	else
		std::cout << "Calling set_v_star: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_w_star(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < N_ && j < M_ && k < (L_+1))
		*(pw_star_ + N_*j + i + N_*M_*k) = value;
	else
		std::cout << "Calling set_w_star: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
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
	if (i < N_ && j < M_ && k < L_)
		*(ppressure_ + N_*j + i + N_*M_*k) = value;
	else
		std::cout << "Calling set_pressure: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_pressure(const Eigen::VectorXd& p) {
	std::copy(p.data(), p.data()+p.size(), ppressure_);
}

//3. Physical properties -----------------------------------------------
void Mac3d::set_solid(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < L_)
		*(psolid_ + N_*j + i + N_*M_*k) = true;
	else
		std::cout << "Calling set_solid: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_fluid(const unsigned i, const unsigned j, const unsigned k){
	if (i < N_ && j < M_ && k < L_)
		*(pfluid_ + N_*j + i + N_*M_*k) = true;
	else
		std::cout << "Calling set_fluid: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::reset_fluid() {
	std::fill(pfluid_, pfluid_ + get_num_cells(), false);
}

//4. Weights -----------------------------------------------------------
void Mac3d::set_weights_u(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_+1) && j < (M_) && k < L_)
		*(pweights_u_ + (N_+1)*j + i + (N_+1)*M_*k) = value;
	else
		std::cout << "Calling set_weights_u: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_weights_v(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_) && j < (M_+1) && k < L_)
		*(pweights_v_ + N_*j + i + N_*(M_+1)*k) = value;
	else
		std::cout << "Calling set_weights_v: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_weights_w(const unsigned i, const unsigned j, const unsigned k, double value){
	if (i < (N_) && j < M_ && k < (L_+1))
		*(pweights_w_ + N_*j + i + N_*M_*k) = value;
	else
		std::cout << "Calling set_weights_w: Index (" << i << ", " << j << ", " << k << ") out of bounds!" << std::endl;
}

void Mac3d::set_weights_to_zero(){
	std::fill(pweights_u_, pweights_u_ + (N_+1)*M_*L_, 0);
	std::fill(pweights_v_, pweights_v_ + N_*(M_+1)*L_, 0);
	std::fill(pweights_w_, pweights_w_ + N_*M_*(L_+1), 0);
}

