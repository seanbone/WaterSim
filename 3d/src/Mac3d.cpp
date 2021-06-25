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
}

void Mac3d::initArrays() {
    const unsigned chunk_size = 32;
	ppressure_ = new (std::align_val_t(chunk_size)) double[N_*M_*L_];
	std::fill(ppressure_, ppressure_+N_*M_*L_, 0.);

	pu_ = new (std::align_val_t(chunk_size)) double[(N_+1)*M_*L_];
	std::fill(pu_, pu_+(N_+1)*M_*L_, 0.);

	pu_star_ = new (std::align_val_t(chunk_size)) double[(N_+1)*M_*L_];
	std::fill(pu_star_, pu_star_+(N_+1)*M_*L_, 0.);

	pv_ = new (std::align_val_t(chunk_size))double[N_*(M_+1)*L_];
	std::fill(pv_, pv_+N_*(M_+1)*L_, 0.);

	pv_star_ = new (std::align_val_t(chunk_size))double[N_*(M_+1)*L_];
	std::fill(pv_star_, pv_star_+N_*(M_+1)*L_, 0.);
	
	pw_ = new (std::align_val_t(chunk_size))double[N_*M_*(L_+1)];
	std::fill(pw_, pw_+N_*M_*(L_+1), 0.);

	pw_star_ = new (std::align_val_t(chunk_size))double[N_*M_*(L_+1)];
	std::fill(pw_star_, pw_star_+N_*M_*(L_+1), 0.);

	psolid_ = new (std::align_val_t(chunk_size))bool[N_*M_*L_];
	std::fill(psolid_, psolid_+N_*M_*L_, 0.);

	pfluid_ = new (std::align_val_t(chunk_size))bool[N_*M_*L_];
	std::fill(pfluid_, pfluid_+N_*M_*L_, 0.);

	pweights_u_ = new (std::align_val_t(chunk_size))double[(N_+1)*M_*L_];
	std::fill(pweights_u_, pweights_u_+(N_+1)*M_*L_, 0.);

	pweights_v_ = new (std::align_val_t(chunk_size))double[N_*(M_+1)*L_];
	std::fill(pweights_v_, pweights_v_+N_*(M_+1)*L_, 0.);
	
	pweights_w_ = new(std::align_val_t(chunk_size)) double[N_*M_*(L_+1)];
	std::fill(pweights_w_, pweights_w_+N_*M_*(L_+1), 0.);

}

void Mac3d::initAdiag() {

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
				A_diag_val.push_back(count);
			}
		}
	}
	// Notice that the elements are sorted by index
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

bool Mac3d::is_fluid(const unsigned i, const unsigned j, const unsigned k) const{
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


//6. Diagonal of A -----------------------------------------------------
const std::vector< Mac3d::Triplet_t >& Mac3d::get_a_diag() const {
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

void Mac3d::set_weights_to_zero(){
	std::fill(pweights_u_, pweights_u_ + (N_+1)*M_*L_, 0);
	std::fill(pweights_v_, pweights_v_ + N_*(M_+1)*L_, 0);
	std::fill(pweights_w_, pweights_w_ + N_*M_*(L_+1), 0);
}

