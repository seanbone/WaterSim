#include "Mac2d.h"

//************************************************************************************
//***********************************GETTERS******************************************
//************************************************************************************
//1) getters for the layout properties of the grid (cell sizes and sizes of the grid);
//2) getters for the velocities (u, v, u*, v*);
//3) getters for the pressures;
//4) getters for the physical properties of the cells (solid, liquid, empty);
//5) getters for the weights for the particle to grid;
//6) getters for the diagonal of the pressure matrix A;
//7) getters for the indices of the cell to which a particle belong.

//1. Layout properties of the grid -------------------------------------
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
double Mac2d::get_u(const int i, const int j){
	assert(i < (N_+1) && j < (M_) && "Index out of bounds!");
	return *(pu_ + (N_+1)*j + i);
}

double Mac2d::get_v(const int i, const int j){
	assert(i < (N_) && j < (M_+1) && "Index out of bounds!");
	return *(pv_ + N_*j + i);
}

double Mac2d::get_u_star(const int i, const int j) {
	assert(i < (N_+1) && j < (M_) && "Index out of bounds!");
	return *(pu_star_ + (N_+1)*j + i);
}

double Mac2d::get_v_star(const int i, const int j) {
	assert(i < (N_) && j < (M_+1) && "Index out of bounds!");
	return *(pv_star_ + N_*j + i);
}

//3. Pressures ---------------------------------------------------------
double Mac2d::get_pressure(const int i, const int j){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	return *(ppressure_ + N_*j + i);
}

//4. Physical properties -----------------------------------------------
bool Mac2d::is_solid(const int i, const int j){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	return *(psolid_ + N_*j + i);
}

bool Mac2d::is_fluid(const int i, const int j){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	return *(pfluid_ + N_*j + i);
}

bool Mac2d::is_empty(const int i, const int j){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	return (!is_fluid(i,j) && !is_solid(i,j));
}

//5. Weights -----------------------------------------------------------
double Mac2d::get_weights_u(const int i, const int j){
	assert(i < (N_+1) && j < (M_) && "Index out of bounds!");
	return *(pweights_u_ + (N_+1)*j + i);
}

double Mac2d::get_weights_v(const int i, const int j){
	assert(i < (N_) && j < (M_+1) && "Index out of bounds!");
	return *(pweights_v_ + N_*j + i);
}

//6. Diagonal of A -----------------------------------------------------
const std::vector< Mac2d::Triplet_t >& Mac2d::get_a_diag() {
	return A_diag_;
}

//7. Indices from coordinate -------------------------------------------
Mac2d::Pair_t Mac2d::index_from_coord(const double x, const double y){
	assert((x < sizex_ || y < sizey_ || x > 0 || y > 0)
			&& "Attention: out of the grid!");
	return Pair_t(int(x), int(y));
}



//************************************************************************************
//***********************************SETTERS******************************************
//************************************************************************************
//1) setters for the velocities (u, v, u*, v*);
//2) setters for the pressures;
//3) setters for the physical properties of the cells (solid, liquid, empty);
//4) setters for the weights for the particle to grid;

//1. Velocities --------------------------------------------------------
void Mac2d::set_u(const int i, const int j, double value){
	assert(i < (N_+1) && j < (M_) && "Index out of bounds!");
	*(pu_ + (N_+1)*j + i) = value;
}

void Mac2d::set_v(const int i, const int j, double value){
	assert(i < (N_) && j < (M_+1) && "Index out of bounds!");
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
void Mac2d::set_pressure(const int i, const int j, double value){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	*(ppressure_ + N_*j + i) = value;
}

void Mac2d::set_pressure(const Eigen::VectorXd& p) {
	std::copy(p.data(), p.data()+p.size(), ppressure_);
}

//3. Physical properties -----------------------------------------------
void Mac2d::set_solid(const int i, const int j){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	*(psolid_ + N_*j + i) = true;
}

void Mac2d::set_fluid(const int i, const int j){
	assert(i < N_ && j < M_ && "Index out of bounds!");
	*(pfluid_ + N_*j + i) = true;
}

void Mac2d::reset_fluid() {
	std::fill(pfluid_, pfluid_ + get_num_cells(), false);
}

//4. Weights -----------------------------------------------------------
void Mac2d::set_weights_u(const int i, const int j, double value){
	assert(i < (N_+1) && j < (M_) && "Index out of bounds!");
	*(pweights_u_ + (N_+1)*j + i) = value;
}

void Mac2d::set_weights_v(const int i, const int j, double value){
	assert(i < (N_) && j < (M_+1) && "Index out of bounds!");
	*(pweights_v_ + N_*j + i) = value;
}

void Mac2d::set_weights_to_zero(){
	std::fill(pweights_u_, pweights_u_ + (N_+1)*M_, 0);
	std::fill(pweights_v_, pweights_v_ + (M_+1)*N_, 0);
}
