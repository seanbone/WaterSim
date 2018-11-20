#include "Mac2d.h"

double Mac2d::get_u(const int i, const int j){
	return *(pu_ + (N_+1)*j + i);
}

double Mac2d::get_v(const int i, const int j){
	return *(pv_ + N_*j + i);
}

double Mac2d::get_pressure(const int i, const int j){
	return *(ppressure_ + N_*j + i);
}

bool Mac2d::is_solid(const int i, const int j){
	return *(psolid_ + N_*j + i);
}

unsigned Mac2d::get_num_cells_x() {
	return N_;
}

unsigned Mac2d::get_num_cells_y() {
	return M_;
}

double Mac2d::get_cell_sizex() {
	return cell_sizex_;
}

double Mac2d::get_cell_sizey() {
	return cell_sizey_;
}

double Mac2d::get_weights_u(const int i, const int j){
	return *(pweights_u_ + (N_+1)*j + i);
}

double Mac2d::get_weights_v(const int i, const int j){
	return *(pweights_v_ + N_*j + i);
}

void Mac2d::set_u(const int i, const int j, double value){
	*(pu_ + (N_+1)*j + i) = value;
}

void Mac2d::set_v(const int i, const int j, double value){
	*(pv_ + N_*j + i) = value;
}

void Mac2d::set_pressure(const int i, const int j, double value){
	*(ppressure_ + N_*j + i) = value;
}

void Mac2d::set_weights_u(const int i, const int j, double value){
	*(pweights_u_ + (N_+1)*j + i) = value;
}

void Mac2d::set_weights_v(const int i, const int j, double value){
	*(pweights_v_ + N_*j + i) = value;
}

void Mac2d::set_velocities_to_zero(){
	for( int j = 0; j < M_; ++j ){
		for( int i = 0; i < N_; ++i ){
			set_u(i, j, 0.);
			set_v(i, j, 0.);
		}	
	}
}

void Mac2d::set_weights_to_zero(){
	for( int j = 0; j < M_+1; ++j ){
		for( int i = 0; i < N_+1; ++i ){
			set_weights_u(i, j, 0.);
			set_weights_v(i, j, 0.);
		}	
	}
}

void Mac2d::set_solid(const int i, const int j){
	*(psolid_ + N_*j + i) = true;
}

Pair_t Mac2d::index_from_coord(const double x, const double y){
	if (x > sizex_ || y > sizey_ || x < 0 || y < 0)
		std::cout << "Attention: out of the grid!" << std::endl;
	return Pair_t(int(x), int(y));
}

//Useful function to print an array
/*std::ostream& operator<< (std::ostream& out, double* list, const int list_dimension){
	for(int i = 0; i < list_dimension; ++i){
		out << *(list + i) << std::endl;
	}
	return out;
}*/

