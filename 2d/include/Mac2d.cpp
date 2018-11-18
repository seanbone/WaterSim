#include "Mac2d.h"


Mac2d::Mac2d()
	: N_(0), M_(0), sizex_(0.), sizey_(0.){}
	

Mac2d::Mac2d(const int n, const int m, const double dx, const double dy)
	: N_(n), M_(m), sizex_(dx), sizey_(dy){
	
	cell_sizex_ = sizex_/(1.*N_);
	cell_sizey_ = sizey_/(1.*M_);
	ppressure_ = new double[N_*M_];
	pu_ = new double[(N_+1)*M_];
	pv_ = new double[N_*(M_+1)];
	psolid_ = new bool[N_*M_];
	//manca da inizializzare la lista A_diag!		
}

~Mac2d(){
	delete ppressure_;
	delete pu_;
	delete pv_;
	delete psolid_;
}

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

Mac2d::Pair_t Mac2d::index_from_coord(const double x, const double y){
	return Pair(0,0); //da modificare!!
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

void Mac2d::set_solid(const int i, const int j){
	*(psolid_ + N_*j + i) = true;
}

