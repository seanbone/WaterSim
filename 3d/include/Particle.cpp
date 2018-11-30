#include "Particle.h"

void Particle::set_position( const double x,
							 const double y,
							 const double z )
{
	xprev_ = x_;
	yprev_ = y_;
	zprev_ = z_;
	x_ = x;
	y_ = y;
	z_ = z;
}

void Particle::set_position( const Eigen::Vector3d& pos ){
	xprev_ = x_;
	yprev_ = y_;
	zprev_ = z_;
	x_ = pos(0);
	y_ = pos(1);
	z_ = pos(2);
}

void Particle::set_velocity( const double u,
							 const double v,
							 const double w )
{
	u_ = u;
	v_ = v;
	w_ = w;
}

void Particle::set_velocity( const Eigen::Vector3d& vel ){
	u_ = vel(0);
	v_ = vel(1);
	w_ = vel(2);
}

Eigen::Vector3d Particle::get_position(){
	Eigen::Vector3d pos (x_, y_, z_);
	return pos;
}
	
Eigen::Vector3d Particle::get_prev_position(){
	Eigen::Vector3d pos (xprev_, yprev_, zprev_);
	return pos;
}

Eigen::Vector3d Particle::get_velocity(){
	Eigen::Vector3d vel (u_, v_, w_);
	return vel;
}
