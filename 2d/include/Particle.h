#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Dense>

struct Particle {
	private:
	
	// Store the position for every coordinate 
	double x_;
	double y_;
	double z_;
	
	// Store the velocity for every direction
	double u_;
	double v_;
	double w_;

	// Store the previous position for every direction
	//  -> always equal to the position before the last time
	//     set_position was called
	double xprev_;
	double yprev_;
	double zprev_;

	public:
	
	// Default Constructor
	Particle() {}
	
	// Copy Constructor
	//~ TODO: Particle copy constructor
	
	// Initialize positions and velocities (zero by default)
	Particle( const double x, 
			  const double y,
			  const double z,
			  const double u = 0.,
			  const double v = 0.,
			  const double w = 0. )
		: x_(x), y_(y), z_(z),
		  u_(u), v_(v), w_(w),
		  xprev_(x), yprev_(y), zprev_(z) {} 
	
	// Destructor
	~Particle() {}
	
	// Setters
	void set_position( const double x,
						const double y,
						const double z );

	void set_position( const Eigen::Vector3d& pos );

	void set_velocity( const double u,
						 const double v,
						 const double w );

	void set_velocity( const Eigen::Vector3d& vel );

	// Getters
	Eigen::Vector3d get_position();

	Eigen::Vector3d get_velocity();

	Eigen::Vector3d get_prev_position();
};

#endif
