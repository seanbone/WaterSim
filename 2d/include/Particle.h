#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Dense>

struct Particle {
	private:
	
	// Store the position for every coordinate 
	double x_;
	double y_;
	double z_;
	
	// Store the velocity for every coordinate
	double u_;
	double v_;
	double w_;

	public:	
	
	// Default Constructor
	Particle() {}
	
	// Copy Constructor
	//~ TODO
	
	// Initialize positions and velocities (zero by default)
	Particle( const double x, 
			  const double y,
			  const double z,
			  const double u = 0.,
			  const double v = 0.,
			  const double w = 0. )
		: x_(x), y_(y), z_(z), u_(u), v_(v), w_(w) {} 
	
	// Destructor
	~Particle() {}
	
	// Setters
	void set_Positions( const double x,
						const double y,
						const double z );
						
	void set_Positions( const Eigen::Vector3d& pos );
	
	void set_Velocities( const double u,
						 const double v,
						 const double w );
						
	void set_Velocities( const Eigen::Vector3d& vel );
						 
	// Getters
	Eigen::Vector3d get_Positions();
	
	Eigen::Vector3d get_Velocities();
	
};

#endif
