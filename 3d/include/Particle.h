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
	
	/** Default Constructor
	 */
	Particle() {}
	
	/** Initialize positions and velocities
	 * Params:
	 * - x, y and z are the coordinates of the particle
	 * - u, v and w are the velocity components, that are set to zero by
	 * 	 default
	 */
	Particle( const double x, 
			  const double y,
			  const double z,
			  const double u = 0.,
			  const double v = 0.,
			  const double w = 0. )
		: x_(x), y_(y), z_(z),
		  u_(u), v_(v), w_(w),
		  xprev_(x), yprev_(y), zprev_(z) {} 
	
	/** Destructor
	 */
	~Particle() {}
	
	/** Set the position of a particle to x, y and z
	 * Params:
	 * - x, y and z are the coordinates to which the particle position 
	 * 	 should be setted
	 */
	void set_position( const double x,
					   const double y,
					   const double z );
	
	/** Set the position of a particle using an Eigen Vector
	 * Params:
	 * - pos is an Eigen Vector containing the coordinates to which the 
	 * 	 particle position should be setted
	 */
	void set_position( const Eigen::Vector3d& pos );

	/** Set the velocity components of a particle to u, v and w
	 * Params:
	 * - u, v and w are the velocity components to which the particle 
	 * 	 velocity should be setted
	 */
	void set_velocity( const double u,
					   const double v,
					   const double w );
	
	/** Set the velocity components of a particle using an Eigen Vector
	 * Params:
	 * - vel is an Eigen Vector containing the velocity components to 
	 * 	 which the particle velocity should be setted
	 */
	void set_velocity( const Eigen::Vector3d& vel );

	/** Get the position of a particle as an Eigen Vector
	 */
	Eigen::Vector3d get_position();

	/** Get the velocity of a particle as an Eigen Vector
	 */
	Eigen::Vector3d get_velocity();

	/** Get the position of a particle, at timestep t-1, as an Eigen 
	 * 	Vector
	 */
	Eigen::Vector3d get_prev_position();
};

#endif
