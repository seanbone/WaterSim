#ifndef WATERSIM_H
#define WATERSIM_H

#include "FLIP.h"
#include "MeshExporter.h"
#include "SimConfig.h"

#include <sys/stat.h> // mkdir
#include <chrono> // accurate timings
#include <list> // std::list

/**
 * This class manages the water simulation.
 * It:
 * - Initializes a MAC grid in the right configuration
 * - Initializes the list of particles
 * - Initializes a FLIP object to handle updating the state
 * - Calls the update method on FLIP at each step
 * - Uses a MeshExporter object to export meshes at each step
 */
class WaterSim {

	// Number of particles in simulation
	unsigned int m_num_particles;

	double m_dt = 0.0;         // length of timestep
	double m_time = 0.0;       // current time
	unsigned long m_step = 0;  // number of performed steps in simulation

public:

	/** Configuration object */
	SimConfig m_cfg;

	/*** Helper class members ***/

	//MAC grid data structure
	Mac3d *p_mac_grid;

	// Simulation Particles
	Particles* flip_particles;

	// FLIP simulator
	FLIP *p_flip;

	//MeshExporter
	MeshExporter *exp;

	/*** Public methods ***/

	/**
	 * Constructor
	 */
	explicit WaterSim(const SimConfig& cfg);

	/**
	 * Destructor
	 */
	~WaterSim() {
		delete p_mac_grid;
		delete p_flip;
		delete flip_particles;
	}

	/**
	 * Update simulation parameters. Requires a reset to take effect.
	 */
	void updateParams(const SimConfig& cfg);

	/**
	 * Reset class variables to reset the simulation.
	 */
	void resetMembers();

	/**
	 * Performs one simulation step of length m_dt. You have to update the time
	 * variables at the end of each step if they are necessary for your
	 * simulation.
	 */
	bool advance();

	unsigned int getNumParticles() const { return m_num_particles; }

	void setTimestep(double t) { m_dt = t; }

	double getTime() const { return m_time; }

	unsigned long getStep() const { return m_step; }

protected:

	/**
	 * Initialize the particle list according to the current settings.
	 * Note: does not delete previous particles!
	 */
	void initParticles();

	/**
	 * Initialize a new instance of the MAC grid
	 * Note: does not delete previous instance!
	 */
	void initMacGrid();

	/**
	 * Initialize a new instance of the Mesh Exporter
	 * Note: does not delete previous instance!
	 */
	void initMeshExp();

	/**
	 * Initialize a new instance of the FLIP simulator
	 * Note: does not delete previous instance!
	 */
	void initFLIP();
};

#endif // WATERSIM_H
