#ifndef WATERSIM_SIMCONFIG_H
#define WATERSIM_SIMCONFIG_H

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "string_format.h"

/**
 * A class to read and write configuration files for the FLIP simulation.
 */
class SimConfig {
		using json = nlohmann::json;
	private:
		// JSON object to store data
		json m_config;

		/** Sort a and b such that a <= b */
		static inline void simpleSort(double& a, double& b) {
			if (a > b) {
				double tmp = a;
				a = b;
				b = tmp;
			}
		}
	public:
		SimConfig() {
			setDefaults(true);
		}

		/** Read JSON file as initial values. */
		explicit SimConfig(const std::string& json_file) {
			readFile(json_file);
		}

		/**
		 * Read a file and parse it as a JSON config file.
		 * @param json_file
		 * @return True on success, false on file read failure.
		 */
		bool readFile(const std::string& json_file);

		/**
		 * Write configuration file.
		 * @param json_file
		 * @return Return false on file write failure, true on success.
		 */
		bool writeFile(const std::string& json_file);

		/**
		 * Give access to JSON serialization.
		 */
		std::string toString() const { return m_config.dump(4); }

		/**
		 * Set default values for missing entries.
		 * @param hard If set to true, will override entries already set with defaults.
		 */
		void setDefaults(bool hard = false);

		/*****  Sim Configuration Parameters  *****/

		/** Whether to export meshes */
		void setExportMeshes(bool v);
		bool getExportMeshes() const;

		/** System size in meters */
		void setSystemSize(double x, double y, double z);
		void getSystemSize(double& x, double& y, double& z) const;

		/** Number of grid-cells on X, Y and Z axes */
		void setGridResolution(int x, int y, int z);
		void getGridResolution(int& x, int& y, int& z) const;

		/** Whether to randomize particle positions on startup */
		void setJitterParticles(bool jitter);
		bool getJitterParticles() const;

		/** Simulation time step */
		void setTimeStep(double dt);
		double getTimeStep() const;

		/** Blending factor between PIC and FLIP.
		 * Pure FLIP: alpha = 0.0
		 * Pure PIC: alpha = 1.0
		 */
		void setAlpha(double alpha);
		double getAlpha() const;

		/** Density of the simulated fluid in kg/m^3 */
		void setDensity(double density);
		double getDensity() const;

		/** Gravity acting on the fluid in m/s^2 */
		void setGravity(double gravity);
		double getGravity() const;

		/** Whether to display the grid. Only affects GUI mode. */
		void setDisplayGrid(bool display);
		bool getDisplayGrid() const;

		/** Whether to display the meshes in the viewport. */
		void setDisplayMeshes(bool displayEdges, bool displayFaces);
		bool getDisplayMeshes() const;
		bool getDisplayMeshEdges() const;
		bool getDisplayMeshFaces() const;

		/** Maximum number of particles to display
		 * Needed to prevent GUI lag for large sims
		 * If smaller than total number of particles,
		 * they will be selected with a stride of num_particles/max_p_disp
		 * to visualize the state of the sim.
		 */
		void setMaxParticlesDisplay(int maxParticles);
		int getMaxParticlesDisplay() const;

		/**
		 * Maximum number of steps to perform in simulation.
		 */
		void setMaxSteps(int maxSteps);
		int getMaxSteps() const;

		/**
		 * Whether or not to apply the force replicating the meteor splash.
		 */
		 void setApplyMeteorForce(bool explode);
		 bool getApplyMeteorForce() const;

		 /**
		  * Define the domain initially filled with fluid. Coords in meters.
		  * Two points (from_x, from_y, from_z) and (to_x, to_y, to_z) define an
		  * axis-aligned bounding box. All cells within that area are flagged as
		  * containing fluid at the start of the simulation.
		  * Note: coordinates are rearranged such that
		  * from_x <= to_x
		  * from_y <= to_y
		  * from_z <= to_z
		  */
		  void setFluidRegion(double from_x, double from_y, double from_z,
						double to_x, double to_y, double to_z);
		  void getFluidRegion(double& from_x, double& from_y, double& from_z,
						double& to_x, double& to_y, double& to_z) const;

		  /**
		   * Set a random seed for the simulation's random number generators.
		   * If < 0, it should be ignored and a non-deterministic seed used instead.
		   */
		  void setRandomSeed(int seed);
		  int getRandomSeed();
};

#endif //WATERSIM_SIMCONFIG_H
