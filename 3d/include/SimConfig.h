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
		 * Set default values for missing entries.
		 * @param hard If set to true, will override entries already set with defaults.
		 */
		void setDefaults(bool hard = false);

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

		/** Maximum number of particles to display
		 * Needed to prevent GUI lag for large sims
		 * If smaller than total number of particles,
		 * they will be selected with a stride of num_particles/max_p_disp
		 * to visualize the state of the sim.
		 */
		void setMaxParticlesDisplay(int maxParticles);
		int getMaxParticlesDisplay() const;
};

#endif //WATERSIM_SIMCONFIG_H
