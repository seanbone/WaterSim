#include <igl/writeOFF.h>
#include <Eigen/Eigen>
#include <vector>
#include <thread>
#include "WaterSim.h"
#include "Gui.h"
#include "Simulator.h"
#include "SimConfig.h"


/**
 * This helper function is a utility to determine which
 * cells should be initialized as containing fluid at the
 * start of the simulation.
 * Note: use this to set the initial state of fluid!
 * Params:
 * - nx, ny and nz are the number of cells on each axis
 */
std::vector<bool> select_fluid_cells(size_t nx, size_t ny, size_t nz) {

	// Vector of flags that determine which cells are fluid or not
	std::vector<bool> is_fluid(nx * ny * nz, false);

	// Iterate over all the cells that have to contain fluid at the 
	// beginning of the simulation
	for (unsigned k = 0; k < nz; k++) {
		for (unsigned j = 0; j < 10; j++) {
			for (unsigned i = 0; i < nx; i++) {
				is_fluid[i + j * nx + nx * ny * k] = true;
			}
		}
	}

	return is_fluid;
}


/**
 * This class is a GUI for our water simulation. It extends the basic GUI
 * defined in Gui.h.
 */
class WaterGui : public Gui {

private:

	std::string configFile = "../config.json";
	SimConfig m_cfg;

	// Simulation parameters
	bool m_export_meshes;

	// X, Y and Z dimension of system in meters
	double m_system_size_x;
	double m_system_size_y;
	double m_system_size_z;

	// Number of grid-cells on X, Y and Z axis
	int m_grid_res_x;
	int m_grid_res_y;
	int m_grid_res_z;

	// Whether to randomize particle positions
	bool m_jitter_particles;

	// Timestep
	double m_dt;

	// FLIP: alpha = 0.
	// PIC: alpha = 1.
	double m_alpha;

	// Fluid density in kg/m^3
	double m_density;

	// Acceleration of gravity in m/s^2
	double m_gravity;

	// Other members
	bool m_display_grid;

	// Maximum number of particles to display
	// Needed to prevent GUI lag for large sims
	// If smaller than total number of particles,
	// they will be selected with a stride of num_particles/max_p_disp
	// to visualize the state of the sim.
	int m_max_p_disp = 4242;

public:

	// Pointer to the simulation
	WaterSim *p_waterSim = nullptr;


	/** Default Constructor
	 */
	WaterGui() {

		m_cfg = SimConfig(configFile);

		readConfig();

		// Initialize fluid flags
		std::vector<bool> is_fluid = select_fluid_cells(m_grid_res_x, m_grid_res_y, m_grid_res_z);


		// Create a new simulation instance
		p_waterSim = new WaterSim(m_viewer, m_cfg, is_fluid);

		// Set this simulation as the simulation that is running in our GUI
		setSimulation(p_waterSim);

		// Start the GUI
		start();
	}

	/**
	 * Update GUI with values read from SimConfig object.
	 */
	void readConfig() {
		m_export_meshes = m_cfg.getExportMeshes();
		m_cfg.getSystemSize(m_system_size_x, m_system_size_y, m_system_size_z);
		m_cfg.getGridResolution(m_grid_res_x, m_grid_res_y, m_grid_res_z);
		m_jitter_particles = m_cfg.getJitterParticles();
		m_dt = m_cfg.getTimeStep();
		m_alpha = m_cfg.getAlpha();
		m_density = m_cfg.getDensity();
		m_gravity = m_cfg.getGravity();
		m_display_grid = m_cfg.getDisplayGrid();
		m_max_p_disp = m_cfg.getMaxParticlesDisplay();
	}

	/**
	 * Update internal SimConfig object with values input to GUI.
	 */
	void updateConfig() {
		m_cfg.setExportMeshes(m_export_meshes);
		m_cfg.setSystemSize(m_system_size_x, m_system_size_y, m_system_size_z);
		m_cfg.setGridResolution(m_grid_res_x, m_grid_res_y, m_grid_res_z);
		m_cfg.setJitterParticles(m_jitter_particles);
		m_cfg.setTimeStep(m_dt);
		m_cfg.setAlpha(m_alpha);
		m_cfg.setDensity(m_density);
		m_cfg.setGravity(m_gravity);
		m_cfg.setDisplayGrid(m_display_grid);
		m_cfg.setMaxParticlesDisplay(m_max_p_disp);
	}

	/**
	* Update the Simulation class with the parameters input to the GUI
	*/
	void updateSimulationParameters() override {
		// Dimensions of system have potentially changed
		//  -> reevaluate fluid population
		auto is_fluid = select_fluid_cells(m_grid_res_x, m_grid_res_y, m_grid_res_z);

		updateConfig();

		p_waterSim->setTimestep(m_dt);
		p_waterSim->updateParams(m_cfg, is_fluid);
	};


	/**
	* Add parameter controls to the GUI
	*/
	virtual void drawSimulationParameterMenu() override {
		ImGui::Checkbox("Display grid", &m_display_grid);
		ImGui::Checkbox("Export meshes", &m_export_meshes);
		ImGui::Checkbox("Randomize particles", &m_jitter_particles);
		ImGui::InputInt("Max particles display", &m_max_p_disp, 0, 0);
		ImGui::InputDouble("Alpha", &m_alpha, 0, 0);
		ImGui::InputDouble("Timestep [s]", &m_dt, 0, 0);
		ImGui::InputDouble("Density [kg/m^3]", &m_density, 0, 0);
		ImGui::InputDouble("Gravity [m/s^2]", &m_gravity, 0, 0);
		ImGui::InputInt("Grid resolution X", &m_grid_res_x, 0, 0);
		ImGui::InputInt("Grid resolution Y", &m_grid_res_y, 0, 0);
		ImGui::InputInt("Grid resolution Z", &m_grid_res_z, 0, 0);
		ImGui::InputDouble("X Size [m]", &m_system_size_x, 0, 0);
		ImGui::InputDouble("Y Size [m]", &m_system_size_y, 0, 0);
		ImGui::InputDouble("Z Size [m]", &m_system_size_z, 0, 0);

		if (ImGui::Button("Save configuration", ImVec2(-1, 0))) {
			updateConfig();
			m_cfg.writeFile(configFile);
		}
		if (ImGui::Button("Reload config file", ImVec2(-1, 0))) {
			m_cfg.readFile(configFile);
			readConfig();
		}
		if (ImGui::Button("Reset to defaults", ImVec2(-1, 0))) {
			m_cfg.setDefaults(true);
			readConfig();
		}
	}


	/**
	* Add maximum and minimum pressure to the GUI
	*/
	virtual void drawSimulationStats() override {

		// Print the maximal and the minimal pressure at every time-step
		double pressure_max = 0;
		double pressure_min = 0;

		// Get total number of cells on each axis
		unsigned nx = (*p_waterSim).p_mac_grid->get_num_cells_x();
		unsigned ny = (*p_waterSim).p_mac_grid->get_num_cells_y();
		unsigned nz = (*p_waterSim).p_mac_grid->get_num_cells_z();

		// Iterate over all grid-cells
		for (unsigned k = 0; k < nz; k++) {
			for (unsigned j = 0; j < ny; j++) {
				for (unsigned i = 0; i < nx; i++) {

					// Find the maximal/minimal pressure
					double temp = (*p_waterSim).p_mac_grid->get_pressure(i, j, k);

					if (i == 0 && j == 0 && k == 0) {
						pressure_max = temp;
						pressure_min = temp;
					} else {

						if (temp > pressure_max)
							pressure_max = temp;

						if (temp < pressure_min)
							pressure_min = temp;
					}
				}
			}
		}

		ImGui::Text("Max pressure: %.5f", pressure_max);
		ImGui::Text("Min pressure: %.5f", pressure_min);
	}


	/**
	* Override Gui::resetSimulation to also update simulation parameters.
	*/
	void resetSimulation() override {
		updateSimulationParameters();
		Gui::resetSimulation();
	}
};

int main() {

	// Tell Eigen we are using multithreading
	// https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
	Eigen::initParallel();

	// Create a new instance of the GUI for the simulation
	new WaterGui();

	return 0;
}

