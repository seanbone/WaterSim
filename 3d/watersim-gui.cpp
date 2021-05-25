#include <igl/writeOFF.h>
#include <vector>
#include <thread>
#include "WaterSimGui.h"
#include "Gui.h"
#include "SimConfig.h"

/**
 * This class is a GUI for our water simulation. It extends the basic GUI
 * defined in Gui.h.
 */
class WaterGui : public Gui {

private:

	std::string configFile = "config.json";
	SimConfig m_cfg;

	// Simulation parameters
	bool m_export_meshes;
	int m_max_steps;
	bool m_meteor_force;
	bool m_display_meshes;
	bool m_display_mesh_edges;
	bool m_display_mesh_faces;
	int m_random_seed;

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

	// Fluid region
	double m_fluid_from_x;
	double m_fluid_from_y;
	double m_fluid_from_z;
	double m_fluid_to_x;
	double m_fluid_to_y;
	double m_fluid_to_z;

public:

	// Pointer to the simulation
	WaterSimGui *p_waterSim = nullptr;


	/** Default Constructor
	 */
	WaterGui() {

		m_cfg = SimConfig(configFile);

		readConfig();

		// Create a new simulation instance
		p_waterSim = new WaterSimGui(m_viewer, m_cfg);

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
		m_display_meshes = m_cfg.getDisplayMeshes();
		m_display_mesh_faces = m_cfg.getDisplayMeshFaces();
		m_display_mesh_edges = m_cfg.getDisplayMeshEdges();
		m_max_p_disp = m_cfg.getMaxParticlesDisplay();
		m_max_steps = m_cfg.getMaxSteps();
		m_meteor_force = m_cfg.getApplyMeteorForce();
		m_cfg.getFluidRegion(m_fluid_from_x, m_fluid_from_y, m_fluid_from_z,
					   m_fluid_to_x, m_fluid_to_y, m_fluid_to_z);
		m_random_seed = m_cfg.getRandomSeed();
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
		m_cfg.setDisplayMeshes(m_display_mesh_edges, m_display_mesh_faces);
		m_cfg.setMaxParticlesDisplay(m_max_p_disp);
		m_cfg.setMaxSteps(m_max_steps);
		m_cfg.setApplyMeteorForce(m_meteor_force);
		m_cfg.setFluidRegion(m_fluid_from_x, m_fluid_from_y, m_fluid_from_z,
		                     m_fluid_to_x, m_fluid_to_y, m_fluid_to_z);
		m_cfg.setRandomSeed(m_random_seed);
	}

	/**
	* Update the Simulation class with the parameters input to the GUI
	*/
	void updateSimulationParameters() override {
		updateConfig();

		p_simulator->setMaxSteps(m_max_steps);
		p_waterSim->setTimestep(m_dt);
		p_waterSim->updateParams(m_cfg);
	};

	/**
	 * Add rendering controls to GUI
	 */
	 void drawRenderOptionsMenu(igl::opengl::glfw::Viewer &viewer) override {
		if (ImGui::Checkbox("Mesh Wireframe", &m_display_mesh_edges) ||
			m_display_mesh_edges != viewer.data().show_lines)
		{
			p_waterSim->m_watersim.m_cfg.setDisplayMeshes(m_display_mesh_edges,
												 m_display_mesh_faces);
			viewer.data().show_lines = m_display_mesh_edges;
			for (size_t i = 0; i < viewer.data_list.size(); i++) {
				viewer.data_list[i].show_lines = viewer.data().show_lines;
			}
		}
		if (ImGui::Checkbox("Mesh Fill", &(m_display_mesh_faces)) ||
			m_display_mesh_faces != viewer.data().show_faces)
		{
			p_waterSim->m_watersim.m_cfg.setDisplayMeshes(m_display_mesh_edges,
			                                              m_display_mesh_faces);
			viewer.data().show_faces = m_display_mesh_faces;
			for (size_t i = 0; i < viewer.data_list.size(); i++) {
				viewer.data_list[i].show_faces = viewer.data().show_faces;
			}
		}
		//ImGui::Checkbox("Show stats", &m_showStats);
		if (ImGui::InputInt("Max particles display", &m_max_p_disp, 0, 0)) {
			p_waterSim->m_watersim.m_cfg.setMaxParticlesDisplay(m_max_p_disp);
			p_waterSim->updateRenderGeometry();
		}
		if (ImGui::Checkbox("Show grid", &m_display_grid)) {
			p_waterSim->m_watersim.m_cfg.setDisplayGrid(m_display_grid);
			p_waterSim->updateRenderGeometry();
		}
	 }

	/**
	* Add parameter controls to the GUI
	*/
	void drawSimulationParameterMenu() override {
		ImGui::Checkbox("Export meshes", &m_export_meshes);
		ImGui::Checkbox("Randomize particles", &m_jitter_particles);
		ImGui::InputInt("Random seed", &m_random_seed);
		ImGui::Checkbox("Meteor force", &m_meteor_force);
		ImGui::InputDouble("Alpha", &m_alpha, 0, 0);
		ImGui::InputDouble("Timestep [s]", &m_dt, 0, 0);
		ImGui::InputInt("Max Steps", &m_max_steps);
		ImGui::InputDouble("Density [kg/m^3]", &m_density, 0, 0);
		ImGui::InputDouble("Gravity [m/s^2]", &m_gravity, 0, 0);
		ImGui::InputInt("Grid resolution X", &m_grid_res_x, 0, 0);
		ImGui::InputInt("Grid resolution Y", &m_grid_res_y, 0, 0);
		ImGui::InputInt("Grid resolution Z", &m_grid_res_z, 0, 0);
		ImGui::InputDouble("X Size [m]", &m_system_size_x, 0, 0);
		ImGui::InputDouble("Y Size [m]", &m_system_size_y, 0, 0);
		ImGui::InputDouble("Z Size [m]", &m_system_size_z, 0, 0);

		if (ImGui::CollapsingHeader("Fluid region",
		                            ImGuiTreeNodeFlags_None)) {
			ImGui::InputDouble("From X [m]", &m_fluid_from_x);
			ImGui::InputDouble("From Y [m]", &m_fluid_from_y);
			ImGui::InputDouble("From Z [m]", &m_fluid_from_z);
			ImGui::InputDouble("To X [m]", &m_fluid_to_x);
			ImGui::InputDouble("To Y [m]", &m_fluid_to_y);
			ImGui::InputDouble("To Z [m]", &m_fluid_to_z);
		}

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
			resetSimulation();
		}
	}


	/**
	* Add maximum and minimum pressure to the GUI
	*/
	void drawSimulationStats() override {

		// Print the maximal and the minimal pressure at every time-step
		double pressure_max = 0;
		double pressure_min = 0;

		// Get total number of cells on each axis
		//unsigned nx = (*p_waterSim).p_mac_grid->get_num_cells_x();
		//unsigned ny = (*p_waterSim).p_mac_grid->get_num_cells_y();
		//unsigned nz = (*p_waterSim).p_mac_grid->get_num_cells_z();
		int nx, ny, nz;
		p_waterSim->m_watersim.m_cfg.getGridResolution(nx, ny, nz);

		// Iterate over all grid-cells
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {

					// Find the maximal/minimal pressure
					double temp = p_waterSim->m_watersim.p_mac_grid->get_pressure(i, j, k);

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
		ImGui::Text("%d particles", p_waterSim->m_watersim.getNumParticles());
		ImGui::Text("%d cells (%d x %d x %d)", nx*ny*nz, nx, ny, nz);
	}


	/**
	* Override Gui::resetSimulation to also update simulation parameters.
	*/
	void resetSimulation() override {
		updateConfig();
		readConfig(); // Read config again to make sure invariants are respected
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

