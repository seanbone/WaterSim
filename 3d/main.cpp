#include <igl/writeOFF.h>
#include <Eigen/Eigen>
#include <vector>
#include <thread>
#include "WaterSim.h"
#include "Gui.h"
#include "Simulator.h"


/*
 * This helper function is a utility to determine which
 * cells should be initialized as containing fluid at the
 * start of the simulation.
 * Note: use this to set the initial state of fluid!
 */
std::vector<bool> select_fluid_cells(size_t nx, size_t ny, size_t nz) {
	std::vector<bool> is_fluid(nx*ny*nz, false);

	//~ std::fill(is_fluid.begin(), is_fluid.end(), true);
	
	for (unsigned k = 3; k < 13; k++) {
		for (unsigned j = 5; j < ny; j++) {
			for (unsigned i = 3; i < 13; i++) {
				is_fluid[i + j*nx + nx*ny*k] = true;
			}
		}
	}

   return is_fluid;
}


/*
 * This class is a GUI for our water simulation. It extends the basic GUI
 * defined in Gui.h.
 */
class WaterGui : public Gui {
private:
	// Simulation parameters
	bool m_export_meshes = false;
	bool m_show_pressures = false;
	bool m_display_velocity_arrows = false;
	int m_export_fps = 30;

	double m_system_size_x = 10; // X dimension of system in m
	double m_system_size_y = 10; // Y dimension of system in m
	double m_system_size_z = 10; // Z dimension of system in m

	int m_grid_res_x = 15; // Number of cells on X axis
	int m_grid_res_y = 15; // Number of cells on Y axis
	int m_grid_res_z = 15; // Number of cells on Z axis
	// Whether to randomize particle positions
	bool m_jitter_particles = false; 

	double m_dt = 0.025;

	double m_alpha = 0.01;

	double m_density = 1000.0;  // Fluid density in kg/m^3
	double m_gravity = 9.81; // Acceleration of gravity in m/s^2

	// PNG export options
	bool m_export_png = false;
	int m_png_sx = 1280;
	int m_png_sy = 800;
	int m_max_pngs = 500;

	// Other members
	bool m_display_grid = true;
	int m_cells_x;
	int m_cells_y;
	int m_cells_z;

public:
	WaterSim *p_waterSim = NULL;  // pointer to the simulation

	WaterGui() {
		// Initialize fluid flags
		std::vector<bool> is_fluid = select_fluid_cells(m_grid_res_x, m_grid_res_y, m_grid_res_z);


		// create a new simulation instance
		p_waterSim = new WaterSim(m_viewer, m_display_grid, m_grid_res_x, m_grid_res_y, m_grid_res_z,
								  m_system_size_x, m_system_size_y, m_system_size_z,
								  m_density, m_gravity, m_alpha,
								  m_show_pressures, m_display_velocity_arrows,
								  std::move(is_fluid), m_jitter_particles,
								  m_export_png, m_png_sx, m_png_sy, m_max_pngs,
								  m_export_meshes);

		// set this simulation as the simulation that is running in our GUI
		setSimulation(p_waterSim);

		// start the GUI
		start();
	}
      

	/**
	* Update the Simulation class with the parameters input to the GUI
	*/
	virtual void updateSimulationParameters() override {
		// Dimensions of system have potentially changed
		//  -> reevaluate fluid population
		auto is_fluid = select_fluid_cells(m_grid_res_x, m_grid_res_y, m_grid_res_z);

		p_waterSim->setTimestep(m_dt);
		p_waterSim->updateParams(m_display_grid, m_grid_res_x, m_grid_res_y, m_grid_res_z,
								 m_system_size_x, m_system_size_y, m_system_size_z,
								 m_density, m_gravity, m_alpha, m_show_pressures, 
								 m_display_velocity_arrows, std::move(is_fluid),
								 m_jitter_particles,
								 m_export_png, m_png_sx, m_png_sy, m_max_pngs,
								 m_export_meshes);
	};

	/**
	* Add parameter controls to the GUI
	*/
	virtual void drawSimulationParameterMenu() override {
		ImGui::Checkbox("Display grid", &m_display_grid);
		ImGui::Checkbox("Export meshes", &m_export_meshes);
		ImGui::Checkbox("Show pressure field", &m_show_pressures);
		ImGui::Checkbox("Display velocity arrows", &m_display_velocity_arrows);      
		ImGui::Checkbox("Randomize particles", &m_jitter_particles);
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
		ImGui::Checkbox("Export PNGs (Warning: lags GUI!)", &m_export_png);
		ImGui::InputInt("PNG size x", &m_png_sx, 0, 0);
		ImGui::InputInt("PNG size y", &m_png_sy, 0, 0);
		ImGui::InputInt("Max PNGs", &m_max_pngs, 0, 0);
	}
      
	/**
	* Add maximum and minimum pressure to the GUI
	*/
	virtual void drawSimulationStats() override{	        
		//Print the maximal and the minimal pressure at every time-step
		double pressure_max;
		double pressure_min;
		unsigned nx = (*p_waterSim).p_mac_grid->get_num_cells_x();
		unsigned ny = (*p_waterSim).p_mac_grid->get_num_cells_y();
		unsigned nz = (*p_waterSim).p_mac_grid->get_num_cells_z();
		for (unsigned k = 0; k < nz; k++) {	
			for (unsigned j = 0; j < nx; j++) {
				for (unsigned i = 0; i < ny; i++) {
					double temp = (*p_waterSim).p_mac_grid->get_pressure(i, j, k);
					if(i == 0 && j == 0 && k == 0){
						pressure_max = temp;
						pressure_min = temp;
					}
					else{
						if(temp > pressure_max)
							pressure_max = temp;
						if(temp < pressure_min)
							pressure_min = temp;
					}
				}
			}
		}
		ImGui::Text("Maximal pressure: %.5f", pressure_max);
		ImGui::Text("Minimal pressure: %.5f", pressure_min);
	}

	/**
	* Override Gui::resetSimulation to also update simulation parameters.
	*/
	void resetSimulation() override {
		updateSimulationParameters();
		Gui::resetSimulation();
	}

};

int main() { //int argc, char *argv[]) {
   // Tell Eigen we're using multithreading
   // https://eigen.tuxfamily.org/dox/TopicMultiThreading.html
   Eigen::initParallel();

   // Create a new instance of the GUI for the simulation
   new WaterGui();

   return 0;
}

