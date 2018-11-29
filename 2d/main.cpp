#include <igl/writeOFF.h>
#include <thread>
#include "WaterSim.h"
#include "Gui.h"
#include "Simulator.h"

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

      double m_system_size_x = 20; // X dimension of system in m
      double m_system_size_y = 20; // Y dimension of system in m
      
      int m_grid_res_x = 20; // Number of cells on X axis
      int m_grid_res_y = 20; // Number of cells on Y axis

      double m_dt = 0.005 * std::sqrt((m_grid_res_x + m_grid_res_y) * 0.5);

      double m_alpha = 0.05;

      double m_density = 1000.0;  // Fluid density in kg/m^3
      double m_gravity = 9.81; // Acceleration of gravity in m/s^2

      // Other members
      bool m_display_grid = true;
      int m_cells_x;
      int m_cells_y;

   public:
      WaterSim *p_waterSim = NULL;  // pointer to the simulation

      WaterGui() {
         // create a new simulation instance
         p_waterSim = new WaterSim(m_viewer, m_grid_res_x, m_grid_res_y,
                                   m_system_size_x, m_system_size_y,
                                   m_density, m_gravity, m_alpha,
                                   m_show_pressures, m_display_velocity_arrows);

         // set this simulation as the simulation that is running in our GUI
         setSimulation(p_waterSim);

         // start the GUI
         start();
      }

      /**
       * Update the Simulation class with the parameters input to the GUI
       */
      virtual void updateSimulationParameters() override {
         p_waterSim->setTimestep(m_dt);
         p_waterSim->updateParams(m_grid_res_x, m_grid_res_y,
                                  m_system_size_x, m_system_size_y,
                                  m_density, m_gravity, m_alpha, m_show_pressures, 
                                  m_display_velocity_arrows);
      };

      /**
       * Add parameter controls to the GUI
       */
      virtual void drawSimulationParameterMenu() override {
            ImGui::Checkbox("Export meshes", &m_export_meshes);
            ImGui::Checkbox("Show pressure field", &m_show_pressures);
            ImGui::InputDouble("Alpha", &m_alpha, 0, 0);
            ImGui::InputDouble("Timestep [s]", &m_dt, 0, 0);
            ImGui::InputDouble("Density [kg/m^3]", &m_density, 0, 0);
            ImGui::InputDouble("Gravity [m/s^2]", &m_gravity, 0, 0);
            ImGui::InputDouble("X Size [m]", &m_system_size_x, 0, 0);
            ImGui::InputDouble("Y Size [m]", &m_system_size_y, 0, 0);
            ImGui::InputInt("Grid resolution X", &m_grid_res_x, 0, 0);
            ImGui::InputInt("Grid resolution Y", &m_grid_res_y, 0, 0);
            ImGui::Checkbox("Display grid", &m_display_grid);
            ImGui::Checkbox("Display velocity arrows", &m_display_velocity_arrows);      
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
	        for (unsigned j = 0; j < nx; j++) {
	            for (unsigned i = 0; i < ny; i++) {
					double temp = (*p_waterSim).p_mac_grid->get_pressure(i, j);
					if(i == 0 && j == 0){
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
   // Create a new instance of the GUI for the simulation
   new WaterGui();

   return 0;
}

