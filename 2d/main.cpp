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
      int m_export_fps = 30;

      double m_system_size_x = 30; // X dimension of system in m
      double m_system_size_y = 30; // Y dimension of system in m

      int m_grid_res_x = 15; // Number of cells on X axis
      int m_grid_res_y = 15; // Number of cells on Y axis
      
      // TODO: proper input
      double m_dt = 0.005 * sqrt((m_grid_res_x + m_grid_res_y) * 0.5);

      // Other members
      bool m_display_grid = true;
      int m_cells_x;
      int m_cells_y;

   public:
      WaterSim *p_waterSim = NULL;  // pointer to the simulation

      WaterGui() {
         // create a new simulation instance
         //WaterSim(viewer_t& viewer, const int res_x, const int res_y, const double len_x, const double len_y);
         p_waterSim = new WaterSim(m_viewer, m_grid_res_x, m_grid_res_y, m_system_size_x, m_system_size_y);

         // set this simulation as the simulation that is running in our GUI
         setSimulation(p_waterSim);

         // start the GUI
         start();
      }

      /**
       * Update the Simulation class with the parameters input to the GUI
       */
      virtual void updateSimulationParameters() override {
         // TODO: copy GUI inputs to WaterSim
         std::cout << "\n*********\n Simulation Parameters:\n\n";
         std::cout << "Export meshes:\t" << m_export_meshes << std::endl;
         std::cout << "Export freq.:\t" << m_export_fps << std::endl;
         std::cout << "X Size [m]:\t" << m_system_size_x << std::endl;
         std::cout << "Y Size [m]:\t" << m_system_size_y << std::endl;
         std::cout << "Resolution X:\t" << m_grid_res_x << std::endl;
         std::cout << "Resolution Y:\t" << m_grid_res_y << std::endl;
         std::cout << "\n*********\n";
         p_waterSim->setTimestep(m_dt);
      };

      /**
       * Add parameter controls to the GUI
       */
      virtual void drawSimulationParameterMenu() override {
            ImGui::Checkbox("Export meshes", &m_export_meshes);
            ImGui::InputInt("Export FPS", &m_export_fps, 0, 0);
            ImGui::InputDouble("X Size [m]", &m_system_size_x, 0, 0);
            ImGui::InputDouble("Y Size [m]", &m_system_size_y, 0, 0);
            ImGui::InputInt("Grid resolution X", &m_grid_res_x, 0, 0);
            ImGui::InputInt("Grid resolution Y", &m_grid_res_y, 0, 0);
            ImGui::Checkbox("Display grid", &m_display_grid);
      }
};

int main() { //int argc, char *argv[]) {
   // Create a new instance of the GUI for the simulation
   new WaterGui();

   return 0;
}

