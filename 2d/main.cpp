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

      double m_system_size_x = 10; // X dimension of system in m
      double m_system_size_y = 10; // Y dimension of system in m

      int m_grid_res_x = 100; // Number of cells on X axis
      int m_grid_res_y = 100; // Number of cells on Y axis

      // Other members

   public:
      WaterSim *p_waterSim = NULL;  // pointer to the simulation

      WaterGui() {
         // create a new simulation instance
         p_waterSim = new WaterSim(m_viewer);

         // set this simulation as the simulation that is running in our GUI
         setSimulation(p_waterSim);

         // start the GUI
         start();
      }

      /**
       * Update the Simulation class with the parameters input to the GUI
       */
      virtual void updateSimulationParameters() override {
         std::cout << "\n*********\n Simulation Parameters:\n\n";
         std::cout << "Export meshes:\t" << m_export_meshes << std::endl;
         std::cout << "Export freq.:\t" << m_export_fps << std::endl;
         std::cout << "X Size [m]:\t" << m_system_size_x << std::endl;
         std::cout << "Y Size [m]:\t" << m_system_size_y << std::endl;
         std::cout << "Resolution X:\t" << m_grid_res_x << std::endl;
         std::cout << "Resolution Y:\t" << m_grid_res_y << std::endl;
         std::cout << "\n*********\n";
      };

      /**
       * Add parameter controls to the GUI
       */
      virtual void drawSimulationParameterMenu() override {
//            if (ImGui::Button("Export Trajectories", ImVec2(-1, 0))) {
//               exportTrajectories();
//            }
//            ImGui::SliderAngle("Angle", &m_angle, -180.0f, 180.0f);
//            ImGui::InputFloat("Force", &m_force, 0, 0);
//            ImGui::InputFloat("Mass", &m_mass, 0, 0);
//            ImGui::InputFloat("dt", &m_dt, 0, 0);
//            ImGui::Combo("Integrator", &m_selected_integrator, m_integrators.data(),
//                         m_integrators.size());
            ImGui::Checkbox("Export meshes", &m_export_meshes);
            ImGui::InputInt("Export FPS", &m_export_fps, 0, 0);
            ImGui::InputDouble("X Size [m]", &m_system_size_x, 0, 0);
            ImGui::InputDouble("Y Size [m]", &m_system_size_y, 0, 0);
            ImGui::InputInt("Grid resolution X", &m_grid_res_x, 0, 0);
            ImGui::InputInt("Grid resolution Y", &m_grid_res_y, 0, 0);
      }
};

int main() { //int argc, char *argv[]) {
   // Create a new instance of the GUI for the simulation
   new WaterGui();

   return 0;
}
