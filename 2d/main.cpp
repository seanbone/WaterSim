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
   public:
      WaterSim *p_waterSim = NULL;  // pointer to the simulation

      WaterGui() {
         // create a new simulation instance
         p_waterSim = new WaterSim();

         // set this simulation as the simulation that is running in our GUI
         setSimulation(p_waterSim);

         // start the GUI
         start();
      }

      virtual void updateSimulationParameters() override{
         // TODO: add simulation parameters to GUI
      };
};

int main() { //int argc, char *argv[]) {
   // Create a new instance of the GUI for the simulation
   new WaterGui();

   return 0;
}
