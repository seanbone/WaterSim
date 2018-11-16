#include <igl/writeOFF.h>
#include <thread>
#include "DummySim.h"
#include "Gui.h"
#include "Simulator.h"

/*
 * This class is a GUI for our dummy simulation. It extends the basic GUI
 * defined in Gui.h. We could add more controls and visuals here, but we don't
 * need any additional functionality for this dummy simulation.
 */
class DummyGui : public Gui {
   public:
    DummySim *p_dummySim = NULL;  // pointer to the dummy simulation

    DummyGui() {
        // create a new dummy simulation
        p_dummySim = new DummySim();

        // set this simulation as the simulation that is running in our GUI
        setSimulation(p_dummySim);

        // start the GUI
        start();
    }

    virtual void updateSimulationParameters() override{
        // We don't have any simulation parameters to update periodically so we
        // don't need to do anything here
    };
};

int main(int argc, char *argv[]) {
    // create a new instance of the GUI for the dummy simulation
    new DummyGui();

    return 0;
}