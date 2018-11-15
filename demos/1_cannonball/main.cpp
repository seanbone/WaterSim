#include <igl/writeOFF.h>
#include "CannonBallSim.h"
#include "Gui.h"

/*
 * GUI for the cannonball simulation. This time we need additional paramters,
 * e.g. which integrator to use for the simulation and the force applied to the
 * cannonball, and we also add some more visualizations (trajectories).
 */
class CannonBallGui : public Gui {
   public:
    // simulation parameters
    float m_angle = 1.047f;
    float m_force = 100.0f;
    float m_dt = 1e-2;
    float m_mass = 1.0;
    int m_log_frequency = 30;

    CannonBallSim *p_cannonBallSim = NULL;

    const vector<char const *> m_integrators = {"Analytic", "Explicit Euler",
                                                "Symplectic Euler"};
    int m_selected_integrator = 0;

    CannonBallGui() {
        // create a new cannonball simulation, set it in the GUI,
        // and start the GUI
        p_cannonBallSim = new CannonBallSim();
        setSimulation(p_cannonBallSim);

        // show vertex velocity instead of normal
        callback_clicked_vertex = [&](int clickedVertexIndex,
                                      int clickedObjectIndex,
                                      Eigen::Vector3d &pos,
                                      Eigen::Vector3d &dir) {
            RigidObject &o = p_cannonBallSim->getObjects()[clickedObjectIndex];
            pos = o.getVertexPosition(clickedVertexIndex);
            dir = o.getVelocity(pos);
        };
        start();
    }

    virtual void updateSimulationParameters() override {
        // change all parameters of the simulation to the values that are set in
        // the GUI
        p_cannonBallSim->setForce(m_force);
        p_cannonBallSim->setAngle(m_angle);
        p_cannonBallSim->setTimestep(m_dt);
        p_cannonBallSim->setMass(m_mass);
        p_cannonBallSim->setMethod(m_selected_integrator);
        p_cannonBallSim->setLogFrequency(m_log_frequency);
    }

    virtual void clearSimulation() override {
        p_cannonBallSim->clearTrajectories();
    }

    /*
     * Writes each trajectory to an individual off-file.
     */
    void exportTrajectories() {
        Eigen::MatrixX3d mat;
        for (int i = 0; i < p_cannonBallSim->getNumTrajectories(); i++) {
            string filename = "trajectory" + to_string(i) + ".off";
            p_cannonBallSim->getTrajectories(i, mat);
            if (mat.rows() <= 1) {
                continue;
            }
            if (igl::writeOFF(filename, mat, Eigen::MatrixXi())) {
                cout << "Wrote trajectory to " << filename << endl;
            } else {
                cout << "Failed to write trajectory to " << filename << endl;
            }
        }
    }

    virtual bool childKeyCallback(igl::opengl::glfw::Viewer &viewer,
                                  unsigned int key, int modifiers) override {
        switch (key) {
            case 'e':
            case 'E':
                exportTrajectories();
                return true;
            // cicle through different integrators
            case '>':
                m_selected_integrator++;
                m_selected_integrator %= m_integrators.size();
                return true;
            case '<':
                m_selected_integrator--;
                m_selected_integrator =
                    (m_integrators.size() + m_selected_integrator) %
                    m_integrators.size();
                return true;
        }
        return false;
    }

    virtual void drawSimulationParameterMenu() override {
        if (ImGui::Button("Export Trajectories", ImVec2(-1, 0))) {
            exportTrajectories();
        }
        ImGui::SliderAngle("Angle", &m_angle, -180.0f, 180.0f);
        ImGui::InputFloat("Force", &m_force, 0, 0);
        ImGui::InputFloat("Mass", &m_mass, 0, 0);
        ImGui::InputFloat("dt", &m_dt, 0, 0);
        ImGui::Combo("Integrator", &m_selected_integrator, m_integrators.data(),
                     m_integrators.size());
        ImGui::InputInt("Log Frequency", &m_log_frequency, 0, 0);
    }
};

int main(int argc, char *argv[]) {
    // create a new instance of the GUI for the cannonball simulation
    new CannonBallGui();

    return 0;
}