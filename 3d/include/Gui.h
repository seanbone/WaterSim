#ifndef GUI_H
#define GUI_H

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include "Arrow.h"
#include "Simulator.h"

/*
 * Base class to open a window with a simple GUI for simulation. Inherit from
 * this class to perform a simulation and customize the menu.
 */
class Gui {
   public:
    Gui() {}

    ~Gui() {}

    /*
     * Set simulation to be performed in the simulator.
     */
    void setSimulation(Simulation *sim);

    /*
     * Initialize all the necessary data structures and callbacks and open the
     * window.
     */
    void start();

    /*
     * Call the setters for your simulation in this method to update the
     * parameters entered in the GUI. This method is called before the
     * simulation is started for the first time.
     */
    virtual void updateSimulationParameters() = 0;

    /*
     * Clear some custom datastructures/visualizations when the user requests
     * it.
     */
    virtual void clearSimulation() {}

    /*
     * Callback to enable custom shortcuts.
     */
    virtual bool childKeyCallback(igl::opengl::glfw::Viewer &viewer,
                                  unsigned int key, int modifiers) {
        return false;
    }

    /*
     * Setup your own (additional) ImGUI components in this method.
     */
    virtual void drawSimulationParameterMenu() {}

    /*
     * Setup your own (additional) ImGUI debugging output.
     */
    virtual void drawSimulationStats() {}

#pragma region ArrowInterface
    /*
     * Create and add an arrow to be displayed in the GUI. Returns the index of
     * the drawn arrow (keep it if you want to delete this arrow later).
     */
    int addArrow(const Eigen::Vector3d &start, const Eigen::Vector3d &end);
    int addArrow(const Eigen::Vector3d &start, const Eigen::Vector3d &end,
                 const Eigen::Vector3d &color);
    /*
     * Delete arrow at given index.
     */
    void removeArrow(size_t index);
#pragma endregion ArrowInterface

    /*
     * Show the standard basis axis (x, y, z)
     */
    void showAxes(bool show_axes);

    /*
     * Callback to show a different vector for a clicked vertex than the normal
     */
    std::function<void(int clickedVertexIndex, int clickedObjectIndex,
                       Eigen::Vector3d &pos, Eigen::Vector3d &dir)>
        callback_clicked_vertex = nullptr;

   private:
    void drawArrow(const Arrow &arrow);

    void showVertexArrow();

    void toggleSimulation();

    void singleStep();

    void exportRecording();

    void clearScreen();

    bool keyCallback(igl::opengl::glfw::Viewer &viewer, unsigned int key,
                     int modifiers);

    bool keyReleasedCallback(igl::opengl::glfw::Viewer &viewer,
                             unsigned int key, int modifiers);

    void drawMenuWindow(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    bool drawMenu(igl::opengl::glfw::Viewer &viewer,
                  igl::opengl::glfw::imgui::ImGuiMenu &menu);

    bool drawCallback(igl::opengl::glfw::Viewer &viewer);

    bool scrollCallback(igl::opengl::glfw::Viewer &viewer, float delta_y);

    bool mouseCallback(igl::opengl::glfw::Viewer &viewer,
                       igl::opengl::glfw::imgui::ImGuiMenu &menu, int button,
                       int modifier);

    Simulator *p_simulator = NULL;
    bool m_request_clear = false;
    int m_simSpeed = 60;
    bool m_fastForward = false;

    int m_clickedVertex = -1;     // index of clicked vertex
    int m_clickedObject = -1;     // id of clicked object
    int m_clickedArrow = -1;      // index of arrow of clicked vertex
    std::vector<Arrow> m_arrows;  // data structure to store all the arrows
                                  // to be rendered
    unsigned long m_numArrows = 0;         // increasing counter for arrows

    int m_axesID = -1;  // (lowest) id of the 3 base axes
    bool m_showAxes = true;

    bool m_showStats = true;
    double m_timerAverage = 0;  // running average of execution time of
                                // one iteration of the simulation
    int m_maxSteps = -1;

    int m_numRecords = 100;  // number of records to keep

   protected:
    igl::opengl::glfw::Viewer m_viewer;

    virtual void resetSimulation();

};

#endif
