#include "Simulation.h"

using namespace std;
/*
 * Example simulation that changes the colors of a cube.
 */
class DummySim : public Simulation {
   public:
    DummySim() : Simulation() { init(); }

    virtual void init() override {
        // create a cube on [-1,1]^3
        // vertices
        m_V.resize(8, 3);
        m_V << -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1,
            -1, 1, 1, 1, 1, 1;

        // faces
        m_F.resize(12, 3);
        m_F << 0, 2, 1, 2, 3, 1, 1, 3, 5, 3, 7, 5, 2, 6, 3, 6, 7, 3, 5, 7, 4, 7,
            6, 4, 4, 6, 0, 6, 2, 0, 0, 4, 1, 4, 5, 1;

        // face colors
        m_C.resize(12, 3);

        reset();
    }

    virtual void resetMembers() override {
        m_C.setZero();
        m_C.col(0).setOnes();
    }

    virtual void updateRenderGeometry() override {
        m_renderV = m_V;
        m_renderF = m_F;
        m_renderC = m_C;
    }

    virtual bool advance() override {
        // do next step of some color animation
        int speed = 60;
        int decColor = (m_step / speed) % 3;
        int incColor = (decColor + 1) % 3;

        for (int i = 0; i < m_C.rows(); i++) {
            m_C(i, decColor) = (m_C(i, decColor) * speed - 1) / speed;
            m_C(i, incColor) = (m_C(i, incColor) * speed + 1) / speed;
        }

        // advance step
        m_step++;
        return false;
    }

    virtual void renderRenderGeometry(
        igl::opengl::glfw::Viewer &viewer) override {
        viewer.data().set_mesh(m_renderV, m_renderF);
        viewer.data().set_colors(m_renderC);
    }

   private:
    Eigen::MatrixXd m_V;  // vertex positions
    Eigen::MatrixXi m_F;  // face indices
    Eigen::MatrixXd m_C;  // colors per face

    Eigen::MatrixXd m_renderV;  // vertex positions for rendering
    Eigen::MatrixXi m_renderF;  // face indices for rendering
    Eigen::MatrixXd m_renderC;  // colors per face for rendering
};
