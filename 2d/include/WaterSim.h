#ifndef WATERSIM_H
#define WATERSIM_H

#include "Simulation.h"
#include "FLIP.h"
//~ #include "MarchCubes.h"

class WaterSim : public Simulation {
	
	// Initialize the grid and the particle positions and velocities
	// for the FLIP algorithm
	//~ TODO: use the FLIP constructor in a 
	//~ 	  modified (inheritance) Simulation::init()?
	// WaterSim currently is a copy of DummySim from 0_dummy
	public:
		WaterSim() : Simulation() { init(); }
		
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
		
		/*
		 * Reset class variables to reset the simulation.
		 */
		virtual void resetMembers() override {
			m_C.setZero();
			m_C.col(0).setOnes();
		}
		
		/*
		 * Update the rendering data structures. This method will be called in
		 * alternation with advance(). This method blocks rendering in the
		 * viewer, so do *not* do extensive computation here (leave it to
		 * advance()).
		 */
		virtual void updateRenderGeometry() override {
			m_renderV = m_V;
			m_renderF = m_F;
			m_renderC = m_C;
		}
		
		/*
		 * Performs one simulation step of length m_dt. This method *must* be
		 * thread-safe with respect to renderRenderGeometry() (easiest is to not
		 * touch any rendering data structures at all). You have to update the time
		 * variables at the end of each step if they are necessary for your
		 * simulation.
		 */
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
		
		/*
		 * Perform any actual rendering here. This method *must* be thread-safe with
		 * respect to advance(). This method runs in the same thread as the
		 * viewer and blocks user IO, so there really should not be any extensive
		 * computation here or the UI will lag/become unresponsive (the whole reason
		 * the simulation itself is in its own thread.)
		 */
		virtual void renderRenderGeometry(
			igl::opengl::glfw::Viewer &viewer) override {
			viewer.data().set_mesh(m_renderV, m_renderF);
			viewer.data().set_colors(m_renderC);
		}
		
	private:
		Eigen::MatrixXd m_V;  // vertex positions, Nx3 for N vertices
		Eigen::MatrixXi m_F;  // face indices
		Eigen::MatrixXd m_C;  // colors per face
		
		Eigen::MatrixXd m_renderV;  // vertex positions for rendering
		Eigen::MatrixXi m_renderF;  // face indices for rendering
		Eigen::MatrixXd m_renderC;  // colors per face for rendering
};

#endif
