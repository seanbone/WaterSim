#ifndef WATERSIM_H
#define WATERSIM_H

#include "Simulation.h"
#include "FLIP.h"
//~ #include "MarchCubes.h"

class WaterSim : public Simulation {
	
	// Initialize the grid and the particle positions and velocities
	// for the FLIP algorithm
	//~ TODO: connect to FLIP algorithm
	
	using viewer_t = igl::opengl::glfw::Viewer;
	public:
		WaterSim(viewer_t& viewer) : Simulation() { init(viewer); }
		
		using Simulation::init;

		virtual void init(viewer_t& viewer) {
			// Create a plane spanning [-1,-1,0]x[1,1,0]
			// Vertices
			m_V.resize(9, 3);
			m_V << -1, -1, 0,
					0, -1, 0, 
					1, -1, 0, 
				   -1,  0, 0, 
					0,  0, 0, 
					1,  0, 0,
				   -1,  1, 0,
					0,  1, 0, 
					1,  1, 0;
			
			// Faces
			m_F.resize(8, 3);
			m_F << 0, 1, 3,
				   3, 1, 4,
				   1, 2, 4,
				   4, 2, 5,
				   3, 4, 6,
				   6, 4, 7,
				   4, 5, 7,
				   7, 5, 8;
			
			// face colors
			m_C.resize(8, 3);

			// Create some dummy particles for rendering
			m_particles.resize(9, 3);
			m_particles << 0, -.5, 0,
					0, 0, 0, 
					0.5, -.5, 0, 
				   -0.5,  0, 0, 
					0.5,  0, 0, 
					0.5, .5, 0,
				   -0.5,-.5, 0,
					  0, .5, 0, 
				   -0.5, .5, 0;

			m_num_particles = m_particles.rows();

			m_particle_colors.resize(m_num_particles, 3);
			m_particle_colors.setZero();
			m_particle_colors.col(2).setOnes();

			m_particles_data_idx = viewer.append_mesh();
			
			reset();
		}
		
		/*
		 * Reset class variables to reset the simulation.
		 */
		virtual void resetMembers() override {
			m_C.setZero();
			m_C.col(0).setOnes();
			m_C.col(1).setOnes();
			m_particle_colors.setZero();
			m_particle_colors.col(2).setOnes();
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
			// TODO: import particle positions from simulation (FLIP)
		}
		
		/*
		 * Performs one simulation step of length m_dt. This method *must* be
		 * thread-safe with respect to renderRenderGeometry() (easiest is to not
		 * touch any rendering data structures at all). You have to update the time
		 * variables at the end of each step if they are necessary for your
		 * simulation.
		 */
		virtual bool advance() override {
//			// do next step of some color animation
//			int speed = 60;
//			int decColor = (m_step / speed) % 3;
//			int incColor = (decColor + 1) % 3;
//			
//			for (int i = 0; i < m_C.rows(); i++) {
//				m_C(i, decColor) = (m_C(i, decColor) * speed - 1) / speed;
//				m_C(i, incColor) = (m_C(i, incColor) * speed + 1) / speed;
//			}


			Eigen::MatrixXd delta;
			delta.resize(m_num_particles, 3);
			delta.col(0).setRandom();
			delta.col(1).setRandom();
			delta.col(2).setZero();

			m_particles += 0.01 * delta;
			
			// advance step
			m_step++;
			m_time += m_dt;
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
			viewer.data_list[m_particles_data_idx].set_points(m_particles, m_particle_colors);
			viewer.data_list[m_particles_data_idx].point_size = 5;
		}
		
	private:
		
		// Index of the ViewerData object containing particles for rendering
		//  in viewer.data_list
		unsigned int m_particles_data_idx;

		unsigned int m_num_particles;
		Eigen::MatrixXd m_particles; // Particle positions for rendering, Nx3
		Eigen::MatrixXd m_particle_colors; // Particle colours, Nx3

		Eigen::MatrixXd m_V;  // vertex positions, Nx3 for N vertices
		Eigen::MatrixXi m_F;  // face indices
		Eigen::MatrixXd m_C;  // colors per face
		
		Eigen::MatrixXd m_renderV;  // vertex positions for rendering
		Eigen::MatrixXi m_renderF;  // face indices for rendering
		Eigen::MatrixXd m_renderC;  // colors per face for rendering
};

#endif // WATERSIM_H
