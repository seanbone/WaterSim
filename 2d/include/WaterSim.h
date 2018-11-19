#ifndef WATERSIM_H
#define WATERSIM_H

#include "Simulation.h"
#include "FLIP.h"

/**
 * This class manages the water simulation.
 * It:
 * - Initializes a MAC grid in the right configuration
 * - Initializes the list of particles
 * - Initializes a FLIP object to handle updating the state
 * - Calls the update method on FLIP at each step
 * - Uses an Exporter object to export meshes at each step
 * - Handles rendering of the visualization (particles & grid)
 */
class WaterSim : public Simulation {

	private:
		
		// MAC grid data structure
		Mac2d* p_mac_grid;

		// List of Particles
		Particle* flip_particles;

		// FLIP simulator
		FLIP* p_flip;

		// Index of the ViewerData object containing particles for rendering
		//  in viewer.data_list
		unsigned int m_particles_data_idx;
		// Index of the ViewerData object containing the mesh representing
		//  the grid MAC grid in viewer.data_list
		unsigned int m_grid_data_idx;

		unsigned int m_num_particles;
		Eigen::MatrixXd m_particles; // Particle positions for rendering, Nx3
		Eigen::MatrixXd m_particle_colors; // Particle colours, Nx3

		Eigen::MatrixXd m_renderV;  // vertex positions for rendering, Nx3
		Eigen::MatrixXi m_renderE;  // MAC grid edges for rendering, Nx2
		Eigen::MatrixXi m_renderF;  // face indices for rendering, Nx3
		Eigen::MatrixXd m_renderC;  // colors per face for rendering, Nx3
		Eigen::MatrixXd m_renderEC; // colors of edges of mac grid, Nx3


	using viewer_t = igl::opengl::glfw::Viewer;
	public:

		WaterSim(viewer_t& viewer, const int res_x, const int res_y, const double len_x, const double len_y);

		~WaterSim() {
			delete p_mac_grid;
			delete p_flip;
			delete [] flip_particles;
		}
		
		using Simulation::init;

		virtual void init(viewer_t& viewer);

		/*
		 * Reset class variables to reset the simulation.
		 */
		virtual void resetMembers() override;

		/*
		 * Update the rendering data structures. This method will be called in
		 * alternation with advance(). This method blocks rendering in the
		 * viewer, so do *not* do extensive computation here (leave it to
		 * advance()). This method need not be thread-safe with
		 * renderRenderGeometry(): mutual exclusion is guaranteed by Simulator.
		 */
		virtual void updateRenderGeometry() override;
		
		/*
		 * Performs one simulation step of length m_dt. This method *must* be
		 * thread-safe with respect to renderRenderGeometry() (easiest is to not
		 * touch any rendering data structures at all). You have to update the time
		 * variables at the end of each step if they are necessary for your
		 * simulation.
		 */
		virtual bool advance() override;
		
		/*
		 * Perform any actual rendering here. This method *must* be thread-safe with
		 * respect to advance(). This method runs in the same thread as the
		 * viewer and blocks user IO, so there really should not be any extensive
		 * computation here or the UI will lag/become unresponsive (the whole reason
		 * the simulation itself is in its own thread.)
		 */
		virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) override;
		
};

#endif // WATERSIM_H
