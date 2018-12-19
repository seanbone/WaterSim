#ifndef WATERSIM_H
#define WATERSIM_H

#include "Simulation.h"
#include "FLIP.h"
#include "MeshExporter.h"

#include <sys/stat.h> // mkdir
#include <igl/opengl/glfw/Viewer.h>
#include <chrono> // accurate timings
#include <list> // std::list

/**
 * This class manages the water simulation.
 * It:
 * - Initializes a MAC grid in the right configuration
 * - Initializes the list of particles
 * - Initializes a FLIP object to handle updating the state
 * - Calls the update method on FLIP at each step
 * - Uses a MeshExporter object to export meshes at each step
 * - Handles rendering of the visualization (particles & grid)
 */
class WaterSim : public Simulation {

	private:
		// Pointer to IGL viewer used by Gui
		using viewer_t = igl::opengl::glfw::Viewer;
		viewer_t* const p_viewer;
		
		// Parameters of MAC grid
		int m_res_x, m_res_y, m_res_z; // no. of cells in X and Y directions
		double m_len_x, m_len_y, m_len_z; // size in [m] of full grid

		// FLIP solver parameters
		double m_fluid_density_;
		double m_gravity_mag_;
		double m_alpha_;

		// Max number of particles to display
		unsigned m_max_p_disp;

		// If we should export meshes for each time step
		bool m_export_meshes;

		// If we should show the grid in the visualization
		bool m_display_grid;

		// If we should randomize particle starting positions
		bool m_jitter_particles;

		// Number of particles in simulation
		unsigned int m_num_particles;
		

		/*** Members relating to visualization ***/

		// Index of the ViewerData object containing particles for rendering
		//  in viewer.data_list
		unsigned int m_particles_data_idx;
		// Index of the ViewerData object containing the mesh representing
		//  the grid MAC grid in viewer.data_list
		unsigned int m_grid_data_idx;
		Eigen::MatrixXd m_particles; // Particle positions for rendering, Nx3
		Eigen::MatrixXd m_particle_colors; // Particle colours, Nx3

		Eigen::MatrixXd m_renderV;  // vertex positions for rendering, Nx3
		Eigen::MatrixXi m_renderE;  // MAC grid edges for rendering, Nx2
		Eigen::MatrixXi m_renderF;  // face indices for rendering, Nx3
		Eigen::MatrixXd m_renderC;  // colors per face for rendering, Nx3
		Eigen::MatrixXd m_renderEC; // colors of edges of mac grid, Nx3
		
	public:
		// Vector of flags for fluid initialization
		std::vector<bool> is_fluid_;
		
		/*** Helper class members ***/

		//MAC grid data structure
		Mac3d* p_mac_grid;
		
		// List of Particles
		Particle* flip_particles;

		// FLIP simulator
		FLIP* p_flip;

		//MeshExporter
		MeshExporter* exp;
		
		/*** Public methods ***/

		/*
		 * Constructor
		 */
		WaterSim(viewer_t& viewer, const bool display_grid,
				 const int res_x, const int res_y, const int res_z,
				 const double len_x, const double len_y, const double len_z,
				 const double density, const double gravity,
				 const double alpha,
				 std::vector<bool> is_fluid, const bool jitter_particles,
				 bool export_meshes, unsigned max_p);

		/*
		 * Destructor
		 */
		~WaterSim() {
			delete p_mac_grid;
			delete p_flip;
			delete [] flip_particles;
		}
		
		/*
		 * Reset class variables to reset the simulation.
		 */
		virtual void resetMembers() override;

		/*
		 * Update simulation parameters. Requires a reset to take effect.
		 */
		void updateParams(const bool display_grid, const int res_x, const int res_y, const int res_z,
						  const double len_x, const double len_y, const double len_z,
						  const double density, const double gravity, const double alpha,
						  std::vector<bool> is_fluid, const bool jitter_particles,
						  bool export_meshes, unsigned max_p);

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
		
	private:

		/*
		 * Initialize the particle list according to the current settings.
		 * Note: does not delete previous particles!
		 */
		void initParticles();

		/*
		 * Initialize a new instance of the MAC grid
		 * Note: does not delete previous instance!
		 */
		void initMacGrid();
		
		/*
		 * Initialize a new instance of the Mesh Exporter
		 * Note: does not delete previous instance!
		 */
		void initMeshExp();

		/*
		 * Initialize a new instance of the FLIP simulator
		 * Note: does not delete previous instance!
		 */
		void initFLIP();

		/*
		 * Initialize a mesh to visualize the MAC grid
		 */
		void initMacViz();
};

#endif // WATERSIM_H
