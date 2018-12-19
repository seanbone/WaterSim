#ifndef WATERSIM_H
#define WATERSIM_H

#include "Simulation.h"
#include "FLIP.h"

#include <sys/stat.h> // mkdir
#include <igl/opengl/glfw/Viewer.h>
#include <igl/png/writePNG.h>

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
		using viewer_t = igl::opengl::glfw::Viewer;

		// Pointer to IGL viewer used by Gui
		viewer_t* const p_viewer;
		
		// MAC grid data structure
		int m_res_x, m_res_y; // no. of cells in X and Y directions
		double m_len_x, m_len_y; // size in [m] of full grid

		// Other params
		double m_fluid_density_;
		double m_gravity_mag_;
		double m_alpha_;
		bool m_show_pressures;
		bool m_show_velocity_arrows;

		// PNG export params
		bool m_export_png_;
		size_t m_png_sx_;
		size_t m_png_sy_;
		std::string m_png_dirname_ = "../PNG_out";
		unsigned m_png_num_ = 0;
		unsigned m_max_pngs_;

		// List of Particles
		Particle* flip_particles;
		std::vector<bool> is_fluid_;

		// FLIP simulator
		FLIP* p_flip;

		// Index of the ViewerData object containing particles for rendering
		//  in viewer.data_list
		unsigned int m_particles_data_idx;
		// Index of the ViewerData object containing the mesh representing
		//  the grid MAC grid in viewer.data_list
		unsigned int m_grid_data_idx;
		unsigned int m_velocity_u_idx;
		unsigned int m_velocity_v_idx;
		

		unsigned int m_num_particles;
		bool m_jitter_particles;
		Eigen::MatrixXd m_particles; // Particle positions for rendering, Nx3
		Eigen::MatrixXd m_particle_colors; // Particle colours, Nx3

		Eigen::MatrixXd m_renderV;  // vertex positions for rendering, Nx3
		Eigen::MatrixXi m_renderE;  // MAC grid edges for rendering, Nx2
		Eigen::MatrixXi m_renderF;  // face indices for rendering, Nx3
		Eigen::MatrixXd m_renderC;  // colors per face for rendering, Nx3
		Eigen::MatrixXd m_renderEC; // colors of edges of mac grid, Nx3
		
		Eigen::MatrixXd m_render_velocity_u_V;
		Eigen::MatrixXi m_render_velocity_u_E;
		Eigen::MatrixXd m_render_velocity_v_V;
		Eigen::MatrixXi m_render_velocity_v_E;

	public:
		//MAC grid data structure
		Mac2d* p_mac_grid;
		
		WaterSim(viewer_t& viewer, 
				 const int res_x, const int res_y,
				 const double len_x, const double len_y,
				 const double density, const double gravity,
				 const double alpha,
				 const bool show_pressures, const bool show_velocity_arrows,
				 std::vector<bool> is_fluid, const bool jitter_particles,
				 bool export_png, int png_sx, int png_sy, int max_pngs);

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
		void updateParams(const int res_x, const int res_y, 
						  const double len_x, const double len_y,
						  const double density, const double gravity, const double alpha,
						  const bool show_pressures, const bool show_velocity_arrows,
						  std::vector<bool> is_fluid, const bool jitter_particles,
						  bool export_png, int png_sx, int png_sy, int max_pngs);

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
		 * Initialize a new instance of the FLIP simulator
		 * Note: does not delete previous instance!
		 */
		void initFLIP();

		/*
		 * Initialize a mesh to visualize the MAC grid
		 */
		void initMacViz();

		/*
		 * Export current viewport as PNG
		 */
		void exportPNG(igl::opengl::glfw::Viewer &viewer);
};

#endif // WATERSIM_H
