#include "WaterSimGui.h"

#include <utility>

WaterSimGui::WaterSimGui(viewer_t &viewer, const SimConfig &cfg)
		: Simulation(), p_viewer(&viewer), m_watersim(cfg) {

	// Index of ViewerData instance dedicated to particles
	m_particles_data_idx = viewer.append_mesh();

	// Generate visualization mesh based on Mac3d
	initMacViz();

	// Index of ViewerData instance dedicated to MAC grid
	m_grid_data_idx = viewer.append_mesh();

	// Index of ViewerData instance dedicated to mesh
	m_mesh_data_idx = viewer.append_mesh();

	m_mesh_renderF = Eigen::MatrixXi::Zero(1, 3);
	m_mesh_renderV = Eigen::MatrixXd::Zero(1, 3);
	m_mesh_renderN = Eigen::MatrixXd::Zero(1, 3);

	// Update rendering geometry
	updateRenderGeometry();
}

/**
 * Reset class variables to reset the simulation.
 */
void WaterSimGui::resetMembers() {
	m_watersim.resetMembers();

	p_viewer->data_list[m_grid_data_idx].clear();
	initMacViz();
}


void WaterSimGui::updateParams(const SimConfig &cfg) {
	m_watersim.updateParams(cfg);

	setTimestep(cfg.getTimeStep());
}

void WaterSimGui::updateRenderGeometry() {
	// Copy particle positions from FLIP's data structure
	unsigned num_particles = m_watersim.getNumParticles();
	unsigned disp_particles = num_particles;
	unsigned particle_step = 1;

	// If there are too many particles to display, only
	// display a subset, using stride particle_set
	if (num_particles > m_watersim.m_cfg.getMaxParticlesDisplay()) {
		disp_particles = m_watersim.m_cfg.getMaxParticlesDisplay();
		particle_step = num_particles / disp_particles;
	}

	m_particles.resize(disp_particles, 3);
	for (unsigned i = 0, j = 0; j < disp_particles && i < num_particles; j++, i += particle_step) {
		m_particles(j, 0) = m_watersim.flip_particles->x[i];
		m_particles(j, 1) = m_watersim.flip_particles->y[i];
		m_particles(j, 2) = m_watersim.flip_particles->z[i];
	}

	m_particle_colors.resize(disp_particles, 3);
	m_particle_colors.setZero();
	m_particle_colors.col(2).setOnes();

	if (m_watersim.m_cfg.getDisplayMeshes()) {
		m_watersim.exp->get_mesh(m_mesh_renderV, m_mesh_renderF);
		igl::per_face_normals(m_mesh_renderV, m_mesh_renderF, m_mesh_renderN);
	} else {
		m_mesh_renderF = Eigen::MatrixXi::Zero(1, 3);
		m_mesh_renderV = Eigen::MatrixXd::Zero(1, 3);
		m_mesh_renderN = Eigen::MatrixXd::Zero(1, 3);
	}
}


bool WaterSimGui::advance() {
	return m_watersim.advance();
}


void WaterSimGui::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) {
	// Display MAC grid
	if (m_watersim.m_cfg.getDisplayGrid()) {
		viewer.data_list[m_grid_data_idx].set_edges(m_renderV, m_renderE, m_renderEC);
	} else {
		Eigen::MatrixXi renderE; // MAC grid edges for rendering, Nx2
		Eigen::MatrixXd renderEC; // colors of edges of mac grid, Nx3
		renderE.resize(0, 2);
		renderEC.resize(0, 3);
		viewer.data_list[m_grid_data_idx].set_edges(renderEC, renderE, renderEC);
	}

	// Display particles
	viewer.data_list[m_particles_data_idx].set_points(m_particles, m_particle_colors);
	viewer.data_list[m_particles_data_idx].point_size = 5;

	// Display mesh
	viewer.data_list[m_mesh_data_idx].clear();
	if (m_watersim.m_cfg.getDisplayMeshes() && m_watersim.getStep() > 0
	    && m_mesh_renderF.cols() == 3) {
		viewer.data_list[m_mesh_data_idx].set_mesh(m_mesh_renderV, m_mesh_renderF);
		viewer.data_list[m_mesh_data_idx].set_normals(m_mesh_renderN);
	}
}


void WaterSimGui::initMacViz() {

	//if (!m_watersim.m_cfg.getDisplayGrid())
	//    return;

	unsigned nx = m_watersim.p_mac_grid->get_num_cells_x();
	unsigned ny = m_watersim.p_mac_grid->get_num_cells_y();
	unsigned nz = m_watersim.p_mac_grid->get_num_cells_z();
	double sx = m_watersim.p_mac_grid->get_cell_sizex();
	double sy = m_watersim.p_mac_grid->get_cell_sizey();
	double sz = m_watersim.p_mac_grid->get_cell_sizez();

	unsigned num_vertices = (nx + 1) * (ny + 1) * (nz + 1);
	unsigned num_edges = (nx + ny + 2) * (nz + 1) + (nx + 1) * (ny + 1);

	// Calculate grid vertex coordinates for rendering
	m_renderV.resize(num_vertices, 3);
	unsigned i = 0;
	for (unsigned z = 0; z <= nz; z++) {
		for (unsigned y = 0; y <= ny; y++) {
			for (unsigned x = 0; x <= nx; x++) {
				// Vertex:
				double cx = x * sx;
				double cy = y * sy;
				double cz = z * sz;
				m_renderV.row(i) << cx - sx / 2., cy - sy / 2., cz - sz / 2.;

				// Increment index
				i++;
			}
		}
	}

	// Calculate edges of MAC grid for rendering
	m_renderE.resize(num_edges, 2);
	m_renderE.setZero();
	m_renderEC.resize(num_edges, 3);
	m_renderEC.setZero();
	i = 0;
	// Edges parallel to x axis <=> normal to yz plane
	for (unsigned z = 0; z <= nz; z++) {
		for (unsigned y = 0; y <= ny; y++, i++) {
			auto a = y * (nx + 1) + z * (nx + 1) * (ny + 1);
			m_renderE.row(i) << a, a + nx;
		}
	}
	// Edges parallel to y axis <=> normal to xz plane
	for (unsigned z = 0; z <= nz; z++) {
		for (unsigned x = 0; x <= nx; x++, i++) {
			auto a = x + z * (nx + 1) * (ny + 1);
			m_renderE.row(i) << a, a + ny * (nx + 1);
		}
	}
	// Edges parallel to z axis <=> normal to xy plane
	for (unsigned y = 0; y <= ny; y++) {
		for (unsigned x = 0; x <= nx; x++, i++) {
			auto a = x + y * (nx + 1);
			m_renderE.row(i) << a, a + nz * (ny + 1) * (nx + 1);
		}
	}
}

