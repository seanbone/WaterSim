#include "WaterSim.h"


WaterSim::WaterSim(WaterSim::viewer_t& viewer,
        const int res_x, const int res_y,
        const double len_x, const double len_y,
        const double density, const double gravity, 
        const double alpha, const bool show_pressures,
        const bool show_velocity_arrows)
        : Simulation(), p_viewer(&viewer), m_res_x(res_x), m_res_y(res_y),
          m_len_x(len_y), m_len_y(len_y), m_fluid_density_(density),
          m_gravity_mag_(gravity), m_alpha_(alpha),
          m_show_pressures(show_pressures),
          m_show_velocity_arrows(show_velocity_arrows) {

    init();
}


void WaterSim::init() {
    WaterSim::viewer_t& viewer = *p_viewer;

    // Initialize MAC grid
    initMacGrid();

    // Initialize particles
    initParticles();

    // Initialize FLIP object
    initFLIP();
    
    // Index of ViewerData instance dedicated to particles
    m_particles_data_idx = viewer.append_mesh();
    
    // Generate visualization mesh based on Mac2d
    initMacViz();

    // Index of ViewerData instance dedicated to MAC grid
    m_grid_data_idx = viewer.append_mesh();
    
    // Update rendering geometry
    updateRenderGeometry();
}


/*
 * Reset class variables to reset the simulation.
 */
void WaterSim::resetMembers() {
    // MAC grid
    delete p_mac_grid;
    initMacGrid();

    p_viewer->data_list[m_grid_data_idx].clear();
    initMacViz();

    // Particles
    delete [] flip_particles;
    initParticles();

    // FLIP simulator
    delete p_flip;
    initFLIP();
}



void WaterSim::updateParams(const int res_x, const int res_y, 
				  const double len_x, const double len_y,
                  const double density, const double gravity, const double alpha,
                  const bool show_pressures, const bool show_velocity_arrows) {
    m_res_x = res_x;
    m_res_y = res_y;
    m_len_x = len_x;
    m_len_y = len_y;
    m_fluid_density_ = density;
    m_gravity_mag_ = gravity;
    m_alpha_ = alpha;
    m_show_pressures = show_pressures;
    m_show_velocity_arrows = show_velocity_arrows;
    std::cout << "\nParams updated\n";
}


void WaterSim::updateRenderGeometry() {
    // Copy particle positions from FLIP's data structure
    m_particles.resize(m_num_particles, 3);
    for (unsigned i = 0; i < m_num_particles; i++) {
        m_particles.row(i) = flip_particles[i].get_position();//.transpose();
    }

    m_particle_colors.resize(m_num_particles, 3);
    m_particle_colors.setZero();
    m_particle_colors.col(2).setOnes();

    if (m_show_pressures) {
        // Render pressure as a scalar field
        unsigned nx = p_mac_grid->get_num_cells_x();
        unsigned ny = p_mac_grid->get_num_cells_y();
        Eigen::VectorXd pressures(2*nx*ny);
        for (unsigned j = 0; j < nx; j++) {
            for (unsigned i = 0; i < ny; i++) {
                pressures(2*(i + j*ny)) = p_mac_grid->get_pressure(i, j);
                pressures(2*(i + j*ny) + 1) = p_mac_grid->get_pressure(i, j);
            }
        }

        igl::jet(pressures, true, m_renderC);
    }

    //~ std::cout << "\n*************\n";
    //~ std::cout << "Pressure at (7, 0): " << p_mac_grid->get_pressure(7, 0) << std::endl;
    //~ std::cout << "Pressure at (7, 1): " << p_mac_grid->get_pressure(7, 1) << std::endl;
    //~ std::cout << "U velocity at (6.5, 1): " << p_mac_grid->get_u(7, 1) << std::endl;
    //~ std::cout << "V velocity at (6.5, 1): " << p_mac_grid->get_v(7, 1) << std::endl;
    //~ std::cout << "U* velocity at (6.5, 1): " << p_mac_grid->get_u_star(7, 1) << std::endl;
    //~ std::cout << "V* velocity at (6.5, 1): " << p_mac_grid->get_v_star(7, 1) << std::endl;
    //~ std::cout << "X of particle 0: " << flip_particles->get_position()(0) << std::endl;
    //~ std::cout << "Y of particle 0: " << flip_particles->get_position()(1) << std::endl;
    //~ std::cout << "U of particle 0: " << flip_particles->get_velocity()(0) << std::endl;
    //~ std::cout << "V of particle 0: " << flip_particles->get_velocity()(1) << std::endl;
    //~ std::cout << "\n*************\n";
}


bool WaterSim::advance() {
    // Perform a FLIP step
    p_flip->step_FLIP(m_dt, m_step);
    
    // advance step
    m_step++;
    m_time += m_dt;
    return false;
}


void WaterSim::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) {
    viewer.data_list[m_grid_data_idx].set_mesh(m_renderV, m_renderF);
    viewer.data_list[m_grid_data_idx].set_colors(m_renderC);
    viewer.data_list[m_grid_data_idx].set_edges(m_renderV, m_renderE, m_renderEC);

    viewer.data_list[m_particles_data_idx].set_points(m_particles, m_particle_colors);
    viewer.data_list[m_particles_data_idx].point_size = 5;
    
    if(m_show_velocity_arrows){
		unsigned nx = p_mac_grid->get_num_cells_x();
        unsigned ny = p_mac_grid->get_num_cells_y();
        double dx = p_mac_grid->get_cell_sizex();
        double dy = p_mac_grid->get_cell_sizey();
        for (unsigned j = 0; j < nx + 1; j++) {
            for (unsigned i = 0; i < ny; i++) {
                Eigen::RowVector3d start((i-0.5)*dx, j*dy, 0);
                double temp = p_mac_grid->get_u(i,j);
                temp *= m_dt;
                Eigen::RowVector3d end((i-0.5)*dx + temp, j*dy, 0);
                viewer.data().add_edges(start, end, Eigen::RowVector3d(0,0,0));
            }
        }
        for (unsigned j = 0; j < nx + 1; j++) {
            for (unsigned i = 0; i < ny + 1; i++) {
                Eigen::RowVector3d start(i*dx, (j - 0.5)*dy, 0);
                double temp = p_mac_grid->get_v(i,j);
                temp *= m_dt;
                Eigen::RowVector3d end(i*dx, (j-0.5)*dy + temp, 0);
                viewer.data().add_edges(start, end, Eigen::RowVector3d(0,0,0));
            }
        }
	}
    
    /*Eigen::RowVector3d start(1,1,0);
	Eigen::RowVector3d end(2,2,0);
	viewer.data().add_edges(start, end, Eigen::RowVector3d(1.0,0,0));*/
	//Arrow probe = new Arrow(start, end);
	//addArrow(start, end);*/
}


void WaterSim::initParticles() {

    double sx = p_mac_grid->get_cell_sizex();
    double sy = p_mac_grid->get_cell_sizey();
    std::cout << "Cell size x: " << sx;
    std::cout << "\nCell size y: " << sy << std::endl;
    unsigned nx = p_mac_grid->get_num_cells_x();
    unsigned ny = p_mac_grid->get_num_cells_y();
    unsigned idx = 0;

    // 4 particles per fluid cell
    m_num_particles = 0;// 4 * (nx) * (ny);
    flip_particles = new Particle[4 * (nx) * (ny)];

    Eigen::VectorXd rnd = Eigen::VectorXd::Random(8*nx*ny);
    //~ Eigen::VectorXd rnd = Eigen::VectorXd::Zero(8*nx*ny);

	//~ flip_particles[0] = Particle(sx*(nx/2), sy*(ny/2), 0.);
	
    for (unsigned x = 5; x < 15; x++) {
        for (unsigned y = 5; y < 15; y++) {
            // Populate cell (x,y)
            double cx = x * sx;
            double cy = y * sy;

            flip_particles[idx]     = Particle(cx - sx/4. + sx*rnd(idx  )/8., cy - sy/4. + sy*rnd(idx+1)/8., 0.);
            flip_particles[idx + 1] = Particle(cx + sx/4. + sx*rnd(idx+2)/8., cy - sy/4. + sy*rnd(idx+3)/8., 0.);
            flip_particles[idx + 2] = Particle(cx - sx/4. + sx*rnd(idx+4)/8., cy + sy/4. + sy*rnd(idx+5)/8., 0.);
            flip_particles[idx + 3] = Particle(cx + sx/4. + sx*rnd(idx+6)/8., cy + sy/4. + sy*rnd(idx+7)/8., 0.);

            idx += 4;
            m_num_particles += 4;
        }
    }
}


void WaterSim::initMacGrid() {
    p_mac_grid = new Mac2d(m_res_x, m_res_y, m_len_x, m_len_y);
}

void WaterSim::initFLIP() {
    p_flip = new FLIP(flip_particles, m_num_particles, p_mac_grid,
                      m_fluid_density_, m_gravity_mag_, m_alpha_);
}


void WaterSim::initMacViz() {

    unsigned nx = p_mac_grid->get_num_cells_x();
    unsigned ny = p_mac_grid->get_num_cells_y();
    double sx = p_mac_grid->get_cell_sizex();
    double sy = p_mac_grid->get_cell_sizey();

    unsigned num_vertices = (nx + 1) * (ny + 1);
    unsigned num_edges = nx + ny + 2;
    unsigned num_tria = 2 * nx * ny;

    // Calculate grid vertex coordinates for rendering
    m_renderV.resize(num_vertices, 3);
    unsigned i = 0;
    for (unsigned y = 0; y <= ny; y++) {
        for (unsigned x = 0; x <= nx; x++) {
            // Vertex:
            double cx = x * sx;
            double cy = y * sy;
            m_renderV.row(i) << cx - sx/2., cy - sy/2., 0;
            
            // Increment index
            i++;
        }
    }

    // Calculate edges of MAC grid for rendering
    m_renderE.resize(num_edges, 2);
    m_renderEC.resize(num_edges, 3);
    m_renderEC.setZero();
    i = 0;
    for (unsigned x = 0; x <= nx; x++, i++) {
        m_renderE.row(i) << x, x + ny * (nx + 1);
    }
    for (unsigned y = 0; y <= ny; y++, i++) {
        m_renderE.row(i) << y * (nx + 1), y * (nx + 1) + nx;
    }

    // Calculate mesh triangles (faces)
    //  & colours for rendering
    m_renderF.resize(num_tria, 3);
    m_renderC.resize(num_tria, 3);
    i = 0;
    for (unsigned y = 0; y < ny; y++) {
        for (unsigned x = 0; x < nx; x++) {
            // Triangles
            unsigned vbl = x + y * (nx + 1); // index of bottom-left vertex
            m_renderF.row(i)   << vbl,   vbl+1,              vbl + (nx + 1);
            m_renderF.row(i+1) << vbl+1, vbl + 1 + (nx + 1), vbl + (nx + 1);
            
            // Colour
            //std::cout << x << " " << y << " " << p_mac_grid->is_solid(x,y) << std::endl;
            if (p_mac_grid->is_solid(x, y)) {
                m_renderC.row(i)   << .5, .25, 0;
                m_renderC.row(i+1) << .5, .25, 0;
            } else {
                m_renderC.row(i)   << 1, 1, 0;
                m_renderC.row(i+1) << 1, 1, 0;
            }
            
            i += 2;
        }
    }

}
