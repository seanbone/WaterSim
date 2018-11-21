#include "WaterSim.h"


// WaterSim(viewer_t& viewer, const int res_x, const int res_y, const double len_x, const double len_y);
WaterSim::WaterSim(WaterSim::viewer_t& viewer,
        const int res_x, const int res_y,
        const double len_x, const double len_y)
        : Simulation(), p_mac_grid(new Mac2d(res_x, res_y, len_x, len_y)) {

    init(viewer);
}


void WaterSim::init(WaterSim::viewer_t& viewer) {
    
    // Initialize particles
    // Fill top left of box with flip particles
    double sx = p_mac_grid->get_cell_sizex();
    double sy = p_mac_grid->get_cell_sizey();
    std::cout << "Cell size x: " << sx;
    std::cout << "\nCell size y: " << sy << std::endl;
    unsigned nx = p_mac_grid->get_num_cells_x();
    unsigned ny = p_mac_grid->get_num_cells_y();
    unsigned idx = 0;

    // 4 particles per fluid cell
    m_num_particles = 4 * (nx/2 - 1) * (ny/2 - 1);
    flip_particles = new Particle[m_num_particles];

    Eigen::VectorXd rnd = Eigen::VectorXd::Random(2*m_num_particles);

    for (unsigned x = 3; x < (nx/2)+2; x++) {
        for (unsigned y = 1; y < ny/2; y++) {
            // Populate cell (x,y)
            double cx = x * sx;
            double cy = y * sy;

            flip_particles[idx]     = Particle(cx - sx/4. + sx*rnd(idx  )/8., cy - sy/4. + sy*rnd(idx+1)/8., 0.);
            flip_particles[idx + 1] = Particle(cx + sx/4. + sx*rnd(idx+2)/8., cy - sy/4. + sy*rnd(idx+3)/8., 0.);
            flip_particles[idx + 2] = Particle(cx - sx/4. + sx*rnd(idx+4)/8., cy + sy/4. + sy*rnd(idx+5)/8., 0.);
            flip_particles[idx + 3] = Particle(cx + sx/4. + sx*rnd(idx+6)/8., cy + sy/4. + sy*rnd(idx+7)/8., 0.);

            idx += 4;
        }
    }

    // Initialize FLIP object
    p_flip = new FLIP(flip_particles, m_num_particles, p_mac_grid);
    
    
    // Index of ViewerData instance dedicated to particles
    m_particles_data_idx = viewer.append_mesh();

    
    // Generate visualization mesh based on Mac2d

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

    m_grid_data_idx = viewer.append_mesh();
    
    // Update rendering geometry
    updateRenderGeometry();
    
    // Initial reset
    reset();
}


void WaterSim::resetMembers() {
    m_particle_colors.setZero();
    m_particle_colors.col(2).setOnes();
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

    //igl::jet(pressures, true, m_renderC);

    std::cout << "\n*************\n";
    std::cout << "Pressure at (7, 0): " << p_mac_grid->get_pressure(7, 0) << std::endl;
    std::cout << "Pressure at (7, 1): " << p_mac_grid->get_pressure(7, 1) << std::endl;
    std::cout << "U velocity at (6.5, 1): " << p_mac_grid->get_u(7, 1) << std::endl;
    std::cout << "V velocity at (6.5, 1): " << p_mac_grid->get_v(7, 1) << std::endl;
    std::cout << "U* velocity at (6.5, 1): " << p_mac_grid->get_u_star(7, 1) << std::endl;
    std::cout << "V* velocity at (6.5, 1): " << p_mac_grid->get_v_star(7, 1) << std::endl;
    std::cout << "X of particle 0: " << flip_particles->get_position()(0) << std::endl;
    std::cout << "Y of particle 0: " << flip_particles->get_position()(1) << std::endl;
    std::cout << "U of particle 0: " << flip_particles->get_velocity()(0) << std::endl;
    std::cout << "V of particle 0: " << flip_particles->get_velocity()(1) << std::endl;
    std::cout << "\n*************\n";
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
}

