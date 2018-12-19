#include "WaterSim.h"


WaterSim::WaterSim(viewer_t& viewer, 
        const int res_x, const int res_y,
        const double len_x, const double len_y,
        const double density, const double gravity,
        const double alpha,
        const bool show_pressures, const bool show_velocity_arrows,
        std::vector<bool> is_fluid, const bool jitter_particles,
        bool export_png, int png_sx, int png_sy, int max_pngs)
        : Simulation(), p_viewer(&viewer), m_res_x(res_x), m_res_y(res_y),
          m_len_x(len_x), m_len_y(len_y), m_fluid_density_(density),
          m_gravity_mag_(gravity), m_alpha_(alpha),
          m_show_pressures(show_pressures),
          m_show_velocity_arrows(show_velocity_arrows),
          is_fluid_(is_fluid), m_jitter_particles(jitter_particles),
          m_export_png_(export_png), m_png_sx_(png_sx), m_png_sy_(png_sy), m_max_pngs_(max_pngs) {


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
    
    // Index of ViewerData instance dedicated to the velocities
    m_velocity_u_idx = viewer.append_mesh();
    m_velocity_v_idx = viewer.append_mesh();
    
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
    
    //Velocity arrows
    p_viewer->data_list[m_velocity_u_idx].clear();
    p_viewer->data_list[m_velocity_v_idx].clear();

    // PNG exporting
    m_png_num_ = 0;
}



void WaterSim::updateParams(const int res_x, const int res_y, 
                  const double len_x, const double len_y,
                  const double density, const double gravity, const double alpha,
                  const bool show_pressures, const bool show_velocity_arrows,
                  std::vector<bool> is_fluid, const bool jitter_particles,
                  bool export_png, int png_sx, int png_sy, int max_pngs) {
    m_res_x = res_x;
    m_res_y = res_y;
    m_len_x = len_x;
    m_len_y = len_y;
    m_fluid_density_ = density;
    m_gravity_mag_ = gravity;
    m_alpha_ = alpha;
    m_show_pressures = show_pressures;
    m_show_velocity_arrows = show_velocity_arrows;
    is_fluid_ = is_fluid;
    m_jitter_particles = jitter_particles;
    m_export_png_ = export_png;
    m_png_sx_ = png_sx;
    m_png_sy_ = png_sy;
    m_max_pngs_ = max_pngs;
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
        for (unsigned j = 0; j < ny; j++) {
            for (unsigned i = 0; i < nx; i++) {
                pressures(2*(i + j*ny)) = p_mac_grid->get_pressure(i, j);
                pressures(2*(i + j*ny) + 1) = p_mac_grid->get_pressure(i, j);
            }
        }

        igl::jet(pressures, true, m_renderC);
    }
    
    if(m_show_velocity_arrows){
		unsigned nx = p_mac_grid->get_num_cells_x();
        unsigned ny = p_mac_grid->get_num_cells_y();
        double dx = p_mac_grid->get_cell_sizex();
        double dy = p_mac_grid->get_cell_sizey();
		m_render_velocity_u_V.resize(2*((nx+1)*ny),3);
		m_render_velocity_u_E.resize((nx+1)*ny,2);
		m_render_velocity_v_V.resize(2*(nx*(ny+1)),3);
		m_render_velocity_v_E.resize(nx*(ny+1),2);
        for (unsigned j = 0; j < ny ; j++) {
            for (unsigned i = 0; i < nx + 1; i++) {
                //insert start vertex
                Eigen::Vector3d start((i-0.5)*dx, j*dy, 0);
                m_render_velocity_u_V.row(2*((nx+1)*j + i)) = start;
                
                //insert end vertex
                double temp = p_mac_grid->get_u(i,j);
                temp *= m_dt;
                Eigen::Vector3d end((i-0.5)*dx + temp, j*dy, 0);
                m_render_velocity_u_V.row(2*((nx+1)*j + i) + 1) = end;
                
                //insert edge
                m_render_velocity_u_E((nx+1)*j + i,0) = 2*((nx+1)*j + i);
                m_render_velocity_u_E((nx+1)*j + i,1) = 2*((nx+1)*j + i) + 1;
            }
        }
        for (unsigned j = 0; j < ny + 1; j++) {
            for (unsigned i = 0; i < nx; i++) {
                //insert start vertex
                Eigen::Vector3d start(i*dx, (j - 0.5)*dy, 0);
                m_render_velocity_v_V.row(2*(nx*j+i)) = start;
                
                //insert end vertex
                double temp = p_mac_grid->get_v(i,j);
                temp *= m_dt;
                Eigen::Vector3d end(i*dx, (j-0.5)*dy + temp, 0);
                m_render_velocity_v_V.row(2*(nx*j+i) + 1) = end;
                
                //insert edge
                m_render_velocity_v_E(nx*j+i,0) = 2*(nx*j+i);
                m_render_velocity_v_E(nx*j+i,1) = 2*(nx*j+i) + 1;
            }
        }
    }

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
        viewer.data_list[m_velocity_u_idx].set_edges(m_render_velocity_u_V, m_render_velocity_u_E, m_renderEC);
        viewer.data_list[m_velocity_v_idx].set_edges(m_render_velocity_v_V, m_render_velocity_v_E, m_renderEC);
    }

    // PNG image export
    if (m_export_png_ and m_png_num_ < m_max_pngs_)
        exportPNG(viewer);
}


void WaterSim::initParticles() {

    double sx = p_mac_grid->get_cell_sizex();
    double sy = p_mac_grid->get_cell_sizey();

    unsigned nx = p_mac_grid->get_num_cells_x();
    unsigned ny = p_mac_grid->get_num_cells_y();

    m_num_particles = 0;
    unsigned idx = 0;
    std::vector<Particle> particles;
    
    Eigen::VectorXd rnd;

    // Random offets to particle positions
    if (m_jitter_particles)
        rnd = Eigen::VectorXd::Random(18*nx*ny);
    else
        rnd = Eigen::VectorXd::Zero(18*nx*ny);

    // Initialize 9 particles per fluid cell
    for (unsigned x = 0; x < nx; x++) {
        for (unsigned y = 0; y < ny; y++) {

            if (!is_fluid_[x + y*nx])
                continue;

            // Populate cell (x,y)
            double cx = x * sx;
            double cy = y * sy;

            particles.push_back(Particle(cx - sx/3. + sx*rnd(idx   )/12., cy - sy/3. + sy*rnd(idx+1 )/12., 0.));
            particles.push_back(Particle(cx         + sx*rnd(idx+2 )/12., cy - sy/3. + sy*rnd(idx+3 )/12., 0.));
            particles.push_back(Particle(cx + sx/3. + sx*rnd(idx+4 )/12., cy - sy/3. + sy*rnd(idx+5 )/12., 0.));

            particles.push_back(Particle(cx - sx/3. + sx*rnd(idx+6 )/12., cy         + sy*rnd(idx+7 )/12., 0.));
            particles.push_back(Particle(cx         + sx*rnd(idx+8 )/12., cy         + sy*rnd(idx+9 )/12., 0.));
            particles.push_back(Particle(cx + sx/3. + sx*rnd(idx+10)/12., cy         + sy*rnd(idx+11)/12., 0.));

            particles.push_back(Particle(cx - sx/3. + sx*rnd(idx+12)/12., cy + sy/3. + sy*rnd(idx+13)/12., 0.));
            particles.push_back(Particle(cx         + sx*rnd(idx+14)/12., cy + sy/3. + sy*rnd(idx+15)/12., 0.));
            particles.push_back(Particle(cx + sx/3. + sx*rnd(idx+16)/12., cy + sy/3. + sy*rnd(idx+17)/12., 0.));

            idx += 9;
            m_num_particles += 9;
        }
    }

    flip_particles = new Particle[m_num_particles];
    std::move(particles.begin(), particles.end(), flip_particles);
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
            
            m_renderC.row(i)   << 1, 1, 0;
            m_renderC.row(i+1) << 1, 1, 0;
            
            i += 2;
        }
    }

}


void WaterSim::exportPNG(igl::opengl::glfw::Viewer &viewer) {
    // Export current frame as PNG
    
    // 1. File name
    const int pad_size = 4; // number of chars to pad number to
    
    std::string file = std::to_string(m_png_num_);
    file.insert(0, pad_size - file.length(), '0');
    file += ".png";
    
    // 2. Check for folder existence, create if necessary
    // http://pubs.opengroup.org/onlinepubs/009695399/functions/mkdir.html
    mkdir(m_png_dirname_.c_str(), S_IRWXU);
    
    // 3. Draw buffer
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(m_png_sx_,m_png_sy_);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(m_png_sx_,m_png_sy_);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(m_png_sx_,m_png_sy_);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(m_png_sx_,m_png_sy_);
    
    viewer.core.draw_buffer(viewer.data_list[m_particles_data_idx], false, R,G,B,A);

    // 4. Write image
    igl::png::writePNG(R,G,B,A, m_png_dirname_ + '/' + file);
    
    m_png_num_++;
    
}

