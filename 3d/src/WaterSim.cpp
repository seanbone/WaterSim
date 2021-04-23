#include "WaterSim.h"


WaterSim::WaterSim(viewer_t& viewer, const bool display_grid,
        const int res_x, const int res_y, const int res_z,
        const double len_x, const double len_y, const double len_z,
        const double density, const double gravity,
        const double alpha,
        std::vector<bool> is_fluid, const bool jitter_particles,
        bool export_meshes, unsigned max_p)
        : Simulation(), p_viewer(&viewer), m_display_grid(display_grid),
          m_res_x(res_x), m_res_y(res_y), m_res_z(res_z),
          m_len_x(len_x), m_len_y(len_y), m_len_z(len_z), m_fluid_density_(density),
          m_gravity_mag_(gravity), m_alpha_(alpha),
          is_fluid_(is_fluid), m_jitter_particles(jitter_particles),
          m_export_meshes(export_meshes), m_max_p_disp(max_p) {


    // Initialize MAC grid
    initMacGrid();

    // Initialize particles
    initParticles();

    // Initialize FLIP object
    initFLIP();
    
    // Index of ViewerData instance dedicated to particles
    m_particles_data_idx = viewer.append_mesh();
    
    // Generate visualization mesh based on Mac3d
    initMacViz();

    // Index of ViewerData instance dedicated to MAC grid
    m_grid_data_idx = viewer.append_mesh();
    
    // Initialize Mesh Exporter
    initMeshExp();
    
    // Update rendering geometry
    updateRenderGeometry();
}


/**
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
    
    // MeshExporter
    delete exp;
    initMeshExp();
}


void WaterSim::updateParams(const bool display_grid,
                  const int res_x, const int res_y, const int res_z,
                  const double len_x, const double len_y, const double len_z,
                  const double density, const double gravity, const double alpha,
                  std::vector<bool> is_fluid, const bool jitter_particles,
                  bool export_meshes, unsigned max_p) {
    m_display_grid = display_grid;
    m_res_x = res_x;
    m_res_y = res_y;
    m_res_z = res_z;
    m_len_x = len_x;
    m_len_y = len_y;
    m_len_z = len_z;
    m_fluid_density_ = density;
    m_gravity_mag_ = gravity;
    m_alpha_ = alpha;
    is_fluid_ = is_fluid;
    m_jitter_particles = jitter_particles;
    m_export_meshes = export_meshes;
    m_max_p_disp = max_p;
    std::cout << "\nParams updated\n";
}


void WaterSim::updateRenderGeometry() {
    // Copy particle positions from FLIP's data structure
    unsigned disp_particles = m_num_particles;
    unsigned particle_step = 1;

    // If there are too many particles to display, only
    // display a subset, using stride particle_set
    if (m_num_particles > m_max_p_disp) {
        disp_particles = m_max_p_disp;
        particle_step = m_num_particles / m_max_p_disp;
    }
     
    m_particles.resize(disp_particles, 3);
    for (unsigned i = 0, j=0; j < disp_particles && i < m_num_particles; j++, i += particle_step) {
        m_particles.row(j) = flip_particles[i].get_position();
    }

    m_particle_colors.resize(disp_particles, 3);
    m_particle_colors.setZero();
    m_particle_colors.col(2).setOnes();
}


bool WaterSim::advance() {

    if (m_step == 0)
        std::cout << "Starting simluation with " << m_num_particles << " particles.\n\n";

    std::cout << "\n\nBegin FLIP step #" << m_step << std::endl;

    // Set up for timings of individual functions
    using timer_t = std::chrono::high_resolution_clock;
    using tpoint_t = timer_t::time_point;
    using namespace std::chrono;
    using ticks_t = std::chrono::microseconds;
    double timer_scale = 1e6;
    auto timer_unit = "s";
    
    // Time at beginning of simulation step
    tpoint_t t1 = timer_t::now();

    // Perform a FLIP step
    p_flip->step_FLIP(m_dt, m_step);
    
    tpoint_t t2 = timer_t::now();
    auto flip_duration = duration_cast<ticks_t>( t2 - t1 ).count() / timer_scale;
    std::cout << "\nFLIP duration: " << flip_duration << timer_unit <<  std::endl;


    // Export mesh if required
    if (m_export_meshes) {
        std::cout << "\nExport mesh..." << std::endl;
        exp->export_mesh();
        tpoint_t t3 = timer_t::now();
        auto export_duration = duration_cast<ticks_t>( t3 - t2 ).count() / timer_scale;
        std::cout << "\nExport duration: " << export_duration << timer_unit <<  std::endl;
    }


    // Time at end of simulation step
    tpoint_t tf = timer_t::now();
    auto tot_duration = duration_cast<ticks_t>( tf - t1 ).count() / timer_scale;
    std::cout << "\nTotal duration: " << tot_duration << timer_unit <<  std::endl;

    // Increase counter and current time
    m_step++;
    m_time += m_dt;

    return false;
}


void WaterSim::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer) {
    // Display MAC grid
    if (m_display_grid){
        viewer.data_list[m_grid_data_idx].set_edges(m_renderV, m_renderE, m_renderEC);
    }

    // Display particles
    viewer.data_list[m_particles_data_idx].set_points(m_particles, m_particle_colors);
    viewer.data_list[m_particles_data_idx].point_size = 5;
}


void WaterSim::initParticles() {

    double sx = p_mac_grid->get_cell_sizex();
    double sy = p_mac_grid->get_cell_sizey();
    double sz = p_mac_grid->get_cell_sizez();

    unsigned nx = p_mac_grid->get_num_cells_x();
    unsigned ny = p_mac_grid->get_num_cells_y();
    unsigned nz = p_mac_grid->get_num_cells_z();

    // Note: this must equal the actual number of particles
    // initialized in the loop below
    const unsigned particles_per_cell = 8;

    m_num_particles = 0;
    unsigned idx = 0;

    // Constant complexity for .push_back
    std::list<Particle> particles;
    
    // Random offets to particle positions
    Eigen::VectorXd rnd;
    if (m_jitter_particles)
        rnd = Eigen::VectorXd::Random(3*particles_per_cell*nx*ny*nz);
    else
        rnd = Eigen::VectorXd::Zero(3*particles_per_cell*nx*ny*nz);

    // Initialize particles_per_cell particles per fluid cell
    for (unsigned z = 0; z < nz; z++) {
        for (unsigned y = 0; y < ny; y++) {
            for (unsigned x = 0; x < nx; x++) {

                // Only populate cells flagged as fluid
                if (!is_fluid_[x + y*nx + z*nx*ny])
                    continue;

                // Center of cell (x,y,z)
                double cx = x * sx;
                double cy = y * sy;
                double cz = z * sz;

                const double f = 4.;

                // Particle positions in a 2x2x2 array around cell center
                double positions[particles_per_cell][3] = {
                    { cx - sx/f + rnd(idx   )*sx/(2*f), cy - sy/f + rnd(idx+1 )*sy/(2*f), cz - sz/f + rnd(idx+2 )*sz/(2*f)},
                    { cx - sx/f + rnd(idx+3 )*sx/(2*f), cy - sy/f + rnd(idx+4 )*sy/(2*f), cz + sz/f + rnd(idx+5 )*sz/(2*f)},
                    { cx - sx/f + rnd(idx+6 )*sx/(2*f), cy + sy/f + rnd(idx+7 )*sy/(2*f), cz - sz/f + rnd(idx+8 )*sz/(2*f)},
                    { cx - sx/f + rnd(idx+9 )*sx/(2*f), cy + sy/f + rnd(idx+10)*sy/(2*f), cz + sz/f + rnd(idx+11)*sz/(2*f)},
                    { cx + sx/f + rnd(idx+12)*sx/(2*f), cy - sy/f + rnd(idx+13)*sy/(2*f), cz - sz/f + rnd(idx+14)*sz/(2*f)},
                    { cx + sx/f + rnd(idx+15)*sx/(2*f), cy - sy/f + rnd(idx+16)*sy/(2*f), cz + sz/f + rnd(idx+17)*sz/(2*f)},
                    { cx + sx/f + rnd(idx+18)*sx/(2*f), cy + sy/f + rnd(idx+19)*sy/(2*f), cz - sz/f + rnd(idx+20)*sz/(2*f)},
                    { cx + sx/f + rnd(idx+21)*sx/(2*f), cy + sy/f + rnd(idx+22)*sy/(2*f), cz + sz/f + rnd(idx+23)*sz/(2*f)}
                };

                for (unsigned p = 0; p < particles_per_cell; p++) {
                    particles.push_back(Particle(positions[p][0], positions[p][1], positions[p][2]));
                }

                idx += particles_per_cell;
                m_num_particles += particles_per_cell;
            }
        }
    }

    flip_particles = new Particle[m_num_particles];
    std::move(particles.begin(), particles.end(), flip_particles);
}

void WaterSim::initMacGrid() {
    p_mac_grid = new Mac3d(m_res_x, m_res_y, m_res_z, m_len_x, m_len_y, m_len_z);
}

void WaterSim::initMeshExp(){
	exp = new MeshExporter(p_mac_grid, flip_particles, m_num_particles);
}

void WaterSim::initFLIP() {
    p_flip = new FLIP(flip_particles, m_num_particles, p_mac_grid,
                      m_fluid_density_, m_gravity_mag_, m_alpha_);
}


void WaterSim::initMacViz() {

    if (!m_display_grid)
        return;

    unsigned nx = p_mac_grid->get_num_cells_x();
    unsigned ny = p_mac_grid->get_num_cells_y();
    unsigned nz = p_mac_grid->get_num_cells_z();
    double sx = p_mac_grid->get_cell_sizex();
    double sy = p_mac_grid->get_cell_sizey();
    double sz = p_mac_grid->get_cell_sizez();

    unsigned num_vertices = (nx + 1) * (ny + 1) * (nz + 1);
    unsigned num_edges = (nx + ny + 2) * (nz + 1) + (nx + 1)*(ny + 1);

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
                m_renderV.row(i) << cx - sx/2., cy - sy/2., cz - sz/2.;
                
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
            auto a = y*(nx+1) + z*(nx+1)*(ny+1);
            m_renderE.row(i) << a, a + nx;
        }
    }
    // Edges parallel to y axis <=> normal to xz plane
    for (unsigned z = 0; z <= nz; z++) {
        for (unsigned x = 0; x <= nx; x++, i++) {
            auto a = x + z*(nx+1)*(ny+1);
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

