#include "WaterSim.h"

WaterSim::WaterSim(const SimConfig& cfg, std::vector<bool> is_fluid)
	: m_cfg(cfg), is_fluid_(std::move(is_fluid)) {

	// Initialize MAC grid
	initMacGrid();

	// Initialize particles
	initParticles();

	// Initialize FLIP object
	initFLIP();

	// Initialize Mesh Exporter
	initMeshExp();
}

void WaterSim::updateParams(const SimConfig& cfg, std::vector<bool> is_fluid) {
	m_cfg = cfg;
	is_fluid_ = std::move(is_fluid);
	setTimestep(cfg.getTimeStep());
}

void WaterSim::resetMembers() {
	// MAC grid
	delete p_mac_grid;
	initMacGrid();

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

bool WaterSim::advance() {

	if (m_step == 0)
        std::cout << "Starting simulation with " << m_num_particles << " particles.\n\n";

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
    if (m_cfg.getExportMeshes()) {
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
    
    // Random offsets to particle positions
    Eigen::VectorXd rnd;
    if (m_cfg.getJitterParticles())
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
	int res_x, res_y, res_z;
	m_cfg.getGridResolution(res_x, res_y, res_z);
	double len_x, len_y, len_z;
	m_cfg.getSystemSize(len_x, len_y, len_z);
    p_mac_grid = new Mac3d(res_x, res_y, res_z, len_x, len_y, len_z);
}

void WaterSim::initMeshExp(){
	exp = new MeshExporter(p_mac_grid, flip_particles, m_num_particles);
}

void WaterSim::initFLIP() {
    p_flip = new FLIP(flip_particles, m_num_particles, p_mac_grid,
                      m_cfg.getDensity(), m_cfg.getGravity(),
                      m_cfg.getAlpha());
}