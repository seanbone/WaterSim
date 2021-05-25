#include "tsc_x86.hpp"
#include "WaterSim.h"

WaterSim::WaterSim(const SimConfig& cfg) : m_cfg(cfg) {

	// Initialize MAC grid
	initMacGrid();

	// Initialize particles
	initParticles();

	// Initialize FLIP object
	initFLIP();

	// Initialize Mesh Exporter
	initMeshExp();
}

void WaterSim::updateParams(const SimConfig& cfg) {
	m_cfg = cfg;
	setTimestep(cfg.getTimeStep());
}

void WaterSim::resetMembers() {
	m_time = 0;
	m_step = 0;

	// MAC grid
	delete p_mac_grid;
	initMacGrid();

	// Particles
	delete flip_particles;
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

	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");
	tsctimer.start_timing("FLIP");
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
	tsctimer.stop_timing("FLIP", true, "");
    auto flip_duration = duration_cast<ticks_t>( t2 - t1 ).count() / timer_scale;
    std::cout << "\nFLIP duration: " << flip_duration << timer_unit <<  std::endl;


    // Compute mesh if required
    if (m_cfg.getDisplayMeshes() || m_cfg.getExportMeshes()) {
		tsctimer.start_timing("compute_mesh");
	    exp->compute_mesh();
		tsctimer.stop_timing("compute_mesh", true, "");
    }

	// Export mesh if required
    if (m_cfg.getExportMeshes()) {
		tsctimer.start_timing("export_mesh");
        exp->export_mesh();
		tsctimer.stop_timing("export_mesh", true, "");
    }


    // Time at end of simulation step
    tpoint_t tf = timer_t::now();
    auto tot_duration = duration_cast<ticks_t>( tf - t1 ).count() / timer_scale;
    std::cout << "\nTotal duration: " << tot_duration << timer_unit <<  std::endl;

    // Increase counter and current time
    m_step++;
    m_time += m_dt;
	tsctimer.step();

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
    std::list<double> particles_x;
	std::list<double> particles_y;
	std::list<double> particles_z;

	// Random offsets to particle positions
	unsigned int n = 3 * particles_per_cell * nx * ny * nz;
	Eigen::VectorXd rnd = Eigen::VectorXd::Zero(n);
    if (m_cfg.getJitterParticles()) {
    	int seed = m_cfg.getRandomSeed();
    	std::srand((unsigned int) ((seed >= 0) ? seed : std::time(nullptr)));
	    rnd = Eigen::VectorXd::Random(n);
    }

    // Get fluid region
	double fluid_from_x, fluid_from_y, fluid_from_z;
	double fluid_to_x, fluid_to_y, fluid_to_z;
	m_cfg.getFluidRegion(fluid_from_x, fluid_from_y, fluid_from_z,
	                     fluid_to_x, fluid_to_y, fluid_to_z);

	auto is_fluid = [&](double cx, double cy, double cz) {
		return cx >= fluid_from_x && cx <= fluid_to_x
			&& cy >= fluid_from_y && cy <= fluid_to_y
			&& cz >= fluid_from_z && cz <= fluid_to_z;
	};

	// Initialize particles_per_cell particles per fluid cell
    for (unsigned z = 0; z < nz; z++) {
        for (unsigned y = 0; y < ny; y++) {
            for (unsigned x = 0; x < nx; x++) {

	            // Center of cell (x,y,z)
	            double cx = x * sx;
	            double cy = y * sy;
	            double cz = z * sz;

	            // Only populate cells flagged as fluid
	            if (!is_fluid(cx, cy, cz))
                    continue;

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

                for (auto & position : positions) {
                    particles_x.push_back(position[0]);
	                particles_y.push_back(position[1]);
	                particles_z.push_back(position[2]);
                }

                idx += particles_per_cell;
                m_num_particles += particles_per_cell;
            }
        }
    }

	flip_particles = new Particles(m_num_particles, *p_mac_grid);
	std::move(particles_x.begin(), particles_x.end(), flip_particles->x);
	std::move(particles_y.begin(), particles_y.end(), flip_particles->y);
	std::move(particles_z.begin(), particles_z.end(), flip_particles->z);
}

void WaterSim::initMacGrid() {
	int res_x, res_y, res_z;
	m_cfg.getGridResolution(res_x, res_y, res_z);
	double len_x, len_y, len_z;
	m_cfg.getSystemSize(len_x, len_y, len_z);
    p_mac_grid = new Mac3d(res_x, res_y, res_z, len_x, len_y, len_z);
}

void WaterSim::initMeshExp() {
	exp = new MeshExporter(p_mac_grid, *flip_particles);
}

void WaterSim::initFLIP() {
    p_flip = new FLIP(*flip_particles, p_mac_grid, m_cfg);
}
