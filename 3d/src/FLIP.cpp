#include "FLIP.h"
#include "tsc_x86.hpp"
//#define WRITE_REFERENCE 340


FLIP::FLIP(Particles& particles, Mac3d* MACGrid, const SimConfig& cfg)
	: cfg_(cfg), particles_(particles), num_particles_(particles.get_num_particles()),
	  fluid_density_(cfg.getDensity()), gravity_mag_(cfg.getGravity()), alpha_(cfg.getAlpha()),
	  MACGrid_(MACGrid), cg_solver(100, *MACGrid_) {
	
    unsigned nx = MACGrid_->get_num_cells_x();
    unsigned ny = MACGrid_->get_num_cells_y();
    unsigned nz = MACGrid_->get_num_cells_z();
    d_ = new (std::align_val_t(32)) double [nx*ny*nz];

#ifdef WRITE_REFERENCE
	ncWriter_ = new NcWriter( "./ref.nc", 
							  7, 
							  num_particles_, 
							  MACGrid_->get_num_cells_x(), 
							  MACGrid_->get_num_cells_y(), 
							  MACGrid_->get_num_cells_z() );

	unsigned timestep = WRITE_REFERENCE;

	ncWriter_->writeScalar("timestep", &timestep, ncUint);

	ncWriter_->addVar("x", "num_particles", ncDouble);
	ncWriter_->addVar("y", "num_particles", ncDouble);
	ncWriter_->addVar("z", "num_particles", ncDouble);
	ncWriter_->addVar("u", "num_particles", ncDouble);
	ncWriter_->addVar("v", "num_particles", ncDouble);
	ncWriter_->addVar("w", "num_particles", ncDouble);

	ncWriter_->addVar("uMAC",  "x_mac_faces", ncDouble);
	ncWriter_->addVar("vMAC",  "y_mac_faces", ncDouble);
	ncWriter_->addVar("wMAC",  "z_mac_faces", ncDouble);
	ncWriter_->addVar("uStar", "x_mac_faces", ncDouble);
	ncWriter_->addVar("vStar", "y_mac_faces", ncDouble);
	ncWriter_->addVar("wStar", "z_mac_faces", ncDouble);
	
	ncWriter_->addVar("pMAC", "mac_centers", ncDouble);
	ncWriter_->addVar("fluid_cells", "mac_centers", ncByte);
	ncWriter_->addVar("solid_cells", "mac_centers", ncByte);
#endif
}


FLIP::~FLIP(){
    delete [] d_;

#ifdef WRITE_REFERENCE
	delete ncWriter_;
#endif
}


/*** PERFORM ONE STEP ***/
void FLIP::step_FLIP(double dt, unsigned long step) {
	/** One FLIP step:
	 * 1. Compute velocity field (particle-to-grid transfer)
	 *    - Particle-to-grid transfer
	 *    - Classify cells (fluid/air)
	 *    - Extrapolate velocity field into air region
	 * 1a. Copy velocity field to intermediate velocity field u^*
	 * 2. Apply external forces (fwd euler on field)
	 * 3. Enforce boundary conditions for grid & solid boundaries
	 * 4. Compute & apply pressure gradients
	 * 5. Update particle velocities from grid-velocities
	 * 6. Subsample time interval to satisfy CFL condition
	 * 7. Update particle positions
	 */

	tsc::TSCTimer& tsctimer = tsc::TSCTimer::get_timer("timings.json");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(0, particles_, MACGrid_);
#endif

	// 1.
	tsctimer.start_timing("particle_to_grid");
    particle_to_grid();

	// 1a.
	MACGrid_->set_uvw_star();
	tsctimer.stop_timing("particle_to_grid", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(1, particles_, MACGrid_);
#endif

	// 2.
	tsctimer.start_timing("apply_forces");
	apply_forces(dt);
	
	if (cfg_.getApplyMeteorForce() && step <= 200 ){
		explode(dt, step, 15, 0, 15, 2., 800.);
	}

	tsctimer.stop_timing("apply_forces", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(2, particles_, MACGrid_);
#endif

	// 3.
	tsctimer.start_timing("apply_boundary_conditions");
	apply_boundary_conditions();
	tsctimer.stop_timing("apply_boundary_conditions", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(3, particles_, MACGrid_);
#endif

	// 4.
	tsctimer.start_timing("apply_pressure_correction");
    apply_pressure_correction(dt);
	tsctimer.stop_timing("apply_pressure_correction", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(4, particles_, MACGrid_);
#endif

	// 5.
	tsctimer.start_timing("grid_to_particle");
	grid_to_particle();
	tsctimer.stop_timing("grid_to_particle", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(5, particles_, MACGrid_);
#endif

	// 6.
	tsctimer.start_timing("advance_particles");
	double dt_new = compute_timestep(dt);
	double num_substeps = std::ceil(dt/dt_new);
	for( int s = 0; s < num_substeps ; ++s ){
		
		// 7.
		advance_particles(dt/num_substeps);
	}
	tsctimer.stop_timing("advance_particles", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE) ncWriter_->writeAll(6, particles_, MACGrid_);
#endif
}
