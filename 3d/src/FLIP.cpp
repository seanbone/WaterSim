#include "FLIP.h"
#include "tsc_x86.hpp"


FLIP::FLIP(Particle* particlesOLD, Particles& particles, unsigned num_particles, Mac3d* MACGrid,
		   const SimConfig& cfg)
	: particlesOLD_(particlesOLD), num_particles_(num_particles), MACGrid_(MACGrid), cfg_(cfg),
      fluid_density_(cfg.getDensity()), gravity_mag_(cfg.getGravity()), alpha_(cfg.getAlpha()),
      particles_(particles) {
	
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
	unsigned n, m, l;
	double* x;
	double* y;
	double* z;
	double* u;
	double* v;
	double* w;
	double* uMAC;
	double* vMAC;
	double* wMAC;
	double* uStar;
	double* vStar;
	double* wStar;
	double* pMAC;
	bool* fluid_cells;
	bool* solid_cells;

	if (step == WRITE_REFERENCE){

		unsigned cacheBlockSize = 64;

		n = MACGrid_->get_num_cells_x();
		m = MACGrid_->get_num_cells_y();
		l = MACGrid_->get_num_cells_z();

		x     = (double*) aligned_alloc(cacheBlockSize, num_particles_ * sizeof(double));
		y     = (double*) aligned_alloc(cacheBlockSize, num_particles_ * sizeof(double));
		z     = (double*) aligned_alloc(cacheBlockSize, num_particles_ * sizeof(double));
		u     = (double*) aligned_alloc(cacheBlockSize, num_particles_ * sizeof(double));
		v     = (double*) aligned_alloc(cacheBlockSize, num_particles_ * sizeof(double));
		w     = (double*) aligned_alloc(cacheBlockSize, num_particles_ * sizeof(double));
		uMAC  = (double*) aligned_alloc(cacheBlockSize, (n+1) * m * l  * sizeof(double));
		vMAC  = (double*) aligned_alloc(cacheBlockSize, n * (m+1) * l  * sizeof(double));
		wMAC  = (double*) aligned_alloc(cacheBlockSize, n * m * (l+1)  * sizeof(double));
		uStar = (double*) aligned_alloc(cacheBlockSize, (n+1) * m * l  * sizeof(double));
		vStar = (double*) aligned_alloc(cacheBlockSize, n * (m+1) * l  * sizeof(double));
		wStar = (double*) aligned_alloc(cacheBlockSize, n * m * (l+1)  * sizeof(double));
		pMAC  = (double*) aligned_alloc(cacheBlockSize, n * m * l      * sizeof(double));

		fluid_cells = (bool*) aligned_alloc(cacheBlockSize, n * m * l * sizeof(bool));
		solid_cells = (bool*) aligned_alloc(cacheBlockSize, n * m * l * sizeof(bool));

		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(0, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
	}
#endif

	// 1.
	tsctimer.start_timing("particle_to_grid");
    particle_to_grid();

	// 1a.
	MACGrid_->set_uvw_star();
	tsctimer.stop_timing("particle_to_grid", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE){
		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(1, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
	}
#endif

	// 2.
	tsctimer.start_timing("apply_forces");
	apply_forces(dt);
	
	if (cfg_.getApplyMeteorForce() && step <= 200 ){
		explode(dt, step, 15, 0, 15, 2, 800);
	}

	tsctimer.stop_timing("apply_forces", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE){
		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(2, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
	}
#endif

	// 3.
	tsctimer.start_timing("apply_boundary_conditions");
	apply_boundary_conditions();
	tsctimer.stop_timing("apply_boundary_conditions", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE){
		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(3, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
	}
#endif

	// 4.
	tsctimer.start_timing("apply_pressure_correction");
    apply_pressure_correction(dt);
	tsctimer.stop_timing("apply_pressure_correction", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE){
		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(4, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
	}
#endif

	// 5.
	tsctimer.start_timing("grid_to_particle");
	grid_to_particle();
	tsctimer.stop_timing("grid_to_particle", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE){
		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(5, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
	}
#endif

	// 6.
	tsctimer.start_timing("advance_particles");
	double dt_new = compute_timestep(dt);
	double num_substeps = std::ceil(dt/dt_new);
	for( int s = 0; s < num_substeps ; ++s ){
		
		// 7.
		advance_particles(dt/num_substeps, step);
	}
	tsctimer.stop_timing("advance_particles", true, "");

#ifdef WRITE_REFERENCE
	if (step == WRITE_REFERENCE){
		ncWriter_->toLinArrays(particlesOLD_, num_particles_, MACGrid_, n, m, l, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);
		ncWriter_->writeAll(6, x, y, z, u, v, w, uMAC, vMAC, wMAC, uStar, vStar, wStar, pMAC, fluid_cells, solid_cells);

		free(x);
		free(y);
		free(z);
		free(u);
		free(v);
		free(w);
		free(uMAC);
		free(vMAC);
		free(wMAC);
		free(uStar);
		free(vStar);
		free(wStar);
		free(pMAC);
		free(fluid_cells);
		free(solid_cells);
	}
#endif
}


void FLIP::particlesOldToNew() {
	Particles::particleIdx_t i;
	for (i = 0; i < particles_.get_num_particles(); i++) {
		particles_.x[i] = particlesOLD_[i].x_;
		particles_.y[i] = particlesOLD_[i].y_;
		particles_.z[i] = particlesOLD_[i].z_;
		particles_.u[i] = particlesOLD_[i].u_;
		particles_.v[i] = particlesOLD_[i].v_;
		particles_.w[i] = particlesOLD_[i].w_;
	}
}

void FLIP::particlesNewToOld() {
	Particles::particleIdx_t i;
	for (i = 0; i < particles_.get_num_particles(); i++) {
		particlesOLD_[i].x_ = particles_.x[i];
		particlesOLD_[i].y_ = particles_.y[i];
		particlesOLD_[i].z_ = particles_.z[i];
		particlesOLD_[i].u_ = particles_.u[i];
		particlesOLD_[i].v_ = particles_.v[i];
		particlesOLD_[i].w_ = particles_.w[i];
	}
}