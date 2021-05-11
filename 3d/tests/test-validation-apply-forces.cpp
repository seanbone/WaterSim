/*
 * A test to check validation of the step where external forces are applied
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(1);
	ncReader->toOldStruct();

	// TODO: initialize Particles struct in ncReader instead.
	Particles particles(ncReader->num_particles, ncReader->cfg, *ncReader->MACGrid);
	FLIP* flip = new FLIP(ncReader->particlesOLD, particles, ncReader->num_particles, ncReader->MACGrid, ncReader->cfg);
	// TODO: remove once new struct is integrated
	flip->particlesOldToNew();

	double dt = ncReader->cfg.getTimeStep();

	flip->apply_forces(dt);
	
	unsigned step = ncReader->timestep;
	if (ncReader->cfg.getApplyMeteorForce() && step <= 200 ){
		flip->explode(dt, step, 15, 0, 15, 2, 800);
	}

	ncReader->readAll(2);
	ncReader->validate();

	delete ncReader;
	delete flip;
}