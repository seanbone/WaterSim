/*
 * A test to check validation of the pressure-correction step
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(5);
	ncReader->toOldStruct();

	// TODO: initialize Particles struct in ncReader instead.
	Particles particles(ncReader->num_particles, ncReader->cfg, *ncReader->MACGrid);
	FLIP* flip = new FLIP(ncReader->particlesOLD, particles, ncReader->num_particles, ncReader->MACGrid, ncReader->cfg);
	// TODO: remove once new struct is integrated
	flip->particlesOldToNew();

	unsigned step = ncReader->timestep;
	double dt = ncReader->cfg.getTimeStep();
	double dt_new = flip->compute_timestep(dt);
	double num_substeps = std::ceil(dt/dt_new);
	for( int s = 0; s < num_substeps ; ++s ){
		
		flip->advance_particles(dt/num_substeps, step);
	}

	ncReader->readAll(6);
	ncReader->validate();

	delete ncReader;
	delete flip;
}