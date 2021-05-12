/*
 * A test to check validation of the particle-to-grid step
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(0);
	ncReader->toOldStruct();

	// TODO: initialize Particles struct in ncReader instead.
	Particles particles(ncReader->num_particles, ncReader->cfg, *ncReader->MACGrid);
	FLIP* flip = new FLIP(ncReader->particlesOLD, particles, ncReader->num_particles, ncReader->MACGrid, ncReader->cfg);
	// TODO: remove once new struct is integrated
    flip->particlesOldToNew();

	flip->particle_to_grid();

	flip->particlesNewToOld();

	ncReader->MACGrid->set_uvw_star();

	ncReader->readAll(1);
	ncReader->validate();

	delete ncReader;
	delete flip;
}