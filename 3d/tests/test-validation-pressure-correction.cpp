/*
 * A test to check validation of the pressure-correction step
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(3);
	ncReader->toOldStruct();

	FLIP* flip = new FLIP(ncReader->particles, ncReader->num_particles, ncReader->MACGrid, ncReader->cfg);

	double dt = ncReader->cfg.getTimeStep();

	flip->do_pressures(dt);

	ncReader->readAll(4);
	ncReader->validate();

	delete ncReader;
	delete flip;
}