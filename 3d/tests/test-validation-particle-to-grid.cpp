/*
 * A test to check validation of the particle-to-grid step
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(0);
	ncReader->toFlipStructures();

	FLIP* flip = new FLIP(*(ncReader->particles), ncReader->MACGrid, ncReader->cfg);

	flip->particle_to_grid();

	ncReader->MACGrid->set_uvw_star();

	ncReader->readAll(1);
	ncReader->validate();

	delete ncReader;
	delete flip;
}