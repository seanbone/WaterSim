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
	
	FLIP* flip = new FLIP(ncReader->particles, ncReader->num_particles, ncReader->MACGrid, ncReader->cfg);
	
	flip->compute_velocity_field();

	ncReader->MACGrid->set_uvw_star();

	ncReader->readAll(1);
	ncReader->validate();

	delete ncReader;
	delete flip;
}