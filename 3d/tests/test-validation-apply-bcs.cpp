/*
 * A test to check validation of the step where boundary conditions are enforced
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(2);
	ncReader->toFlipStructures();

	FLIP* flip = new FLIP(*(ncReader->particles), ncReader->MACGrid, ncReader->cfg);

	flip->apply_boundary_conditions();

	ncReader->readAll(3);
	ncReader->validate();

	delete ncReader;
	delete flip;
}