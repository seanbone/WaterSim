/*
 * A test to check validation of the pressure-correction step
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(4);
	ncReader->toFlipStructures();

	FLIP* flip = new FLIP(*(ncReader->particles), ncReader->MACGrid, ncReader->cfg);

	flip->grid_to_particle();

	ncReader->readAll(5);
	ncReader->validate();

	delete ncReader;
	delete flip;
}