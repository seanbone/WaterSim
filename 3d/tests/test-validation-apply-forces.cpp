/*
 * A test to check validation of the step where external forces are applied
 */
#include "includes/watersim-test-common.h"
#include "NcReader.h"
#include "FLIP.h"


int main(){

	NcReader* ncReader = new NcReader(validation_data_ref, validation_data_cfg);

	ncReader->readAll(1);
	ncReader->toFlipStructures();

	FLIP* flip = new FLIP(*(ncReader->particles), ncReader->MACGrid, ncReader->cfg);

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