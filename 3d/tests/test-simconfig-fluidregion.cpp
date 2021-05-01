/*
 * A test to check the reordering property of SimConfig::setFluidRegion
 */
#include "includes/watersim-test-common.h"
#include "SimConfig.h"

int main() {
	SimConfig cfg;
	double fx, fy, fz;
	fx = 10;
	fy = 20;
	fz = 30;
	double tx, ty, tz;
	tx = 40;
	ty = 50;
	tz = 60;
	cfg.setFluidRegion(tx, ty, tz, fx, fy, fz);
	double fx2, fy2, fz2, tx2, ty2, tz2;
	cfg.getFluidRegion(fx2, fy2, fz2, tx2, ty2, tz2);
	assert(fx2 == fx);
	assert(fy2 == fy);
	assert(fz2 == fz);
	assert(tx2 == tx);
	assert(ty2 == ty);
	assert(tz2 == tz);
}
