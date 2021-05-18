/*
 * A test to check the linear interpolation routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

int main() {
	double v0 = 10;
	double v1 = 20;
	double minX = 100;
	double maxX = 200;
	double pos = 175;
	double expected = 17.5;

	double returned = Mac3d::linear_interpolation(v0, v1, minX, 1./(maxX-minX), pos);

	assert(std::abs(expected - returned) < interpolation_tolerance);
}
