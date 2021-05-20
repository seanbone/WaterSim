/*
 * A test to check the normalized linear interpolation routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

int main() {
	double v0 = 10;
	double v1 = 20;
	double pos = 0.75;
	double expected = 17.5;

	double returned = Mac3d::linear_interpolation_normalized(v0, v1, pos);

	assert(std::abs(expected - returned) < interpolation_tolerance);
}
