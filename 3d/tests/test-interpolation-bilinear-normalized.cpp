/*
 * A test to check the normalized bilinear interpolation routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

int main() {
	double v00 = 1;
	double v01 = 2;
	double v10 = 3;
	double v11 = 4;

	double posx = 0.55;
	double posy = 0.65;

	double v0 = Mac3d::linear_interpolation_normalized(v00, v10, posx);
	double v1 = Mac3d::linear_interpolation_normalized(v01, v11, posx);

	double expected = Mac3d::linear_interpolation_normalized(v0, v1, posy);

	double returned = Mac3d::bilinear_interpolation_normalized(v00, v01, v10, v11, posx, posy);

	assert(std::abs(expected - returned) < interpolation_tolerance);
}
