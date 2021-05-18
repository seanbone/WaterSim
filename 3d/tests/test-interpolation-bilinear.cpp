/*
 * A test to check the bilinear interpolation routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

int main() {
	double v00 = 1;
	double v01 = 2;
	double v10 = 3;
	double v11 = 4;

	double min_x = 100;
	double max_x = 200;
	double min_y = 300;
	double max_y = 500;

	double r_size_x = 1/(max_x - min_x);
	double r_size_y = 1/(max_y - min_y);

	double posx = 155;
	double posy = 422;

	double v0 = Mac3d::linear_interpolation(v00, v10, min_x, r_size_x, posx);
	double v1 = Mac3d::linear_interpolation(v01, v11, min_x, r_size_x, posx);

	double expected = Mac3d::linear_interpolation(v0, v1, min_y, r_size_y, posy);

	double returned = Mac3d::bilinear_interpolation(v00, v01, v10, v11,
												    min_x, min_y,
												    r_size_x, r_size_y,
												    posx, posy);

	assert(std::abs(expected - returned) < interpolation_tolerance);
}
