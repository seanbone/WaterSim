/*
 * A test to check the normalized trilinear interpolation routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

void test_all_axes() {
	double v000 = 1;
	double v001 = 2;
	double v010 = 3;
	double v011 = 4;
	double v100 = 5;
	double v101 = 6;
	double v110 = 7;
	double v111 = 8;

	double posx = 0.45;
	double posy = 0.98;
	double posz = 0.13;

	double v00 = Mac3d::linear_interpolation_normalized(v000, v100, posx);
	double v10 = Mac3d::linear_interpolation_normalized(v010, v110, posx);
	double v01 = Mac3d::linear_interpolation_normalized(v001, v101, posx);
	double v11 = Mac3d::linear_interpolation_normalized(v011, v111, posx);

	double v0 = Mac3d::linear_interpolation_normalized(v00, v10, posy);
	double v1 = Mac3d::linear_interpolation_normalized(v01, v11, posy);

	double expected = Mac3d::linear_interpolation_normalized(v0, v1, posz);

	double returned = Mac3d::trilinear_interpolation_normalized(v000, v001, v010, v011,
	                                                            v100, v101, v110, v111,
	                                                            posx, posy, posz);

	assert(std::abs(expected - returned) < interpolation_tolerance);
}

int main() {
	double v = 0.65;
	double expected = v;

	// Test X axis alone
	{
		double returned = Mac3d::trilinear_interpolation_normalized(0, 0, 0, 0, 1, 1, 1, 1, v, v, v);
		assert("Checking X-axis only" && std::abs(expected - returned) < interpolation_tolerance);
	}
	// Test Y axis alone
	{
		double returned = Mac3d::trilinear_interpolation_normalized(0, 0, 1, 1, 0, 0, 1, 1, v, v, v);
		assert("Checking Y-axis only" && std::abs(expected - returned) < interpolation_tolerance);
	}
	// Test Z axis alone
	{
		double returned = Mac3d::trilinear_interpolation_normalized(0, 1, 0, 1, 0, 1, 0, 1, v, v, v);
		assert("Checking Z-axis only" && std::abs(expected - returned) < interpolation_tolerance);
	}

	test_all_axes();
}
