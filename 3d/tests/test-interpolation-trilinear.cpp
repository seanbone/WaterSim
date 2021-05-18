/*
 * A test to check the trilinear interpolation routine of Mac3d.
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

	double min_x = 101;
	double max_x = 202;
	double min_y = 303;
	double max_y = 504;
	double min_z = 605;
	double max_z = 906;

	double r_size_x = 1/(max_x - min_x);
	double r_size_y = 1/(max_y - min_y);
	double r_size_z = 1/(max_z - min_z);

	posx = min_x + posx * (max_x - min_x);
	posy = min_y + posy * (max_y - min_y);
	posz = min_z + posz * (max_z - min_z);

	double v00 = Mac3d::linear_interpolation(v000, v100, min_x, r_size_x, posx);
	double v10 = Mac3d::linear_interpolation(v010, v110, min_x, r_size_x, posx);
	double v01 = Mac3d::linear_interpolation(v001, v101, min_x, r_size_x, posx);
	double v11 = Mac3d::linear_interpolation(v011, v111, min_x, r_size_x, posx);

	double v0 = Mac3d::linear_interpolation(v00, v10, min_y, r_size_y, posy);
	double v1 = Mac3d::linear_interpolation(v01, v11, min_y, r_size_y, posy);

	double expected = Mac3d::linear_interpolation(v0, v1, min_z, r_size_z, posz);

	double returned = Mac3d::trilinear_interpolation(v000, v001, v010, v011,
												  v100, v101, v110, v111,
												  min_x, min_y, min_z,
												  r_size_x, r_size_y, r_size_z,
												  posx, posy, posz);

	assert(std::abs(expected - returned) < interpolation_tolerance);
}

int main() {
	double v = 0.65;
	double expected = v;

	double min_x = 100;
	double max_x = 200;
	double min_y = 300;
	double max_y = 500;
	double min_z = 600;
	double max_z = 900;

	double r_size_x = 1/(max_x - min_x);
	double r_size_y = 1/(max_y - min_y);
	double r_size_z = 1/(max_z - min_z);

	double pos_x = min_x + v * (max_x - min_x);
	double pos_y = min_y + v * (max_y - min_y);
	double pos_z = min_z + v * (max_z - min_z);

	// Test X axis alone
	{
		double returned = Mac3d::trilinear_interpolation(0, 0, 0, 0, 1, 1, 1, 1,
												   min_x, min_y, min_z,
												   r_size_x, r_size_y, r_size_z,
												   pos_x, pos_y, pos_z);
		assert("Checking X-axis only" && std::abs(expected - returned) < interpolation_tolerance);
	}
	// Test Y axis alone
	{
		double returned = Mac3d::trilinear_interpolation(0, 0, 1, 1, 0, 0, 1, 1,
												   min_x, min_y, min_z,
												   r_size_x, r_size_y, r_size_z,
												   pos_x, pos_y, pos_z);
		assert("Checking Y-axis only" && std::abs(expected - returned) < interpolation_tolerance);
	}
	// Test Z axis alone
	{
		double returned = Mac3d::trilinear_interpolation(0, 1, 0, 1, 0, 1, 0, 1,
												   min_x, min_y, min_z,
												   r_size_x, r_size_y, r_size_z,
												   pos_x, pos_y, pos_z);
		assert("Checking Z-axis only" && std::abs(expected - returned) < interpolation_tolerance);
	}

	test_all_axes();
}
