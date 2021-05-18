/*
 * A test to check the grid_interpolate routine of Mac3d.
 */
#include "includes/watersim-test-common.h"
#include "Mac3d.h"

int main() {
	unsigned nx = 5;
	unsigned ny = 5;
	unsigned nz = 5;
	double sx = 10;
	double sy = 10;
	double sz = 10;
	Mac3d mac(nx, ny, nz, sx, sy, sz);

	for (unsigned k = 0; k < nz; k++) {
		for (unsigned j = 0; j < ny; j++) {
			for (unsigned i = 0; i < nx; i++) {
				mac.ppressure_[i + nx * j + nx * ny * k] = i + nx * j + nx * ny * k;
			}
		}
	}

	// Actual tests
	/*
	{
		double norm_px = 0;
		double norm_py = 0;
		double norm_pz = 1.5;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = 37.5;

		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}
	 */
	{
		// Inner domain
		double norm_px = 1.25;
		double norm_py = 3.3;
		double norm_pz = 2.95;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		unsigned cell_x = px * mac.rcell_sizex_;
		unsigned cell_y = py * mac.rcell_sizey_;
		unsigned cell_z = pz * mac.rcell_sizez_;

		double* p = mac.ppressure_;

		unsigned i000 = (cell_x    ) + nx * (cell_y    ) + ny * nx * (cell_z    );
		unsigned i001 = (cell_x    ) + nx * (cell_y    ) + ny * nx * (cell_z + 1);
		unsigned i010 = (cell_x    ) + nx * (cell_y + 1) + ny * nx * (cell_z    );
		unsigned i011 = (cell_x    ) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);
		unsigned i100 = (cell_x + 1) + nx * (cell_y    ) + ny * nx * (cell_z    );
		unsigned i101 = (cell_x + 1) + nx * (cell_y    ) + ny * nx * (cell_z + 1);
		unsigned i110 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z    );
		unsigned i111 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);

		/*
		std::cout << "here:" << std::endl;
		std::cout << p[i000] << std::endl;
		std::cout << p[i100] << std::endl;
		std::cout << p[i010] << std::endl;
		std::cout << p[i110] << std::endl;
		std::cout << p[i001] << std::endl;
		std::cout << p[i101] << std::endl;
		std::cout << p[i011] << std::endl;
		std::cout << p[i111] << std::endl;

		std::cout << "**********" << std::endl;
		 */

		double expected = Mac3d::trilinear_interpolation_normalized(p[i000], p[i001], p[i010], p[i011],
															        p[i100], p[i101], p[i110], p[i111],
															        norm_px - cell_x, norm_py - cell_y, norm_pz - cell_z);

		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		//std::cout << "***\nexpected: " << expected << std::endl;
		//std::cout << "returned: " << returned << std::endl;

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}

	{
		// Bilinear on left face
		double norm_px = -0.5;
		double norm_py = 0.2;
		double norm_pz = 0.7;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = Mac3d::bilinear_interpolation_normalized(0, 25, 5, 30, norm_pz, norm_py);
		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}

	{
		// Bilinear on bottom face
		double norm_px = 1.5;
		double norm_py = -0.2;
		double norm_pz = 1.7;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = Mac3d::bilinear_interpolation_normalized(26, 27, 51, 52, norm_px-1, norm_pz-1);
		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}

	{
		// Bilinear on front face
		double norm_px = 1.5;
		double norm_py = 1.2;
		double norm_pz = -1.7;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = Mac3d::bilinear_interpolation_normalized(6, 11, 7, 12, norm_px-1, norm_py-1);
		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}

	/*
	{
		// Bilinear on back face
		double norm_px = 1.5;
		double norm_py = 1.2;
		double norm_pz = 5.7;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = Mac3d::bilinear_interpolation_normalized(0, 11, 7, 12, norm_px-1, norm_py-1);
		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}
	 */




	{
		// Linear on bottom-left edge
		double norm_px = -0.5;
		double norm_py = -1.5;
		double norm_pz = 1.5;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = Mac3d::linear_interpolation_normalized(0, 25, norm_pz); //Mac3d::bilinear_interpolation_normalized(0, 25, 5, 30, norm_pz, norm_py);
		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}

	{
		// Linear on front-left edge
		double norm_px = -1;
		double norm_py = 1.25;
		double norm_pz = -1;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = Mac3d::linear_interpolation_normalized(5, 10, norm_py-1);

		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}

	{
		// Linear on front-bottom edge
		double norm_px = 1.33;
		double norm_py = -0.5;
		double norm_pz = -0.1;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		double expected = norm_px;

		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}
	return 0;

	{
		double norm_px = -0.5;
		double norm_py = -0.5;
		double norm_pz = 0.5;
		double px = norm_px * mac.cell_sizex_;
		double py = norm_py * mac.cell_sizey_;
		double pz = norm_pz * mac.cell_sizez_;

		int cell_x = 0; //px * mac.rcell_sizex_;
		int cell_y = 0; //py * mac.rcell_sizey_;
		int cell_z = 0; //pz * mac.rcell_sizez_;

		double* p = mac.ppressure_;

		unsigned i000 = (cell_x    ) + nx * (cell_y    ) + ny * nx * (cell_z    );
		unsigned i001 = (cell_x    ) + nx * (cell_y    ) + ny * nx * (cell_z + 1);
		unsigned i010 = (cell_x    ) + nx * (cell_y + 1) + ny * nx * (cell_z    );
		unsigned i011 = (cell_x    ) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);
		unsigned i100 = (cell_x + 1) + nx * (cell_y    ) + ny * nx * (cell_z    );
		unsigned i101 = (cell_x + 1) + nx * (cell_y    ) + ny * nx * (cell_z + 1);
		unsigned i110 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z    );
		unsigned i111 = (cell_x + 1) + nx * (cell_y + 1) + ny * nx * (cell_z + 1);

		std::cout << "here:" << std::endl;
		std::cout << p[i000] << std::endl;
		std::cout << p[i100] << std::endl;
		std::cout << p[i010] << std::endl;
		std::cout << p[i110] << std::endl;
		std::cout << p[i001] << std::endl;
		std::cout << p[i101] << std::endl;
		std::cout << p[i011] << std::endl;
		std::cout << p[i111] << std::endl;

		std::cout << "**********" << std::endl;

		//double expected = Mac3d::trilinear_interpolation_normalized(p[i000], p[i001], p[i010], p[i011],
		//                                                            p[i100], p[i101], p[i110], p[i111],
		//                                                            norm_px - cell_x, norm_py - cell_y, norm_pz - cell_z);
		double expected = 15; // = Mac3d::bilinear_interpolation_normalized(0, 5, 25, 30, 0.5, 0.5);

		double returned = mac.grid_interpolate<Mac3d::GRID_P>(px, py, pz);

		std::cout << "***\nexpected: " << expected << std::endl;
		std::cout << "returned: " << returned << std::endl;

		assert(std::abs(expected - returned) < interpolation_tolerance);
	}
}
